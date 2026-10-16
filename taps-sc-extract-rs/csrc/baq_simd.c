/*
 * Batched SIMD BAQ: the same result as htslib 1.19.1
 * sam_prob_realn(b, ref, ref_len, BAQ_APPLY|BAQ_EXTEND) for each read, with the
 * profile HMM (probaln_glocal) run on 4 (AVX2) or 8 (AVX-512) reads at once.
 *
 * taps_baq_batch() takes a batch of reads. For each read it does what
 * sam_prob_realn does before the HMM (alignment span, band, translated
 * sequences), groups reads whose HMM matrices have the same shape
 * (l_query, l_ref, band), runs the lane-parallel HMM on each group, then does
 * what sam_prob_realn does after it (BAQ_EXTEND running min/max, cap QUAL).
 * The ZQ tag is not appended: nothing downstream reads it.
 *
 * Reads the fast path does not cover go to htslib's sam_prob_realn unchanged:
 * an existing BQ/ZQ tag, reads over 1000 bp (long-read HMM parameters), other
 * flags than BAQ_APPLY|BAQ_EXTEND, or an empty reference span.
 *
 * Wrapper code follows htslib realn.c: Copyright (C) 2009-2011, 2014-2016,
 * 2018-2021, 2023 Genome Research Ltd. Portions copyright (C) 2009-2011
 * Broad Institute. MIT licence.
 */

#define _GNU_SOURCE /* qsort_r */
#include <errno.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include <htslib/hts.h>
#include <htslib/sam.h>

#define EI .25
#define EM .33333333333
#define set_u(u, b, i, k) { int x=(i)-(b); x=x>0?x:0; (u)=((k)-x+1)*3; }

typedef struct {
    const uint8_t *ref;   /* translated reference window, codes 0-4 */
    const uint8_t *query; /* translated read, codes 0-4 */
    const uint8_t *iqual; /* base qualities */
} taps_lane_in;

static float taps_qual2prob[256];

/* Same table as probaln.c's g_qual2prob. */
__attribute__((constructor)) static void taps_qual2prob_init(void)
{
    for (int i = 0; i < 256; ++i) taps_qual2prob[i] = pow(10, -i/10.);
}

/* ---- kernels (x86 only; elsewhere every read goes to sam_prob_realn) ---- */

#if (defined(__x86_64__) || defined(__i386__)) && !defined(TAPS_BAQ_NO_SIMD)
#define TAPS_BAQ_SIMD 1

#pragma GCC push_options
#pragma GCC target("avx2")
#define VW 4
#define KFN taps_probaln_x4_avx2
#include "baq_kernel.h"
#undef KFN
#undef VW
/* 8 lanes on AVX2 instructions: lets an AVX2-only host test the 8-lane logic. */
#define VW 8
#define KFN taps_probaln_x8_avx2
#include "baq_kernel.h"
#undef KFN
#undef VW
#pragma GCC pop_options

#pragma GCC push_options
#pragma GCC target("avx512f")
#define VW 8
#define KFN taps_probaln_x8_avx512
#include "baq_kernel.h"
#undef KFN
#undef VW
#pragma GCC pop_options
#endif

/* ---- widths ---- */

enum {
    TAPS_BAQ_HTSLIB = 0,     /* htslib sam_prob_realn per read */
    TAPS_BAQ_AVX2 = 4,       /* 4 lanes, AVX2 */
    TAPS_BAQ_AVX512 = 8,     /* 8 lanes, AVX-512F */
    TAPS_BAQ_AVX2_X8 = 108,  /* 8 lanes on AVX2 (testing) */
};

/* Widest kernel this CPU runs: 8, 4, or 0. */
int taps_baq_cpu_width(void)
{
#ifndef TAPS_BAQ_SIMD
    return TAPS_BAQ_HTSLIB;
#endif
    __builtin_cpu_init();
    if (__builtin_cpu_supports("avx512f")) return TAPS_BAQ_AVX512;
    if (__builtin_cpu_supports("avx2")) return TAPS_BAQ_AVX2;
    return TAPS_BAQ_HTSLIB;
}

int taps_baq_width_ok(int width)
{
#ifndef TAPS_BAQ_SIMD
    return width == TAPS_BAQ_HTSLIB;
#endif
    __builtin_cpu_init();
    switch (width) {
    case TAPS_BAQ_HTSLIB: return 1;
    case TAPS_BAQ_AVX2: case TAPS_BAQ_AVX2_X8: return __builtin_cpu_supports("avx2");
    case TAPS_BAQ_AVX512: return __builtin_cpu_supports("avx512f");
    default: return 0;
    }
}

/* ---- workspace ---- */

typedef struct {
    bam1_t *b;
    const char *ref;
    hts_pos_t ref_len;
} taps_baq_job;

typedef struct {
    int job;            /* index into jobs[] */
    int l_query, l_ref, bw_key;
    hts_pos_t xb;       /* reference start of the HMM window */
    size_t tseq, tref, state, q;  /* arena offsets */
} prep_t;

typedef struct {
    uint8_t *arena; size_t arena_cap;
    prep_t *prep; size_t prep_cap;
    int *order; size_t order_cap;
    void *mat; size_t mat_cap;   /* f, b, s, R, Q, QL for one group */
    uint8_t *fin; size_t fin_cap; /* bq, left, rght for finish() */
} taps_baq_ws;

void *taps_baq_ws_new(void)
{
    return calloc(1, sizeof(taps_baq_ws));
}

void taps_baq_ws_free(void *p)
{
    taps_baq_ws *ws = p;
    if (!ws) return;
    free(ws->arena); free(ws->prep); free(ws->order); free(ws->mat); free(ws->fin);
    free(ws);
}

static int grow(void **p, size_t *cap, size_t want, size_t align)
{
    if (*cap >= want) return 0;
    size_t n = *cap ? *cap : 4096;
    while (n < want) n *= 2;
    void *q;
    if (align) {
        q = aligned_alloc(align, (n + align - 1) / align * align);
        if (!q) return -1;
        free(*p);
    } else {
        q = realloc(*p, n);
        if (!q) return -1;
    }
    *p = q; *cap = n;
    return 0;
}

static size_t al16(size_t n) { return (n + 15) & ~(size_t)15; }

/*
 * Before the HMM: sam_prob_realn up to the probaln_glocal call.
 * Returns 0 queued, 1 no-op (htslib returns without changes), 2 use htslib,
 * -1 out of memory.
 */
static int prepare(taps_baq_ws *ws, size_t *used, const taps_baq_job *j, int job, prep_t *p)
{
    bam1_t *b = j->b;
    bam1_core_t *c = &b->core;
    uint32_t *cigar = bam_get_cigar(b);
    uint8_t *qual = bam_get_qual(b);
    int k, bw, y, yb, ye, xb, xe;
    hts_pos_t i, x;

    if (c->l_qseq > 1000) return 2;
    if ((c->flag & BAM_FUNMAP) || c->l_qseq == 0 || qual[0] == (uint8_t)-1) return 1;
    if (bam_aux_get(b, "BQ") || bam_aux_get(b, "ZQ")) return 2;

    x = c->pos, y = 0, yb = ye = xb = xe = -1;
    for (k = 0; k < c->n_cigar; ++k) {
        int op, l;
        op = cigar[k]&0xf; l = cigar[k]>>4;
        if (op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF) {
            if (yb < 0) yb = y;
            if (xb < 0) xb = x;
            ye = y + l; xe = x + l;
            x += l; y += l;
        } else if (op == BAM_CSOFT_CLIP || op == BAM_CINS) y += l;
        else if (op == BAM_CDEL) x += l;
        else if (op == BAM_CREF_SKIP) return 1;
    }
    if (xb == -1) return 1;
    bw = 7;
    if (abs((xe - xb) - (ye - yb)) > bw)
        bw = abs((xe - xb) - (ye - yb)) + 3;

    xb -= yb + bw/2; if (xb < 0) xb = 0;
    xe += c->l_qseq - ye + bw/2;
    if (xe - xb - c->l_qseq > bw)
        xb += (xe - xb - c->l_qseq - bw) / 2, xe -= (xe - xb - c->l_qseq - bw) / 2;

    size_t lq = c->l_qseq;
    size_t lref_cap = xe > xb ? (size_t)(xe - xb) : 1;
    size_t o_tseq = *used;
    size_t o_tref = o_tseq + al16(lq);
    size_t o_state = o_tref + al16(lref_cap);
    size_t o_q = o_state + al16(lq * sizeof(int));
    size_t end = o_q + al16(lq);
    if (grow((void **)&ws->arena, &ws->arena_cap, end, 0) < 0) return -1;
    uint8_t *tseq = ws->arena + o_tseq, *tref = ws->arena + o_tref;
    uint8_t *seq = bam_get_seq(b);
    for (i = 0; i < (hts_pos_t)lq; ++i)
        tseq[i] = seq_nt16_int[bam_seqi(seq, i)];
    for (i = xb; i < xe; ++i) {
        if (i >= j->ref_len || j->ref[i] == '\0') { xe = i; break; }
        tref[i-xb] = seq_nt16_int[seq_nt16_table[(unsigned char)j->ref[i]]];
    }
    /* probaln_glocal returns 0 without writing state/q here; leave to htslib. */
    if (xe - xb <= 0) return 2;

    int l_ref = xe - xb, l_query = (int)lq;
    int bwk = l_ref > l_query ? l_ref : l_query;
    if (bwk > bw) bwk = bw;
    if (bwk < abs(l_ref - l_query)) bwk = abs(l_ref - l_query);

    p->job = job;
    p->l_query = l_query;
    p->l_ref = l_ref;
    p->bw_key = bwk;
    p->xb = xb;
    p->tseq = o_tseq; p->tref = o_tref; p->state = o_state; p->q = o_q;
    *used = end;
    return 0;
}

/* After the HMM: sam_prob_realn's BAQ_EXTEND block and BAQ_APPLY, without ZQ. */
static void finish(bam1_t *b, hts_pos_t xb0, const int *state, const uint8_t *q,
                   uint8_t *bq, uint8_t *left, uint8_t *rght)
{
    bam1_core_t *c = &b->core;
    uint32_t *cigar = bam_get_cigar(b);
    uint8_t *qual = bam_get_qual(b);
    int k, y, len = 0;
    hts_pos_t i, x, xb = xb0;

    memcpy(bq, qual, c->l_qseq); bq[c->l_qseq] = 0;
    for (k = 0, x = c->pos, y = 0; k < c->n_cigar; ++k) {
        int op = cigar[k]&0xf, l = cigar[k]>>4;
        if (op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF) {
            if ((k + 1) < c->n_cigar) {
                int next_op = bam_cigar_op(cigar[k + 1]);
                if (next_op == BAM_CMATCH || next_op == BAM_CEQUAL || next_op == BAM_CDIFF) {
                    len += l;
                    continue;
                }
            }
            l += len;
            len = 0;
        }
        if (l == 0) continue;
        if (op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF) {
            if (l > c->l_qseq - y)
                l = c->l_qseq - y;
            for (i = y; i < y + l; ++i)
                bq[i] = ((state[i]&3) != 0 || state[i]>>2 != x - xb + (i - y))? 0 : q[i];
            for (left[y] = bq[y], i = y + 1; i < y + l; ++i)
                left[i] = bq[i] > left[i-1]? bq[i] : left[i-1];
            for (rght[y+l-1] = bq[y+l-1], i = y + l - 2; i >= y; --i)
                rght[i] = bq[i] > rght[i+1]? bq[i] : rght[i+1];
            for (i = y; i < y + l; ++i)
                bq[i] = left[i] < rght[i]? left[i] : rght[i];
            x += l; y += l;
        } else if (op == BAM_CSOFT_CLIP || op == BAM_CINS) {
            if (l > c->l_qseq - y)
                l = c->l_qseq - y;
            y += l;
        } else if (op == BAM_CDEL) {
            x += l;
        }
    }
    for (i = 0; i < c->l_qseq; ++i) bq[i] = 64 + (qual[i] <= bq[i]? 0 : qual[i] - bq[i]);
    for (i = 0; i < c->l_qseq; ++i) qual[i] -= bq[i] - 64;
}

static int cmp_prep(const void *a, const void *b, void *arg)
{
    const prep_t *preps = arg;
    const prep_t *x = &preps[*(const int *)a], *y = &preps[*(const int *)b];
    if (x->l_query != y->l_query) return x->l_query - y->l_query;
    if (x->l_ref != y->l_ref) return x->l_ref - y->l_ref;
    if (x->bw_key != y->bw_key) return x->bw_key - y->bw_key;
    return x->job - y->job;
}

static int run_group(taps_baq_ws *ws, const taps_baq_job *jobs, const int *idx, int n_real,
                     int width)
{
    int VWn = width == TAPS_BAQ_AVX2 ? 4 : 8;
    const prep_t *p0 = &ws->prep[idx[0]];
    int lq = p0->l_query, lr = p0->l_ref, bwk = p0->bw_key;
    int bw2 = bwk * 2 + 1;
    size_t i_dim = bw2 < lr ? (size_t)bw2 * 3 + 6 : (size_t)lr * 3 + 6;
    size_t vb = (size_t)VWn * 8;
    size_t fb = vb * (size_t)(lq + 1) * i_dim;
    size_t mat = 2 * fb + vb * ((size_t)(lq + 2) + lr + 2 * (size_t)lq);
    if (grow(&ws->mat, &ws->mat_cap, mat, 64) < 0) return -1;
    char *mp = ws->mat;
    void *f = mp, *bb = mp + fb, *s = mp + 2 * fb;
    void *R = (char *)s + vb * (lq + 2);
    void *Q = (char *)R + vb * lr;
    void *QL = (char *)Q + vb * lq;

    taps_lane_in in[8];
    int *st[8];
    uint8_t *qq[8];
    for (int l = 0; l < VWn; ++l) {
        /* A short group is padded with copies of its last read. */
        const prep_t *p = &ws->prep[idx[l < n_real ? l : n_real - 1]];
        in[l].query = ws->arena + p->tseq;
        in[l].ref = ws->arena + p->tref;
        in[l].iqual = bam_get_qual(jobs[p->job].b);
        st[l] = (int *)(ws->arena + p->state);
        qq[l] = ws->arena + p->q;
    }
    /* Illumina parameters (sam_prob_realn, l_qseq <= 1000). probaln_par_t
     * holds them as float, and htslib evaluates e.g. (1 - c->d) / l_ref in
     * float; doubles here shift HMM values in about the 7th digit. */
    const float cd = 0.001f, ce = 0.1f;
#ifdef TAPS_BAQ_SIMD
    switch (width) {
    case TAPS_BAQ_AVX2:
        taps_probaln_x4_avx2(in, lr, lq, bwk, cd, ce, f, bb, s, R, Q, QL, st, qq); break;
    case TAPS_BAQ_AVX2_X8:
        taps_probaln_x8_avx2(in, lr, lq, bwk, cd, ce, f, bb, s, R, Q, QL, st, qq); break;
    default:
        taps_probaln_x8_avx512(in, lr, lq, bwk, cd, ce, f, bb, s, R, Q, QL, st, qq); break;
    }
#else
    (void)cd; (void)ce; (void)f; (void)bb; (void)s; (void)R; (void)Q; (void)QL; (void)in; (void)st; (void)qq;
    return -1;
#endif
    if (grow((void **)&ws->fin, &ws->fin_cap, 3 * al16((size_t)lq + 1), 0) < 0) return -1;
    for (int l = 0; l < n_real; ++l) {
        const prep_t *p = &ws->prep[idx[l]];
        uint8_t *bq = ws->fin, *left = bq + al16(lq + 1), *rght = left + al16(lq + 1);
        finish(jobs[p->job].b, p->xb, (int *)(ws->arena + p->state), ws->arena + p->q,
               bq, left, rght);
    }
    return 0;
}

/*
 * BAQ (BAQ_APPLY|BAQ_EXTEND) for jobs[0..n). width: 0 htslib, 4 AVX2,
 * 8 AVX-512, 108 eight lanes on AVX2. Every read is processed exactly once.
 * Returns the number of reads that went through the SIMD kernel.
 */
int taps_baq_batch(void *wsp, const taps_baq_job *jobs, int n, int width)
{
    taps_baq_ws *ws = wsp;
    if (width == TAPS_BAQ_HTSLIB || !taps_baq_width_ok(width)
        || grow((void **)&ws->prep, &ws->prep_cap, (size_t)n * sizeof(prep_t), 0) < 0
        || grow((void **)&ws->order, &ws->order_cap, (size_t)n * sizeof(int), 0) < 0) {
        for (int j = 0; j < n; ++j)
            sam_prob_realn(jobs[j].b, jobs[j].ref, jobs[j].ref_len, BAQ_APPLY | BAQ_EXTEND);
        return 0;
    }
    size_t used = 0;
    int np = 0;
    for (int j = 0; j < n; ++j) {
        int r = prepare(ws, &used, &jobs[j], j, &ws->prep[np]);
        if (r == 0) {
            ws->order[np] = np;
            ++np;
        } else if (r != 1) {
            sam_prob_realn(jobs[j].b, jobs[j].ref, jobs[j].ref_len, BAQ_APPLY | BAQ_EXTEND);
        }
    }
    qsort_r(ws->order, np, sizeof(int), cmp_prep, ws->prep);

    int VWn = width == TAPS_BAQ_AVX2 ? 4 : 8, simd = 0;
    for (int g = 0; g < np; ) {
        const prep_t *p0 = &ws->prep[ws->order[g]];
        int h = g + 1;
        while (h < np && h - g < VWn) {
            const prep_t *ph = &ws->prep[ws->order[h]];
            if (ph->l_query != p0->l_query || ph->l_ref != p0->l_ref || ph->bw_key != p0->bw_key) break;
            ++h;
        }
        if (run_group(ws, jobs, &ws->order[g], h - g, width) == 0) {
            simd += h - g;
        } else {
            for (int t = g; t < h; ++t) {
                const taps_baq_job *j = &jobs[ws->prep[ws->order[t]].job];
                sam_prob_realn(j->b, j->ref, j->ref_len, BAQ_APPLY | BAQ_EXTEND);
            }
        }
        g = h;
    }
    return simd;
}
