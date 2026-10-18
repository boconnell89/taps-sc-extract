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
#include <stdio.h>
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

/*
 * One read to BAQ. ref points at reference bases [ref_lo, ref_hi) of the
 * read's contig (uppercase, unmasked); ref_hi is the contig length when the
 * slice reaches the contig end. The slice must cover taps_baq_need(b), clipped
 * to the contig. mask is the base the HMM sees as N ('C' for OT reads, 'G'
 * for OB reads in the TAPS modes, 0 for none), so no masked copy of the
 * reference is ever stored.
 */
typedef struct {
    bam1_t *b;
    const char *ref;
    int64_t ref_lo;
    int64_t ref_hi;
    int mask;
} taps_baq_job;

/*
 * Reference span sam_prob_realn can read for b: [pos - slack, end + slack)
 * with slack = 2*l_qseq + ref span + 16 (realn.c widens the window by at most
 * l_qseq plus half a band, and the band is at most ref span + l_qseq + 3).
 * Not clipped to the contig. Must match `baq_need` in extract.rs.
 */
static void taps_baq_need(const bam1_t *b, int64_t *lo, int64_t *hi)
{
    const uint32_t *cigar = bam_get_cigar(b);
    int64_t span = 0;
    for (uint32_t k = 0; k < b->core.n_cigar; ++k) {
        int op = cigar[k] & 0xf;
        if (op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF || op == BAM_CDEL || op == BAM_CREF_SKIP)
            span += cigar[k] >> 4;
    }
    int64_t slack = 2 * (int64_t)b->core.l_qseq + span + 16;
    *lo = b->core.pos - slack < 0 ? 0 : b->core.pos - slack;
    *hi = b->core.pos + span + slack;
}

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
    char *fb; size_t fb_cap;      /* masked reference span for htslib_one() */
    int gpu_np;                   /* reads prepared for the device (taps_baq_gpu_*) */
} taps_baq_ws;

void *taps_baq_ws_new(void)
{
    return calloc(1, sizeof(taps_baq_ws));
}

void taps_baq_ws_free(void *p)
{
    taps_baq_ws *ws = p;
    if (!ws) return;
    free(ws->arena); free(ws->prep); free(ws->order); free(ws->mat); free(ws->fin); free(ws->fb);
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

/* Per-read lookup tables (filled once at load): SEQ byte -> two HMM base
 * codes, and reference byte -> HMM code with the TAPS mask (none, C, G)
 * already applied, as seq_nt16_int[seq_nt16_table[ch]]. */
static uint8_t taps_seq_pair[256][2];
static uint8_t taps_ref_code[3][256];

__attribute__((constructor)) static void taps_tables_init(void)
{
    for (int v = 0; v < 256; ++v) {
        taps_seq_pair[v][0] = seq_nt16_int[v >> 4];
        taps_seq_pair[v][1] = seq_nt16_int[v & 0xf];
        int code = seq_nt16_int[seq_nt16_table[v]];
        int n = seq_nt16_int[seq_nt16_table['N']];
        taps_ref_code[0][v] = code;
        taps_ref_code[1][v] = v == 'C' ? n : code;
        taps_ref_code[2][v] = v == 'G' ? n : code;
    }
}

/* bam_aux_get(b, "BQ") || bam_aux_get(b, "ZQ") in one pass over the aux
 * data. Malformed aux data counts as present (htslib then handles the read). */
static int has_bq_or_zq(const bam1_t *b)
{
    const uint8_t *s = bam_get_aux(b), *end = b->data + b->l_data;
    while (end - s >= 3) {
        if ((s[0] == 'B' || s[0] == 'Z') && s[1] == 'Q') return 1;
        uint8_t t = s[2];
        s += 3;
        switch (t) {
        case 'A': case 'c': case 'C': s += 1; break;
        case 's': case 'S': s += 2; break;
        case 'i': case 'I': case 'f': s += 4; break;
        case 'd': s += 8; break;
        case 'Z': case 'H': {
            const uint8_t *z = memchr(s, 0, end - s);
            if (!z) return 1;
            s = z + 1;
            break;
        }
        case 'B': {
            if (end - s < 5) return 1;
            int sz;
            switch (s[0]) {
            case 'c': case 'C': sz = 1; break;
            case 's': case 'S': sz = 2; break;
            case 'i': case 'I': case 'f': sz = 4; break;
            default: return 1;
            }
            uint32_t cnt;
            memcpy(&cnt, s + 1, 4);
            if ((uint64_t)cnt * sz > (uint64_t)(end - s - 5)) return 1;
            s += 5 + (size_t)cnt * sz;
            break;
        }
        default: return 1;
        }
    }
    return s != end; /* trailing bytes: malformed */
}

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
    if (has_bq_or_zq(b)) return 2;

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
    for (i = 0; i + 1 < (hts_pos_t)lq; i += 2)
        memcpy(tseq + i, taps_seq_pair[seq[i >> 1]], 2);
    if (lq & 1) tseq[lq - 1] = taps_seq_pair[seq[lq >> 1]][0];
    if (xb < j->ref_lo) return 2; /* slice does not cover it (should not happen) */
    /* TAPS mask applied here instead of stored (taps_ref_code). */
    const uint8_t *rc = taps_ref_code[j->mask == 'C' ? 1 : j->mask == 'G' ? 2 : 0];
    if (j->mask && j->mask != 'C' && j->mask != 'G') return 2;
    if (xe > j->ref_hi) xe = j->ref_hi;
    const unsigned char *rs = (const unsigned char *)j->ref - j->ref_lo;
    for (i = xb; i < xe; ++i) {
        unsigned char ch = rs[i];
        if (ch == '\0') { xe = i; break; }
        tref[i-xb] = rc[ch];
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
/* state is int (CPU kernels) or, when state16 is set, the device's int16. */
static inline __attribute__((always_inline)) void finish_t(bam1_t *b, hts_pos_t xb0, const int *state32,
                                                           const int16_t *state16, const uint8_t *q,
                                                           uint8_t *bq, uint8_t *left, uint8_t *rght)
{
#define STATE(i) (state16 ? (int)state16[i] : state32[i])
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
                bq[i] = ((STATE(i)&3) != 0 || STATE(i)>>2 != x - xb + (i - y))? 0 : q[i];
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
#undef STATE
}

static void finish(bam1_t *b, hts_pos_t xb0, const int *state, const uint8_t *q,
                   uint8_t *bq, uint8_t *left, uint8_t *rght)
{
    finish_t(b, xb0, state, NULL, q, bq, left, rght);
}

static void finish16(bam1_t *b, hts_pos_t xb0, const int16_t *state, const uint8_t *q,
                     uint8_t *bq, uint8_t *left, uint8_t *rght)
{
    finish_t(b, xb0, NULL, state, q, bq, left, rght);
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
 * htslib's sam_prob_realn on one read. It indexes the reference by contig
 * coordinate, so it gets a masked copy of only the span it can read
 * (taps_baq_need), addressed as if the copy started at coordinate 0.
 */
static void htslib_one(taps_baq_ws *ws, const taps_baq_job *j)
{
    int64_t lo, hi;
    taps_baq_need(j->b, &lo, &hi);
    if (lo < j->ref_lo) lo = j->ref_lo;
    if (hi > j->ref_hi) hi = j->ref_hi;
    size_t n = hi > lo ? (size_t)(hi - lo) : 0;
    if (grow((void **)&ws->fb, &ws->fb_cap, n + 1, 0) < 0) {
        fprintf(stderr, "taps BAQ: out of memory\n");
        abort();
    }
    const char *src = j->ref + (lo - j->ref_lo);
    for (size_t i = 0; i < n; ++i)
        ws->fb[i] = (j->mask && src[i] == j->mask) ? 'N' : src[i];
    ws->fb[n] = '\0';
    /* Reads stay within [lo, hi); sam_prob_realn stops at ref_len = hi, which
     * is the contig end whenever the read's span reaches it. */
    const char *base = (const char *)((uintptr_t)ws->fb - (uintptr_t)lo);
    sam_prob_realn(j->b, base, hi, BAQ_APPLY | BAQ_EXTEND);
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
            htslib_one(ws, &jobs[j]);
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
            htslib_one(ws, &jobs[j]);
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
                htslib_one(ws, &jobs[ws->prep[ws->order[t]].job]);
            }
        }
        g = h;
    }
    return simd;
}

/* ---- GPU path (csrc/gpu/): prepare and finish here, HMM on the device ---- */

#include "gpu/taps_gpu.h"

typedef struct {
    int64_t n_dev;     /* reads for the device */
    int64_t ref_bytes; /* sum of their l_ref */
    int64_t bases;     /* sum of their l_query */
    int64_t n_cpu;     /* reads too long or too wide for the device kernel */
} taps_gpu_sizes;

/*
 * Step 1: everything sam_prob_realn does before the HMM, for jobs[0..n).
 * Nothing is modified yet. Reads the device does not take are listed in
 * cpu_jobs (job indices; caller array of n) for the CPU path
 * (taps_baq_batch), which must run after taps_baq_gpu_finish (they share the
 * workspace): reads htslib must handle, and reads beyond the device kernel's
 * limits (l_query > TAPS_GPU_MAXLQ or i_dim > TAPS_GPU_MAXIDIM). If the device
 * fails, the caller can still run every job on the CPU.
 */
void taps_baq_gpu_prepare(void *wsp, const taps_baq_job *jobs, int n, taps_gpu_sizes *sz, int *cpu_jobs)
{
    taps_baq_ws *ws = wsp;
    memset(sz, 0, sizeof *sz);
    if (grow((void **)&ws->prep, &ws->prep_cap, (size_t)n * sizeof(prep_t), 0) < 0) {
        for (int j = 0; j < n; ++j) cpu_jobs[sz->n_cpu++] = j;
        ws->gpu_np = 0;
        return;
    }
    size_t used = 0;
    int np = 0;
    for (int j = 0; j < n; ++j) {
        int r = prepare(ws, &used, &jobs[j], j, &ws->prep[np]);
        if (r == 0) {
            const prep_t *p = &ws->prep[np];
            int bw2 = p->bw_key * 2 + 1;
            int i_dim = bw2 < p->l_ref ? bw2 * 3 + 6 : p->l_ref * 3 + 6;
            if (p->l_query > TAPS_GPU_MAXLQ || i_dim > TAPS_GPU_MAXIDIM) {
                cpu_jobs[sz->n_cpu++] = j;
                continue;
            }
            sz->ref_bytes += p->l_ref;
            sz->bases += p->l_query;
            ++np;
        } else if (r != 1) {
            cpu_jobs[sz->n_cpu++] = j;
        }
    }
    ws->gpu_np = np;
    sz->n_dev = np;
}

/*
 * Step 2: pack the prepared reads for the device. Offsets start at 0; the
 * HMM parameters are computed here with probaln.c's float expressions.
 */
void taps_baq_gpu_export(void *wsp, const taps_baq_job *jobs, taps_gpu_read *rd, uint8_t *refb,
                         uint8_t *qry, uint8_t *qual)
{
    taps_baq_ws *ws = wsp;
    const float cd = 0.001f, ce = 0.1f; /* probaln_par_t {d, e} for reads <= 1000 bp */
    size_t ro = 0, bo = 0;
    for (int t = 0; t < ws->gpu_np; ++t) {
        const prep_t *p = &ws->prep[t];
        taps_gpu_read *g = &rd[t];
        int l_query = p->l_query, l_ref = p->l_ref;
        double sM, sI;
        memset(g, 0, sizeof *g);
        sM = sI = 1. / (2 * l_query + 2);
        g->m[0*3+0] = (1 - cd - cd) * (1 - sM);
        g->m[0*3+1] = g->m[0*3+2] = cd * (1 - sM);
        g->m[1*3+0] = (1 - ce) * (1 - sI);
        g->m[1*3+1] = ce * (1 - sI);
        g->m[1*3+2] = 0.;
        g->m[2*3+0] = 1 - ce;
        g->m[2*3+1] = 0.;
        g->m[2*3+2] = ce;
        g->bM = (1 - cd) / l_ref;
        g->bI = cd / l_ref;
        g->sM = sM;
        g->l_ref = l_ref;
        g->l_query = l_query;
        g->bw = p->bw_key;
        {
            int bw2 = p->bw_key * 2 + 1;
            g->i_dim = bw2 < l_ref ? bw2 * 3 + 6 : l_ref * 3 + 6;
        }
        g->ref_off = (unsigned)ro;
        g->base_off = (unsigned)bo;
        memcpy(refb + ro, ws->arena + p->tref, l_ref);
        memcpy(qry + bo, ws->arena + p->tseq, l_query);
        memcpy(qual + bo, bam_get_qual(jobs[p->job].b), l_query);
        ro += l_ref;
        bo += l_query;
    }
}

/* q from a posterior maximum, exactly as probaln_glocal and the CPU kernels. */
static uint8_t q_of_max(double mx)
{
    int k = (int)(-4.343 * log(1. - mx) + .499);
    return k > 100 ? 99 : k;
}

/*
 * Step 3: the device's state[] and q[] (q == 0xFF: recompute from maxv with
 * glibc) for each read, then sam_prob_realn's BAQ_EXTEND and QUAL cap.
 */
void taps_baq_gpu_finish(void *wsp, const taps_baq_job *jobs, const int16_t *state16, const uint8_t *q,
                         const double *maxv)
{
    taps_baq_ws *ws = wsp;
    size_t bo = 0;
    for (int t = 0; t < ws->gpu_np; ++t) {
        const prep_t *p = &ws->prep[t];
        int lq = p->l_query;
        /* q == 0xFF (rare): recomputed from maxv into the arena copy. */
        const uint8_t *qq = q + bo;
        if (maxv && memchr(qq, 0xFF, lq)) {
            uint8_t *fx = ws->arena + p->q;
            for (int i = 0; i < lq; ++i) fx[i] = qq[i] != 0xFF ? qq[i] : q_of_max(maxv[bo + i]);
            qq = fx;
        }
        if (grow((void **)&ws->fin, &ws->fin_cap, 3 * al16((size_t)lq + 1), 0) < 0) {
            fprintf(stderr, "taps BAQ: out of memory\n");
            abort();
        }
        uint8_t *bq = ws->fin, *left = bq + al16(lq + 1), *rght = left + al16(lq + 1);
        finish16(jobs[p->job].b, p->xb, state16 + bo, qq, bq, left, rght);
        bo += lq;
    }
}

/* Fill the qual -> probability table exactly as probaln.c (float). */
void taps_baq_gpu_q2p(float *out)
{
    for (int i = 0; i < 256; ++i) out[i] = taps_qual2prob[i];
}
