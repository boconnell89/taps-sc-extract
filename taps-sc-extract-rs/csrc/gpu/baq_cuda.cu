/*
 * GPU BAQ backend for NVIDIA (CUDA) and, compiled with hipcc and TAPS_HIP,
 * AMD (HIP). The kernel is baq_hmm_body.h; see taps_gpu.h for the API.
 * Build flags keep the arithmetic exact: CUDA -fmad=false --prec-div=true,
 * HIP -ffp-contract=off.
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>

#ifdef TAPS_HIP
#include <hip/hip_runtime.h>
#define cudaError_t hipError_t
#define cudaSuccess hipSuccess
#define cudaGetDeviceCount hipGetDeviceCount
#define cudaSetDevice hipSetDevice
#define cudaDeviceProp hipDeviceProp_t
#define cudaGetDeviceProperties hipGetDeviceProperties
#define cudaMalloc hipMalloc
#define cudaFree hipFree
#define cudaMemcpy hipMemcpy
#define cudaMemcpyHostToDevice hipMemcpyHostToDevice
#define cudaMemcpyDeviceToHost hipMemcpyDeviceToHost
#define cudaMemset hipMemset
#define cudaDeviceSynchronize hipDeviceSynchronize
#define cudaGetLastError hipGetLastError
#define cudaGetErrorString hipGetErrorString
#define PFX(x) taps_hip_##x
#else
#include <cuda_runtime.h>
#define PFX(x) taps_cuda_##x
#endif

#define TAPS_FN __device__
#define TAPS_G
#define taps_log(x) log(x)
#define taps_atomic_inc(p) atomicAdd((p), 1)
#include "baq_hmm_body.h"

__global__ void taps_baq_hmm_k(const taps_gpu_read *rd, int n, int idim_pad,
                               const unsigned char *refb, const unsigned char *qry,
                               const unsigned char *qual8, const float *q2p, double *f, double *scr,
                               short *state_out, unsigned char *q_out, double *max_out, int *nflag)
{
    int r = blockIdx.x * blockDim.x + threadIdx.x;
    if (r < n)
        taps_hmm_read(r, n, idim_pad, rd, refb, qry, qual8, q2p, f, scr, state_out, q_out, max_out, nflag);
}

typedef struct {
    void *p;
    size_t cap;
} dbuf;

typedef struct {
    int dev;
    size_t f_budget;
    dbuf rd, refb, qry, qual, q2p, f, scr, state, q, maxv, nflag;
} ctx_t;

static int dgrow(dbuf *b, size_t want)
{
    if (b->cap >= want && b->p) return 0;
    if (b->p) cudaFree(b->p);
    b->p = NULL;
    size_t n = want < 4096 ? 4096 : want + want / 4;
    if (cudaMalloc(&b->p, n) != cudaSuccess) { b->cap = 0; return -1; }
    b->cap = n;
    return 0;
}

extern "C" void *PFX(open)(const char *want, char *name, int name_len, size_t f_budget)
{
    int nd = 0;
    if (cudaGetDeviceCount(&nd) != cudaSuccess || nd == 0) {
        snprintf(name, name_len, "no GPU device");
        return NULL;
    }
    int pick = -1;
    cudaDeviceProp prop;
    for (int d = 0; d < nd; ++d) {
        if (cudaGetDeviceProperties(&prop, d) != cudaSuccess) continue;
        if (!want || !*want || strstr(prop.name, want)) { pick = d; break; }
    }
    if (pick < 0) {
        snprintf(name, name_len, "no GPU device matching '%s'", want);
        return NULL;
    }
    cudaGetDeviceProperties(&prop, pick);
    if (cudaSetDevice(pick) != cudaSuccess) {
        snprintf(name, name_len, "cannot select device %d", pick);
        return NULL;
    }
    ctx_t *c = (ctx_t *)calloc(1, sizeof(ctx_t));
    c->dev = pick;
    c->f_budget = f_budget;
    snprintf(name, name_len, "%s", prop.name);
    return c;
}

extern "C" void PFX(close)(void *cp)
{
    ctx_t *c = (ctx_t *)cp;
    if (!c) return;
    dbuf *all[] = {&c->rd, &c->refb, &c->qry, &c->qual, &c->q2p, &c->f, &c->scr, &c->state, &c->q, &c->maxv, &c->nflag};
    for (size_t i = 0; i < sizeof(all) / sizeof(all[0]); ++i)
        if (all[i]->p) cudaFree(all[i]->p);
    free(c);
}

#define CK(x) do { cudaError_t e_ = (x); if (e_ != cudaSuccess) { fprintf(stderr, "taps GPU BAQ: %s: %s\n", #x, cudaGetErrorString(e_)); return 1; } } while (0)

extern "C" int PFX(run)(void *cp, const taps_gpu_read *rd, int n, const uint8_t *refb, size_t nref,
                        const uint8_t *qry, const uint8_t *qual, size_t nbase, const float *q2p,
                        int16_t *state, uint8_t *q, double *maxv, int *nflag)
{
    ctx_t *c = (ctx_t *)cp;
    *nflag = 0;
    if (n == 0) return 0;
    CK(cudaSetDevice(c->dev));
    if (dgrow(&c->rd, (size_t)n * sizeof(taps_gpu_read)) || dgrow(&c->refb, nref + 1) || dgrow(&c->qry, nbase + 1)
        || dgrow(&c->qual, nbase + 1) || dgrow(&c->q2p, 256 * sizeof(float)) || dgrow(&c->state, nbase * 2 + 2)
        || dgrow(&c->q, nbase + 1) || dgrow(&c->maxv, nbase * 8 + 8) || dgrow(&c->nflag, sizeof(int)))
        return 2;
    CK(cudaMemcpy(c->rd.p, rd, (size_t)n * sizeof(taps_gpu_read), cudaMemcpyHostToDevice));
    CK(cudaMemcpy(c->refb.p, refb, nref, cudaMemcpyHostToDevice));
    CK(cudaMemcpy(c->qry.p, qry, nbase, cudaMemcpyHostToDevice));
    CK(cudaMemcpy(c->qual.p, qual, nbase, cudaMemcpyHostToDevice));
    CK(cudaMemcpy(c->q2p.p, q2p, 256 * sizeof(float), cudaMemcpyHostToDevice));
    CK(cudaMemset(c->nflag.p, 0, sizeof(int)));
    /* Launch in chunks whose f fits the budget: n_chunk * (max l_query + 1) * idim_pad doubles. */
    int start = 0;
    while (start < n) {
        int lq = 0, id = 0, cnt = 0;
        while (start + cnt < n) {
            int nlq = rd[start + cnt].l_query > lq ? rd[start + cnt].l_query : lq;
            int nid = rd[start + cnt].i_dim > id ? rd[start + cnt].i_dim : id;
            size_t bytes = (size_t)(cnt + 1) * (nlq + 1) * (nid + 3) * 8;
            if (cnt > 0 && bytes > c->f_budget) break;
            lq = nlq; id = nid; ++cnt;
        }
        size_t fb = (size_t)cnt * (lq + 1) * (id + 3) * 8;
        size_t sb = (size_t)cnt * (4 * (id + 3) + lq + 2) * 8;
        if (dgrow(&c->f, fb) || dgrow(&c->scr, sb)) return 2;
        int tpb = 128, blocks = (cnt + tpb - 1) / tpb;
        taps_baq_hmm_k<<<blocks, tpb>>>((const taps_gpu_read *)c->rd.p + start, cnt, id + 3,
                                        (const unsigned char *)c->refb.p, (const unsigned char *)c->qry.p,
                                        (const unsigned char *)c->qual.p, (const float *)c->q2p.p, (double *)c->f.p, (double *)c->scr.p,
                                        (short *)c->state.p, (unsigned char *)c->q.p, (double *)c->maxv.p,
                                        (int *)c->nflag.p);
        CK(cudaGetLastError());
        start += cnt;
    }
    CK(cudaDeviceSynchronize());
    CK(cudaMemcpy(state, c->state.p, nbase * 2, cudaMemcpyDeviceToHost));
    CK(cudaMemcpy(q, c->q.p, nbase, cudaMemcpyDeviceToHost));
    CK(cudaMemcpy(nflag, c->nflag.p, sizeof(int), cudaMemcpyDeviceToHost));
    if (*nflag > 0)
        CK(cudaMemcpy(maxv, c->maxv.p, nbase * 8, cudaMemcpyDeviceToHost));
    return 0;
}
