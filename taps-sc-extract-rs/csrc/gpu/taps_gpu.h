/* Host-side declarations shared by the GPU BAQ backends. */
#ifndef TAPS_GPU_H
#define TAPS_GPU_H
#include <stddef.h>
#include <stdint.h>

/* Must match baq_hmm_body.h (128 bytes). */
#ifndef TAPS_GPU_READ_DEFINED
#define TAPS_GPU_READ_DEFINED
typedef struct {
    double m[9];
    double sM, bM, bI;
    int l_ref, l_query, bw, i_dim;
    unsigned int ref_off, base_off;
    int pad0, pad1;
} taps_gpu_read;
#endif

#define TAPS_GPU_MAXLQ 192
#define TAPS_GPU_MAXIDIM 72

/*
 * Every backend exposes (prefix taps_ocl_, taps_cuda_, taps_hip_):
 *   void *open(const char *want, char *name, int name_len, size_t f_budget)
 *     want: device name substring or NULL; name: chosen device; NULL on failure
 *     (name then holds the reason).
 *   int run(void *ctx, const taps_gpu_read *rd, int n, const uint8_t *refb,
 *           size_t nref, const uint8_t *qry, const uint8_t *qual, size_t nbase,
 *           const float *q2p, int16_t *state, uint8_t *q, double *maxv,
 *           int *nflag)
 *     rd must be sorted so that consecutive reads have similar l_query;
 *     returns 0 or an error code. maxv is filled only when *nflag > 0.
 *   void close(void *ctx)
 */
#endif
