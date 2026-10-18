#ifdef TAPS_SGS /* timing experiments: force the sub-group (SIMD) width */
__attribute__((intel_reqd_sub_group_size(TAPS_SGS)))
#endif
__kernel void taps_baq_hmm(__global const taps_gpu_read *rd, int n, int idim_pad,
                           __global const unsigned char *refb,
                           __global const unsigned char *qry,
                           __global const unsigned char *qual8,
                           __global const float *q2p,
                           __global double *f,
                           __global double *scr,
                           __global short *state_out,
                           __global unsigned char *q_out,
                           __global double *max_out,
                           __global int *nflag)
{
    int r = (int)get_global_id(0);
    if (r < n)
        taps_hmm_read(r, n, idim_pad, rd, refb, qry, qual8, q2p, f, scr, state_out, q_out, max_out, nflag);
}
