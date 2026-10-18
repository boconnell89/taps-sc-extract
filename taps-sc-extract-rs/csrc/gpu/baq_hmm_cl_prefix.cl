/* OpenCL C prelude for baq_hmm_body.h (the Rust side concatenates prelude,
 * body and baq_hmm_cl_kernel.cl and builds them at run time). */
#pragma OPENCL EXTENSION cl_khr_fp64 : enable
#pragma OPENCL FP_CONTRACT OFF
#define TAPS_FN
#define TAPS_G __global
#ifdef TAPS_TIMING_NOLOG /* timing experiments only: wrong results */
#define taps_log(x) (x)
#else
#define taps_log(x) log(x)
#endif
#define taps_atomic_inc(p) atomic_inc(p)
