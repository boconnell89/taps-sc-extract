/*
 * GPU BAQ backend for OpenCL (Intel compute runtime on Arc / Data Center GPUs;
 * also AMD and NVIDIA OpenCL, and PoCL). libOpenCL is loaded at run time, so
 * the binary needs no OpenCL at build time and runs without it. The kernel
 * source (prelude + baq_hmm_body.h + kernel) is passed in by the caller and
 * built with the vendor's compiler; FP_CONTRACT is off and fp64 is required.
 * API: see taps_gpu.h (prefix taps_ocl_).
 */
#define _GNU_SOURCE
#include <dlfcn.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "taps_gpu.h"
#include <time.h>

/* Cumulative timings (seconds) for TAPS_GPU_PROFILE: upload, kernel, download. */
static double prof_up, prof_kern, prof_down;
static double now_s(void)
{
    struct timespec t;
    clock_gettime(CLOCK_MONOTONIC, &t);
    return t.tv_sec + t.tv_nsec * 1e-9;
}
void taps_ocl_prof(double *out)
{
    out[0] = prof_up;
    out[1] = prof_kern;
    out[2] = prof_down;
}

/* Minimal OpenCL 1.2 declarations (the API is stable; avoids needing headers). */
typedef int32_t cl_int;
typedef uint32_t cl_uint;
typedef uint64_t cl_ulong;
typedef cl_ulong cl_bitfield;
typedef cl_bitfield cl_device_type;
typedef cl_uint cl_device_info;
typedef cl_uint cl_platform_info;
typedef cl_uint cl_program_build_info;
typedef cl_bitfield cl_mem_flags;
typedef cl_bitfield cl_command_queue_properties;
typedef intptr_t cl_context_properties;
typedef struct _cl_platform_id *cl_platform_id;
typedef struct _cl_device_id *cl_device_id;
typedef struct _cl_context *cl_context;
typedef struct _cl_command_queue *cl_command_queue;
typedef struct _cl_program *cl_program;
typedef struct _cl_kernel *cl_kernel;
typedef struct _cl_mem *cl_mem;
typedef struct _cl_event *cl_event;
#define CL_SUCCESS 0
#define CL_DEVICE_TYPE_GPU (1 << 2)
#define CL_DEVICE_TYPE_ACCELERATOR (1 << 3)
#define CL_DEVICE_NAME 0x102B
#define CL_DEVICE_VENDOR 0x102C
#define CL_DEVICE_DOUBLE_FP_CONFIG 0x1032
#define CL_PROGRAM_BUILD_LOG 0x1183
#define CL_MEM_READ_WRITE (1 << 0)
#define CL_TRUE 1

#define CLF(X) \
    X(cl_int, clGetPlatformIDs, (cl_uint, cl_platform_id *, cl_uint *)) \
    X(cl_int, clGetDeviceIDs, (cl_platform_id, cl_device_type, cl_uint, cl_device_id *, cl_uint *)) \
    X(cl_int, clGetDeviceInfo, (cl_device_id, cl_device_info, size_t, void *, size_t *)) \
    X(cl_context, clCreateContext, (const cl_context_properties *, cl_uint, const cl_device_id *, void *, void *, cl_int *)) \
    X(cl_command_queue, clCreateCommandQueue, (cl_context, cl_device_id, cl_command_queue_properties, cl_int *)) \
    X(cl_program, clCreateProgramWithSource, (cl_context, cl_uint, const char **, const size_t *, cl_int *)) \
    X(cl_int, clBuildProgram, (cl_program, cl_uint, const cl_device_id *, const char *, void *, void *)) \
    X(cl_int, clGetProgramBuildInfo, (cl_program, cl_device_id, cl_program_build_info, size_t, void *, size_t *)) \
    X(cl_kernel, clCreateKernel, (cl_program, const char *, cl_int *)) \
    X(cl_mem, clCreateBuffer, (cl_context, cl_mem_flags, size_t, void *, cl_int *)) \
    X(cl_int, clReleaseMemObject, (cl_mem)) \
    X(cl_int, clSetKernelArg, (cl_kernel, cl_uint, size_t, const void *)) \
    X(cl_int, clEnqueueWriteBuffer, (cl_command_queue, cl_mem, cl_uint, size_t, size_t, const void *, cl_uint, const cl_event *, cl_event *)) \
    X(cl_int, clEnqueueReadBuffer, (cl_command_queue, cl_mem, cl_uint, size_t, size_t, void *, cl_uint, const cl_event *, cl_event *)) \
    X(cl_int, clEnqueueNDRangeKernel, (cl_command_queue, cl_kernel, cl_uint, const size_t *, const size_t *, const size_t *, cl_uint, const cl_event *, cl_event *)) \
    X(cl_int, clFinish, (cl_command_queue)) \
    X(cl_int, clReleaseKernel, (cl_kernel)) \
    X(cl_int, clReleaseProgram, (cl_program)) \
    X(cl_int, clReleaseCommandQueue, (cl_command_queue)) \
    X(cl_int, clReleaseContext, (cl_context))

#define DECL(ret, name, args) static ret (*p_##name) args;
CLF(DECL)

static int load_cl(char *err, int err_len)
{
    static int loaded = 0;
    if (loaded) return 0;
    void *h = dlopen("libOpenCL.so.1", RTLD_NOW | RTLD_LOCAL);
    if (!h) h = dlopen("libOpenCL.so", RTLD_NOW | RTLD_LOCAL);
    if (!h) {
        snprintf(err, err_len, "OpenCL runtime not found (libOpenCL.so.1)");
        return -1;
    }
#define LOAD(ret, name, args) \
    if (!(p_##name = (ret(*) args)dlsym(h, #name))) { snprintf(err, err_len, "OpenCL symbol %s missing", #name); return -1; }
    CLF(LOAD)
    loaded = 1;
    return 0;
}

typedef struct {
    cl_mem m;
    size_t cap;
} cbuf;

typedef struct {
    cl_context ctx;
    cl_command_queue q;
    cl_program prog;
    cl_kernel k;
    size_t f_budget;
    int q2p_done; /* the qual -> probability table is uploaded once */
    cbuf rd, refb, qry, qual, q2p, f, scr, state, qo, maxv, nflag;
} octx;

static int cgrow(octx *c, cbuf *b, size_t want)
{
    if (b->cap >= want && b->m) return 0;
    if (b->m) p_clReleaseMemObject(b->m);
    b->m = NULL;
    size_t n = want < 4096 ? 4096 : want + want / 4;
    cl_int e;
    b->m = p_clCreateBuffer(c->ctx, CL_MEM_READ_WRITE, n, NULL, &e);
    if (e != CL_SUCCESS || !b->m) { b->m = NULL; b->cap = 0; return -1; }
    b->cap = n;
    return 0;
}

void *taps_ocl_open(const char *src, const char *want, char *name, int name_len, size_t f_budget)
{
    if (load_cl(name, name_len)) return NULL;
    cl_platform_id plats[16];
    cl_uint np = 0;
    if (p_clGetPlatformIDs(16, plats, &np) != CL_SUCCESS || np == 0) {
        snprintf(name, name_len, "no OpenCL platform");
        return NULL;
    }
    cl_device_id dev = NULL;
    char dname[256] = "";
    for (cl_uint pi = 0; pi < np && !dev; ++pi) {
        cl_device_id devs[16];
        cl_uint nd = 0;
        if (p_clGetDeviceIDs(plats[pi], CL_DEVICE_TYPE_GPU | CL_DEVICE_TYPE_ACCELERATOR, 16, devs, &nd) != CL_SUCCESS)
            continue;
        for (cl_uint di = 0; di < nd; ++di) {
            cl_ulong fp64 = 0;
            char nm[256] = "";
            p_clGetDeviceInfo(devs[di], CL_DEVICE_NAME, sizeof nm, nm, NULL);
            p_clGetDeviceInfo(devs[di], CL_DEVICE_DOUBLE_FP_CONFIG, sizeof fp64, &fp64, NULL);
            if (!fp64) continue; /* exact BAQ needs native double precision */
            if (want && *want && !strstr(nm, want)) continue;
            dev = devs[di];
            snprintf(dname, sizeof dname, "%s", nm);
            break;
        }
    }
    if (!dev) {
        snprintf(name, name_len, "no OpenCL GPU with fp64%s%s", want && *want ? " matching " : "", want && *want ? want : "");
        return NULL;
    }
    octx *c = calloc(1, sizeof(octx));
    cl_int e;
    c->f_budget = f_budget;
    c->ctx = p_clCreateContext(NULL, 1, &dev, NULL, NULL, &e);
    if (e != CL_SUCCESS) { snprintf(name, name_len, "clCreateContext failed (%d)", e); free(c); return NULL; }
    c->q = p_clCreateCommandQueue(c->ctx, dev, 0, &e);
    if (e != CL_SUCCESS) { snprintf(name, name_len, "clCreateCommandQueue failed (%d)", e); free(c); return NULL; }
    size_t len = strlen(src);
    c->prog = p_clCreateProgramWithSource(c->ctx, 1, &src, &len, &e);
    if (e != CL_SUCCESS) { snprintf(name, name_len, "clCreateProgramWithSource failed (%d)", e); free(c); return NULL; }
    /* No -cl-fast-relaxed-math / -cl-mad-enable: the arithmetic must stay IEEE. */
    /* Intel GPUs: 256 registers per thread (the kernel keeps two HMM rows in
     * registers; measured fastest on Arc Pro B70). TAPS_OPENCL_OPTS replaces
     * the options (testing; never fast-math options). */
    char vendor[256] = "";
    p_clGetDeviceInfo(dev, CL_DEVICE_VENDOR, sizeof vendor, vendor, NULL);
    const char *opts = getenv("TAPS_OPENCL_OPTS");
    if (!opts) opts = strstr(vendor, "Intel") ? "-cl-intel-256-GRF-per-thread" : "";
    e = p_clBuildProgram(c->prog, 1, &dev, opts, NULL, NULL);
    if (e != CL_SUCCESS) {
        char log[4096] = "";
        p_clGetProgramBuildInfo(c->prog, dev, CL_PROGRAM_BUILD_LOG, sizeof log - 1, log, NULL);
        snprintf(name, name_len, "OpenCL build failed (%d): %.400s", e, log);
        free(c);
        return NULL;
    }
    c->k = p_clCreateKernel(c->prog, "taps_baq_hmm", &e);
    if (e != CL_SUCCESS) { snprintf(name, name_len, "clCreateKernel failed (%d)", e); free(c); return NULL; }
    snprintf(name, name_len, "%s", dname);
    return c;
}

void taps_ocl_close(void *cp)
{
    octx *c = cp;
    if (!c) return;
    cbuf *all[] = {&c->rd, &c->refb, &c->qry, &c->qual, &c->q2p, &c->f, &c->scr, &c->state, &c->qo, &c->maxv, &c->nflag};
    for (size_t i = 0; i < sizeof(all) / sizeof(all[0]); ++i)
        if (all[i]->m) p_clReleaseMemObject(all[i]->m);
    p_clReleaseKernel(c->k);
    p_clReleaseProgram(c->prog);
    p_clReleaseCommandQueue(c->q);
    p_clReleaseContext(c->ctx);
    free(c);
}

#define CK(x) do { cl_int e_ = (x); if (e_ != CL_SUCCESS) { fprintf(stderr, "taps GPU BAQ (OpenCL): %s failed (%d)\n", #x, e_); return 1; } } while (0)

int taps_ocl_run(void *cp, const taps_gpu_read *rd, int n, const uint8_t *refb, size_t nref,
                 const uint8_t *qry, const uint8_t *qual, size_t nbase, const float *q2p,
                 int16_t *state, uint8_t *q, double *maxv, int *nflag)
{
    octx *c = cp;
    *nflag = 0;
    if (n == 0) return 0;
    if (cgrow(c, &c->rd, (size_t)n * sizeof(taps_gpu_read)) || cgrow(c, &c->refb, nref + 1)
        || cgrow(c, &c->qry, nbase + 1) || cgrow(c, &c->qual, nbase + 1) || cgrow(c, &c->q2p, 256 * sizeof(float))
        || cgrow(c, &c->state, nbase * 2 + 2) || cgrow(c, &c->qo, nbase + 1) || cgrow(c, &c->maxv, nbase * 8 + 8)
        || cgrow(c, &c->nflag, sizeof(int)))
        return 2;
    int zero = 0;
    double t0 = now_s();
    /* Non-blocking: one wait per launch instead of one per transfer (each
     * blocking call costs about a millisecond on Intel's runtime). The host
     * arrays outlive the launch. */
    CK(p_clEnqueueWriteBuffer(c->q, c->rd.m, 0, 0, (size_t)n * sizeof(taps_gpu_read), rd, 0, NULL, NULL));
    CK(p_clEnqueueWriteBuffer(c->q, c->refb.m, 0, 0, nref, refb, 0, NULL, NULL));
    CK(p_clEnqueueWriteBuffer(c->q, c->qry.m, 0, 0, nbase, qry, 0, NULL, NULL));
    CK(p_clEnqueueWriteBuffer(c->q, c->qual.m, 0, 0, nbase, qual, 0, NULL, NULL));
    if (!c->q2p_done) {
        CK(p_clEnqueueWriteBuffer(c->q, c->q2p.m, CL_TRUE, 0, 256 * sizeof(float), q2p, 0, NULL, NULL));
        c->q2p_done = 1;
    }
    CK(p_clEnqueueWriteBuffer(c->q, c->nflag.m, 0, 0, sizeof(int), &zero, 0, NULL, NULL));
    double t1 = now_s();
    prof_up += t1 - t0;
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
        if (cgrow(c, &c->f, fb) || cgrow(c, &c->scr, sb)) return 2;
        /* The chunk's records start at rd + start: pass that offset through a
         * sub-copy (OpenCL 1.2 has no pointer offsets on buffers). */
        cbuf sub = {0};
        if (start > 0) {
            if (cgrow(c, &sub, (size_t)cnt * sizeof(taps_gpu_read))) return 2;
            CK(p_clEnqueueWriteBuffer(c->q, sub.m, CL_TRUE, 0, (size_t)cnt * sizeof(taps_gpu_read), rd + start, 0, NULL, NULL));
            /* (rare: only when one launch exceeds the f budget) */
        }
        cl_mem rdm = start > 0 ? sub.m : c->rd.m;
        int idim_pad = id + 3;
        CK(p_clSetKernelArg(c->k, 0, sizeof(cl_mem), &rdm));
        CK(p_clSetKernelArg(c->k, 1, sizeof(int), &cnt));
        CK(p_clSetKernelArg(c->k, 2, sizeof(int), &idim_pad));
        CK(p_clSetKernelArg(c->k, 3, sizeof(cl_mem), &c->refb.m));
        CK(p_clSetKernelArg(c->k, 4, sizeof(cl_mem), &c->qry.m));
        CK(p_clSetKernelArg(c->k, 5, sizeof(cl_mem), &c->qual.m));
        CK(p_clSetKernelArg(c->k, 6, sizeof(cl_mem), &c->q2p.m));
        CK(p_clSetKernelArg(c->k, 7, sizeof(cl_mem), &c->f.m));
        CK(p_clSetKernelArg(c->k, 8, sizeof(cl_mem), &c->scr.m));
        CK(p_clSetKernelArg(c->k, 9, sizeof(cl_mem), &c->state.m));
        CK(p_clSetKernelArg(c->k, 10, sizeof(cl_mem), &c->qo.m));
        CK(p_clSetKernelArg(c->k, 11, sizeof(cl_mem), &c->maxv.m));
        CK(p_clSetKernelArg(c->k, 12, sizeof(cl_mem), &c->nflag.m));
        const char *ls = getenv("TAPS_OPENCL_LOCAL");
        size_t local = ls ? (size_t)atoi(ls) : 64, global = ((size_t)cnt + local - 1) / local * local;
        CK(p_clEnqueueNDRangeKernel(c->q, c->k, 1, NULL, &global, &local, 0, NULL, NULL));
        if (sub.m) {
            CK(p_clFinish(c->q));
            p_clReleaseMemObject(sub.m);
        }
        start += cnt;
    }
    CK(p_clEnqueueReadBuffer(c->q, c->state.m, 0, 0, nbase * 2, state, 0, NULL, NULL));
    CK(p_clEnqueueReadBuffer(c->q, c->qo.m, 0, 0, nbase, q, 0, NULL, NULL));
    CK(p_clEnqueueReadBuffer(c->q, c->nflag.m, 0, 0, sizeof(int), nflag, 0, NULL, NULL));
    CK(p_clFinish(c->q));
    double t2 = now_s();
    prof_kern += t2 - t1; /* now kernel plus queued transfers */
    if (*nflag > 0)
        CK(p_clEnqueueReadBuffer(c->q, c->maxv.m, CL_TRUE, 0, nbase * 8, maxv, 0, NULL, NULL));
    prof_down += now_s() - t2;
    return 0;
}
