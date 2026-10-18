/*
 * GPU BAQ: one read per work-item, shared by every backend (OpenCL C in
 * baq_hmm.cl, CUDA and HIP in baq_cuda.cu).
 *
 * The arithmetic is htslib 1.19.1 probaln_glocal (probaln.c, default
 * non-PROBALN_ORIG path), operation for operation, as in ../baq_kernel.h:
 * no FMA contraction, IEEE double division, band sums in order. The HMM
 * parameters (m[], sM, bM, bI) are computed on the host with the same float
 * expressions htslib uses (taps_baq_gpu_export), and qualities arrive as the
 * same float probabilities (q2p), so the device only does double + - * / and
 * one log per base.
 *
 * Differences in layout only, not arithmetic:
 *  - f is stored per read in global memory, interleaved across the launch
 *    (cell (row, u) of read r at f[(row * idim_pad + u) * n + r]) so adjacent
 *    work-items touch adjacent addresses.
 *  - The MAP step runs inside the backward pass (row i needs f[i] and b[i]),
 *    so b needs only two rows.
 *  - Each read has its own shape (length, reference span, band): no grouping.
 *  - Rows are zeroed before use (htslib callocs f and b); the code reads at
 *    most 3 cells past a row's i_dim, which htslib finds zero in the next row.
 *
 * q is computed here with the device log. If -4.343*log(1-max)+.499 lies
 * within 1e-9 of an integer (q's rounding point; device and glibc logs differ
 * by a few ulp, ~1e-14 here), or is not finite, q is set to 0xFF and max is
 * stored: the host recomputes q with glibc exactly like the CPU kernel.
 *
 * Before including, define: TAPS_FN (function qualifier), TAPS_G (global
 * address space), taps_log(x), taps_atomic_inc(p).
 */

#define TAPS_GPU_MAXLQ 192
#define TAPS_GPU_MAXIDIM 72
#define TAPS_GPU_EI .25
#define TAPS_GPU_EM .33333333333
#define TAPS_GPU_SET_U(u, b, i, k) { int x_ = (i) - (b); x_ = x_ > 0 ? x_ : 0; (u) = ((k) - x_ + 1) * 3; }

/* 128 bytes; layout identical on host (taps_gpu.h) and device. */
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

/* Precision of the HMM arithmetic. double is exact (htslib); anything else
 * is for timing experiments / an opt-in fast mode (buffers stay double). */
#ifndef TAPS_REAL
#define TAPS_REAL double
#endif

TAPS_FN void taps_hmm_read_generic(int r, int n, int idim_pad,
                           TAPS_G const taps_gpu_read *rd,
                           TAPS_G const unsigned char *refb,
                           TAPS_G const unsigned char *qry,
                           TAPS_G const unsigned char *qual8,
                           TAPS_G const float *q2p,
                           TAPS_G double *f,
                           TAPS_G double *scr,
                           TAPS_G short *state_out,
                           TAPS_G unsigned char *q_out,
                           TAPS_G double *max_out,
                           TAPS_G int *nflag)
{
    TAPS_G const taps_gpu_read *p = rd + r;
    int l_ref = p->l_ref, l_query = p->l_query, bw = p->bw, i_dim = p->i_dim;
    TAPS_G const unsigned char *ref = refb + p->ref_off;
    TAPS_G const unsigned char *query = qry + p->base_off;
    TAPS_G const unsigned char *iq = qual8 + p->base_off;
    unsigned int base = p->base_off;
    TAPS_REAL m[9];
    TAPS_REAL sM = p->sM, sI = p->sM, bM = p->bM, bI = p->bI;
    int i, k, t;
    unsigned nu = (unsigned)n, ru = (unsigned)r;
/* Global f keeps only what MAP reads: the M and I cells of each row's band.
 * Per-read scratch (two forward rows, two backward rows, the scaling factors
 * s[]) lives in scr, interleaved across the launch like f, so neighbouring
 * work-items touch neighbouring addresses. */
/* 32-bit index arithmetic: 64-bit integer multiplies are emulated on some
 * GPUs (Intel Xe). Every index is < 2^32: f holds at most n * (l_query + 1) *
 * idim_pad doubles (65,536 * 193 * 75 < 1e9), scr far fewer. */
#define TIDX(row, u) (((unsigned)(row) * (unsigned)idim_pad + (unsigned)(u)) * nu + ru)
#define TF(row, u) f[TIDX(row, u)]
#define TP(row, u) scr[TIDX((row) & 1, u)]
#define TB(row, u) scr[TIDX(2 + ((row) & 1), u)]
#define S(i) scr[TIDX(4, 0) + (unsigned)(i) * nu]
#define QL(i) ((TAPS_REAL)q2p[iq[(i)]])

    for (t = 0; t < 9; ++t) m[t] = p->m[t];

    /*** forward (rows in fp; band M/I cells copied to global f for MAP) ***/
    for (t = 0; t < i_dim + 3; ++t) TP(0, t) = 0.;
    TAPS_GPU_SET_U(k, bw, 0, 0);
    TP(0, k) = 1.;
    S(0) = 1.;
    {
        TAPS_REAL sum = 0.;
        int beg = 1, end = l_ref < bw + 1 ? l_ref : bw + 1;
        TAPS_REAL q0 = QL(0);
        int qy = query[0];
        for (t = 0; t < i_dim + 3; ++t) TP(1, t) = 0.;
        for (k = beg; k <= end; ++k) {
            int u, rk = ref[k - 1];
            TAPS_REAL e = (rk > 3 || qy > 3) ? 1. : (rk == qy ? 1. - q0 : q0 * TAPS_GPU_EM);
            TAPS_GPU_SET_U(u, bw, 1, k);
            TP(1, u + 0) = e * bM;
            TP(1, u + 1) = TAPS_GPU_EI * bI;
            sum += TP(1, u) + TP(1, u + 1);
            TF(1, u + 0) = TP(1, u + 0);
            TF(1, u + 1) = TP(1, u + 1);
        }
        S(1) = sum;
    }
    for (i = 2; i <= l_query; ++i) {
        TAPS_REAL qli = QL(i - 1);
        int qyi = query[i - 1];
        int beg = 1, end = l_ref, x, u, v11;
        TAPS_REAL sum = 0., M, xm0, xm1, xm2, xm3, xm4, E0, E1, l_x0, l_x2;
        x = i - bw; beg = beg > x ? beg : x;
        x = i + bw; end = end < x ? end : x;
        for (t = 0; t < i_dim + 3; ++t) TP(i, t) = 0.;
        E0 = qli * TAPS_GPU_EM;
        E1 = 1. - qli;
        M = 1. / S(i - 1);
        xm0 = M * m[0];
        xm1 = M * m[3];
        xm2 = M * m[6];
        xm3 = TAPS_GPU_EI * M * m[1];
        xm4 = TAPS_GPU_EI * M * m[4];
        TAPS_GPU_SET_U(u, bw, i, beg);
        TAPS_GPU_SET_U(v11, bw, i - 1, beg - 1);
        l_x0 = m[2] * TP(i, u + 0);
        l_x2 = m[8] * TP(i, u + 2);
        for (k = beg; k <= end; ++k, u += 3, v11 += 3) {
            int rk = ref[k - 1];
            TAPS_REAL e = (rk > 3 || qyi > 3) ? 1. : (rk == qyi ? E1 : E0);
            TAPS_REAL z0 = xm0 * TP(i - 1, v11 + 0);
            TAPS_REAL z1 = xm1 * TP(i - 1, v11 + 1);
            TAPS_REAL z2 = xm2 * TP(i - 1, v11 + 2);
            TAPS_REAL z3 = xm3 * TP(i - 1, v11 + 3);
            TAPS_REAL z4 = xm4 * TP(i - 1, v11 + 4);
            TAPS_REAL a = e * (z0 + z1 + z2);
            TAPS_REAL b1 = z3 + z4;
            TAPS_REAL c = l_x0 + l_x2;
            TP(i, u + 0) = a;
            TP(i, u + 1) = b1;
            TP(i, u + 2) = c;
            TF(i, u + 0) = a;
            TF(i, u + 1) = b1;
            sum += a + b1 + c;
            l_x0 = m[2] * a;
            l_x2 = m[8] * c;
        }
        S(i) = sum;
    }
    {
        TAPS_REAL sum = 0., M = 1. / S(l_query);
        for (k = 1; k <= l_ref; ++k) {
            int u;
            TAPS_GPU_SET_U(u, bw, l_query, k);
            if (u < 3 || u >= i_dim) continue;
            sum += M * TP(l_query, u + 0) * sM + M * TP(l_query, u + 1) * sI;
        }
        S(l_query + 1) = sum;
    }

    /*** backward, with MAP for each finished row ***/
    for (i = l_query; i >= 1; --i) {
        int beg = 1, end = l_ref, x, u;
        x = i - bw; beg = beg > x ? beg : x;
        x = i + bw; end = end < x ? end : x;
        for (t = 0; t < i_dim + 3; ++t) TB(i, t) = 0.;
        if (i == l_query) {
            for (k = 1; k <= l_ref; ++k) {
                TAPS_GPU_SET_U(u, bw, l_query, k);
                if (u < 3 || u >= i_dim) continue;
                TB(i, u + 0) = sM / S(l_query) / S(l_query + 1);
                TB(i, u + 1) = sI / S(l_query) / S(l_query + 1);
            }
        } else {
            TAPS_REAL y = (i > 1);
            TAPS_REAL qli1 = QL(i);
            int qyi1 = query[i];
            TAPS_REAL E0 = qli1 * TAPS_GPU_EM, E1 = 1. - qli1;
            TAPS_REAL e1 = TAPS_GPU_EI * m[1], e4 = TAPS_GPU_EI * m[4];
            TAPS_REAL nrm = 1. / S(i);
            TAPS_REAL xi_5;
            int v10;
            TAPS_GPU_SET_U(u, bw, i, end);
            TAPS_GPU_SET_U(v10, bw, i + 1, end);
            xi_5 = TB(i, u + 5);
            for (k = end; k >= beg; --k, u -= 3, v10 -= 3) {
                TAPS_REAL e, x0, x1, x2;
                if (k >= l_ref) {
                    e = 0.;
                } else {
                    int rk = ref[k];
                    e = ((rk > 3 || qyi1 > 3) ? 1. : (rk == qyi1 ? E1 : E0)) * TB(i + 1, v10 + 3);
                }
                x1 = e * m[3] + e4 * TB(i + 1, v10 + 1);
                x0 = e * m[0] + e1 * TB(i + 1, v10 + 1) + m[2] * xi_5;
                x2 = (e * m[6] + m[8] * xi_5) * y;
                xi_5 = x2;
                TB(i, u + 1) = x1 * nrm;
                TB(i, u + 0) = x0 * nrm;
                TB(i, u + 2) = x2 * nrm;
            }
        }
        /* MAP for row i */
        {
            TAPS_REAL sum = 0., mx = 0., M = 1. / S(i), v;
            int max_k = -1;
            TAPS_GPU_SET_U(u, bw, i, beg);
            for (k = beg; k <= end; ++k, u += 3) {
                TAPS_REAL z1 = M * TF(i, u + 0) * TB(i, u + 0);
                TAPS_REAL z2 = M * TF(i, u + 1) * TB(i, u + 1);
                int which = z2 > z1;
                TAPS_REAL zm = which ? z2 : z1;
                if (zm > mx) {
                    mx = zm;
                    max_k = (k - 1) << 2 | which;
                }
                sum += z1 + z2;
            }
            mx /= sum;
            state_out[base + i - 1] = (short)max_k;
            v = -4.343 * taps_log(1. - mx) + .499;
            if (!(fabs(v - rint(v)) >= 1e-9)) {
                q_out[base + i - 1] = 0xFF;
                max_out[base + i - 1] = mx;
                taps_atomic_inc(nflag);
            } else {
                int kq = (int)v;
                q_out[base + i - 1] = (unsigned char)(kq > 100 ? 99 : kq);
            }
        }
    }
#undef TF
#undef TIDX
#undef TB
#undef TP
#undef S
#undef QL
}


/*
 * Fast path for the common band (bw == 7, i_dim == 51, ~99% of TAPS reads):
 * the same arithmetic in the same order as taps_hmm_read_generic, but the
 * band loops are written out cell by cell so every row index is a constant
 * and the working rows live in registers instead of memory. Cells are
 * band-relative (c = k - x_i + 1, 1..15; c = 0 and 16 stay zero); from one row
 * to the next the band slides by one cell (rows beyond 7) or not at all.
 */
#define TAPS_BW7_ROW 54
/* Written out so every index is a constant (registers, no loop). */
#define Z54(A) { (A)[0] = 0.; (A)[1] = 0.; (A)[2] = 0.; (A)[3] = 0.; (A)[4] = 0.; (A)[5] = 0.; (A)[6] = 0.; (A)[7] = 0.; (A)[8] = 0.; (A)[9] = 0.; (A)[10] = 0.; (A)[11] = 0.; (A)[12] = 0.; (A)[13] = 0.; (A)[14] = 0.; (A)[15] = 0.; (A)[16] = 0.; (A)[17] = 0.; (A)[18] = 0.; (A)[19] = 0.; (A)[20] = 0.; (A)[21] = 0.; (A)[22] = 0.; (A)[23] = 0.; (A)[24] = 0.; (A)[25] = 0.; (A)[26] = 0.; (A)[27] = 0.; (A)[28] = 0.; (A)[29] = 0.; (A)[30] = 0.; (A)[31] = 0.; (A)[32] = 0.; (A)[33] = 0.; (A)[34] = 0.; (A)[35] = 0.; (A)[36] = 0.; (A)[37] = 0.; (A)[38] = 0.; (A)[39] = 0.; (A)[40] = 0.; (A)[41] = 0.; (A)[42] = 0.; (A)[43] = 0.; (A)[44] = 0.; (A)[45] = 0.; (A)[46] = 0.; (A)[47] = 0.; (A)[48] = 0.; (A)[49] = 0.; (A)[50] = 0.; (A)[51] = 0.; (A)[52] = 0.; (A)[53] = 0.; }
#define CP54(D, S_) { (D)[0] = (S_)[0]; (D)[1] = (S_)[1]; (D)[2] = (S_)[2]; (D)[3] = (S_)[3]; (D)[4] = (S_)[4]; (D)[5] = (S_)[5]; (D)[6] = (S_)[6]; (D)[7] = (S_)[7]; (D)[8] = (S_)[8]; (D)[9] = (S_)[9]; (D)[10] = (S_)[10]; (D)[11] = (S_)[11]; (D)[12] = (S_)[12]; (D)[13] = (S_)[13]; (D)[14] = (S_)[14]; (D)[15] = (S_)[15]; (D)[16] = (S_)[16]; (D)[17] = (S_)[17]; (D)[18] = (S_)[18]; (D)[19] = (S_)[19]; (D)[20] = (S_)[20]; (D)[21] = (S_)[21]; (D)[22] = (S_)[22]; (D)[23] = (S_)[23]; (D)[24] = (S_)[24]; (D)[25] = (S_)[25]; (D)[26] = (S_)[26]; (D)[27] = (S_)[27]; (D)[28] = (S_)[28]; (D)[29] = (S_)[29]; (D)[30] = (S_)[30]; (D)[31] = (S_)[31]; (D)[32] = (S_)[32]; (D)[33] = (S_)[33]; (D)[34] = (S_)[34]; (D)[35] = (S_)[35]; (D)[36] = (S_)[36]; (D)[37] = (S_)[37]; (D)[38] = (S_)[38]; (D)[39] = (S_)[39]; (D)[40] = (S_)[40]; (D)[41] = (S_)[41]; (D)[42] = (S_)[42]; (D)[43] = (S_)[43]; (D)[44] = (S_)[44]; (D)[45] = (S_)[45]; (D)[46] = (S_)[46]; (D)[47] = (S_)[47]; (D)[48] = (S_)[48]; (D)[49] = (S_)[49]; (D)[50] = (S_)[50]; (D)[51] = (S_)[51]; (D)[52] = (S_)[52]; (D)[53] = (S_)[53]; }
TAPS_FN void taps_hmm_read_bw7(int r, int n, int idim_pad,
                               TAPS_G const taps_gpu_read *rd,
                               TAPS_G const unsigned char *refb,
                               TAPS_G const unsigned char *qry,
                               TAPS_G const unsigned char *qual8,
                               TAPS_G const float *q2p,
                               TAPS_G double *f,
                               TAPS_G double *scr,
                               TAPS_G short *state_out,
                               TAPS_G unsigned char *q_out,
                               TAPS_G double *max_out,
                               TAPS_G int *nflag)
{
    TAPS_G const taps_gpu_read *p = rd + r;
    int l_ref = p->l_ref, l_query = p->l_query;
    TAPS_G const unsigned char *ref = refb + p->ref_off;
    TAPS_G const unsigned char *query = qry + p->base_off;
    TAPS_G const unsigned char *iq = qual8 + p->base_off;
    unsigned int base = p->base_off;
    TAPS_REAL m0 = p->m[0], m1 = p->m[1], m2 = p->m[2], m3 = p->m[3], m4 = p->m[4], m6 = p->m[6], m8 = p->m[8];
    TAPS_REAL sM = p->sM, sI = p->sM, bM = p->bM, bI = p->bI;
    TAPS_REAL P[TAPS_BW7_ROW], Cr[TAPS_BW7_ROW];
    int i, xi, c_beg, c_end, pxi;
    /* Reference codes for the band, kept in registers and slid one base per
     * row (one load per row instead of one per cell): RW[C] = ref[k - 1]
     * (forward), RB[C] = ref[k] (backward), k = C + xi - 1. */
    int RW[17], RB[17];
#define LDREF(j) (((j) >= 0 && (j) < l_ref) ? (int)ref[(j)] : 4)
#define SHL(A_) A_[1] = A_[2]; A_[2] = A_[3]; A_[3] = A_[4]; A_[4] = A_[5]; A_[5] = A_[6]; A_[6] = A_[7]; A_[7] = A_[8]; \
    A_[8] = A_[9]; A_[9] = A_[10]; A_[10] = A_[11]; A_[11] = A_[12]; A_[12] = A_[13]; A_[13] = A_[14]; A_[14] = A_[15];
#define SHR(A_) A_[15] = A_[14]; A_[14] = A_[13]; A_[13] = A_[12]; A_[12] = A_[11]; A_[11] = A_[10]; A_[10] = A_[9]; A_[9] = A_[8]; \
    A_[8] = A_[7]; A_[7] = A_[6]; A_[6] = A_[5]; A_[5] = A_[4]; A_[4] = A_[3]; A_[3] = A_[2]; A_[2] = A_[1];
    unsigned nu = (unsigned)n, ru = (unsigned)r;
#ifdef TAPS_ABL_NOTF /* timing experiments only: wrong results */
    TAPS_REAL tf_dummy_ = 0.5;
#define TF(row, u) tf_dummy_
#else
#define TF(row, u) f[(((unsigned)(row) * (unsigned)idim_pad + (unsigned)(u)) * nu + ru)]
#endif
#define S(i_) scr[(((unsigned)4 * (unsigned)idim_pad) * nu + ru) + (unsigned)(i_) * nu]
#define QL(i_) ((TAPS_REAL)q2p[iq[(i_)]])
#define BAND(row) { int b_ = 1, e_ = l_ref, x_; xi = (row) - 7; xi = xi > 0 ? xi : 0; \
    x_ = (row) - 7; b_ = b_ > x_ ? b_ : x_; x_ = (row) + 7; e_ = e_ < x_ ? e_ : x_; \
    c_beg = b_ - xi + 1; c_end = e_ - xi + 1; }

    /*** forward ***/
    S(0) = 1.;
    {
        TAPS_REAL sum = 0., q0 = QL(0);
        int qy = query[0];
        Z54(Cr);
        BAND(1);
        RW[1] = 4; RW[2] = LDREF(0); RW[3] = LDREF(1); RW[4] = LDREF(2); RW[5] = LDREF(3); RW[6] = LDREF(4);
        RW[7] = LDREF(5); RW[8] = LDREF(6); RW[9] = LDREF(7); RW[10] = LDREF(8); RW[11] = LDREF(9);
        RW[12] = LDREF(10); RW[13] = LDREF(11); RW[14] = LDREF(12); RW[15] = LDREF(13);
        pxi = 0;
#define R1CELL(C) if ((C) >= c_beg && (C) <= c_end) { \
            int rk = RW[(C)]; \
            TAPS_REAL e = (rk > 3 || qy > 3) ? 1. : (rk == qy ? 1. - q0 : q0 * TAPS_GPU_EM); \
            Cr[3 * (C) + 0] = e * bM; Cr[3 * (C) + 1] = TAPS_GPU_EI * bI; \
            sum += Cr[3 * (C)] + Cr[3 * (C) + 1]; \
            TF(1, 3 * (C) + 0) = Cr[3 * (C) + 0]; TF(1, 3 * (C) + 1) = Cr[3 * (C) + 1]; }
        R1CELL(2) R1CELL(3) R1CELL(4) R1CELL(5) R1CELL(6) R1CELL(7) R1CELL(8) R1CELL(9)
#undef R1CELL
        S(1) = sum;
    }
    for (i = 2; i <= l_query; ++i) {
        TAPS_REAL qli = QL(i - 1);
        int qyi = query[i - 1];
        TAPS_REAL sum = 0., M, xm0, xm1, xm2, xm3, xm4, E0, E1, l_x0, l_x2;
        int t_;
        CP54(P, Cr);
        Z54(Cr);
        BAND(i);
        if (xi > pxi) { SHL(RW) RW[15] = LDREF(15 + xi - 2); pxi = xi; }
        E0 = qli * TAPS_GPU_EM;
        E1 = 1. - qli;
        M = 1. / S(i - 1);
        xm0 = M * m0;
        xm1 = M * m3;
        xm2 = M * p->m[6];
        xm3 = TAPS_GPU_EI * M * m1;
        xm4 = TAPS_GPU_EI * M * m4;
        l_x0 = m2 * 0.;
        l_x2 = m8 * 0.;
#define FCELL(C, PC) if ((C) >= c_beg && (C) <= c_end) { \
            int rk = RW[(C)]; \
            TAPS_REAL e = (rk > 3 || qyi > 3) ? 1. : (rk == qyi ? E1 : E0); \
            TAPS_REAL z0 = xm0 * P[3 * (PC) + 0]; \
            TAPS_REAL z1 = xm1 * P[3 * (PC) + 1]; \
            TAPS_REAL z2 = xm2 * P[3 * (PC) + 2]; \
            TAPS_REAL z3 = xm3 * P[3 * (PC) + 3]; \
            TAPS_REAL z4 = xm4 * P[3 * (PC) + 4]; \
            TAPS_REAL a = e * (z0 + z1 + z2); \
            TAPS_REAL b1 = z3 + z4; \
            TAPS_REAL c3 = l_x0 + l_x2; \
            Cr[3 * (C) + 0] = a; Cr[3 * (C) + 1] = b1; Cr[3 * (C) + 2] = c3; \
            TF(i, 3 * (C) + 0) = a; TF(i, 3 * (C) + 1) = b1; \
            sum += a + b1 + c3; \
            l_x0 = m2 * a; \
            l_x2 = m8 * c3; }
        if (i > 7) {
            FCELL(1, 1) FCELL(2, 2) FCELL(3, 3) FCELL(4, 4) FCELL(5, 5) FCELL(6, 6) FCELL(7, 7) FCELL(8, 8)
            FCELL(9, 9) FCELL(10, 10) FCELL(11, 11) FCELL(12, 12) FCELL(13, 13) FCELL(14, 14) FCELL(15, 15)
        } else {
            FCELL(2, 1) FCELL(3, 2) FCELL(4, 3) FCELL(5, 4) FCELL(6, 5) FCELL(7, 6) FCELL(8, 7)
            FCELL(9, 8) FCELL(10, 9) FCELL(11, 10) FCELL(12, 11) FCELL(13, 12) FCELL(14, 13) FCELL(15, 14)
        }
#undef FCELL
        S(i) = sum;
    }
    {
        /* Cr is row l_query. c = 1..16 covers every k htslib sums (u in [3, i_dim)). */
        TAPS_REAL sum = 0., M = 1. / S(l_query);
        BAND(l_query);
#define LCELL(C) { int k = (C) + xi - 1; if (k >= 1 && k <= l_ref) sum += M * Cr[3 * (C) + 0] * sM + M * Cr[3 * (C) + 1] * sI; }
        LCELL(1) LCELL(2) LCELL(3) LCELL(4) LCELL(5) LCELL(6) LCELL(7) LCELL(8)
        LCELL(9) LCELL(10) LCELL(11) LCELL(12) LCELL(13) LCELL(14) LCELL(15) LCELL(16)
#undef LCELL
        S(l_query + 1) = sum;
    }

    /*** backward, with MAP for each finished row ***/
    for (i = l_query; i >= 1; --i) {
        int t_;
        if (i == l_query) {
            TAPS_REAL v0 = sM / S(l_query) / S(l_query + 1), v1 = sI / S(l_query) / S(l_query + 1);
            Z54(Cr);
            BAND(i);
            RB[1] = LDREF(1 + xi - 1); RB[2] = LDREF(2 + xi - 1); RB[3] = LDREF(3 + xi - 1); RB[4] = LDREF(4 + xi - 1);
            RB[5] = LDREF(5 + xi - 1); RB[6] = LDREF(6 + xi - 1); RB[7] = LDREF(7 + xi - 1); RB[8] = LDREF(8 + xi - 1);
            RB[9] = LDREF(9 + xi - 1); RB[10] = LDREF(10 + xi - 1); RB[11] = LDREF(11 + xi - 1); RB[12] = LDREF(12 + xi - 1);
            RB[13] = LDREF(13 + xi - 1); RB[14] = LDREF(14 + xi - 1); RB[15] = LDREF(15 + xi - 1);
            pxi = xi;
#define ICELL(C) { int k = (C) + xi - 1; if (k >= 1 && k <= l_ref) { Cr[3 * (C) + 0] = v0; Cr[3 * (C) + 1] = v1; } }
            ICELL(1) ICELL(2) ICELL(3) ICELL(4) ICELL(5) ICELL(6) ICELL(7) ICELL(8)
            ICELL(9) ICELL(10) ICELL(11) ICELL(12) ICELL(13) ICELL(14) ICELL(15) ICELL(16)
#undef ICELL
        } else {
            TAPS_REAL y = (i > 1);
            TAPS_REAL qli1 = QL(i);
            int qyi1 = query[i];
            TAPS_REAL E0 = qli1 * TAPS_GPU_EM, E1 = 1. - qli1;
            TAPS_REAL e1 = TAPS_GPU_EI * m1, e4 = TAPS_GPU_EI * m4;
            TAPS_REAL nrm = 1. / S(i);
            TAPS_REAL xi_5 = 0.;
            CP54(P, Cr);
            Z54(Cr);
            BAND(i);
            if (xi < pxi) { SHR(RB) RB[1] = LDREF(xi); pxi = xi; }
#define BCELL(C, NC) if ((C) >= c_beg && (C) <= c_end) { \
                int k = (C) + xi - 1; \
                TAPS_REAL e, x0, x1, x2; \
                if (k >= l_ref) { \
                    e = 0.; \
                } else { \
                    int rk = RB[(C)]; \
                    e = ((rk > 3 || qyi1 > 3) ? 1. : (rk == qyi1 ? E1 : E0)) * P[3 * (NC) + 3]; \
                } \
                x1 = e * m3 + e4 * P[3 * (NC) + 1]; \
                x0 = e * m0 + e1 * P[3 * (NC) + 1] + m2 * xi_5; \
                x2 = (e * m6 + m8 * xi_5) * y; \
                xi_5 = x2; \
                Cr[3 * (C) + 1] = x1 * nrm; \
                Cr[3 * (C) + 0] = x0 * nrm; \
                Cr[3 * (C) + 2] = x2 * nrm; }
#ifdef TAPS_ABL_NOBWD /* timing experiments only */
            if (0) {
#else
            if (i >= 7) {
#endif
                BCELL(15, 14) BCELL(14, 13) BCELL(13, 12) BCELL(12, 11) BCELL(11, 10) BCELL(10, 9) BCELL(9, 8)
                BCELL(8, 7) BCELL(7, 6) BCELL(6, 5) BCELL(5, 4) BCELL(4, 3) BCELL(3, 2) BCELL(2, 1) BCELL(1, 0)
            } else {
                BCELL(15, 15) BCELL(14, 14) BCELL(13, 13) BCELL(12, 12) BCELL(11, 11) BCELL(10, 10) BCELL(9, 9)
                BCELL(8, 8) BCELL(7, 7) BCELL(6, 6) BCELL(5, 5) BCELL(4, 4) BCELL(3, 3) BCELL(2, 2)
            }
#undef BCELL
        }
        /* MAP for row i */
        {
            TAPS_REAL sum = 0., mx = 0., M = 1. / S(i), v;
            int max_k = -1;
            BAND(i);
#define MCELL(C) if ((C) >= c_beg && (C) <= c_end) { \
                int k = (C) + xi - 1; \
                TAPS_REAL z1 = M * TF(i, 3 * (C) + 0) * Cr[3 * (C) + 0]; \
                TAPS_REAL z2 = M * TF(i, 3 * (C) + 1) * Cr[3 * (C) + 1]; \
                int which = z2 > z1; \
                TAPS_REAL zm = which ? z2 : z1; \
                if (zm > mx) { mx = zm; max_k = (k - 1) << 2 | which; } \
                sum += z1 + z2; }
            MCELL(1) MCELL(2) MCELL(3) MCELL(4) MCELL(5) MCELL(6) MCELL(7) MCELL(8)
            MCELL(9) MCELL(10) MCELL(11) MCELL(12) MCELL(13) MCELL(14) MCELL(15)
#undef MCELL
            mx /= sum;
            state_out[base + i - 1] = (short)max_k;
            v = -4.343 * taps_log(1. - mx) + .499;
            if (!(fabs(v - rint(v)) >= 1e-9)) {
                q_out[base + i - 1] = 0xFF;
                max_out[base + i - 1] = mx;
                taps_atomic_inc(nflag);
            } else {
                int kq = (int)v;
                q_out[base + i - 1] = (unsigned char)(kq > 100 ? 99 : kq);
            }
        }
    }
#undef TF
#undef S
#undef QL
#undef BAND
#undef LDREF
#undef SHL
#undef SHR
}

TAPS_FN void taps_hmm_read(int r, int n, int idim_pad,
                           TAPS_G const taps_gpu_read *rd,
                           TAPS_G const unsigned char *refb,
                           TAPS_G const unsigned char *qry,
                           TAPS_G const unsigned char *qual8,
                           TAPS_G const float *q2p,
                           TAPS_G double *f,
                           TAPS_G double *scr,
                           TAPS_G short *state_out,
                           TAPS_G unsigned char *q_out,
                           TAPS_G double *max_out,
                           TAPS_G int *nflag)
{
#ifndef TAPS_GPU_GENERIC_ONLY /* timing comparisons */
    if (rd[r].bw == 7 && rd[r].i_dim == 51)
#else
    if (0)
#endif
        taps_hmm_read_bw7(r, n, idim_pad, rd, refb, qry, qual8, q2p, f, scr, state_out, q_out, max_out, nflag);
    else
        taps_hmm_read_generic(r, n, idim_pad, rd, refb, qry, qual8, q2p, f, scr, state_out, q_out, max_out, nflag);
}
