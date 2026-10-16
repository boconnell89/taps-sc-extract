/*
 * Lane-parallel copy of htslib 1.19.1 probaln_glocal (probaln.c, the default
 * non-PROBALN_ORIG path), for VW reads that share l_ref, l_query and band.
 *
 * Every lane performs exactly the scalar code's IEEE double operations in the
 * same order (built with -ffp-contract=off, like htslib), so each lane's
 * state[] and q[] equal htslib's for that read. Only what sam_prob_realn uses
 * is computed: the likelihood Pr and b[0] are skipped.
 *
 * Included by baq_simd.c with VW (lanes) and KFN (function name) defined,
 * inside a `#pragma GCC target` region.
 *
 * Original code: Copyright (C) 2003-2006, 2008-2010 Genome Research Ltd.,
 * Heng Li; 2013, 2018 Genome Research Ltd. MIT licence (see probaln.c).
 */

#define KCAT_(a, b) a##b
#define KCAT(a, b) KCAT_(a, b)
#define vd KCAT(KFN, _vd)
#define vi KCAT(KFN, _vi)
#define vsel KCAT(KFN, _vsel)
#define visel KCAT(KFN, _visel)

typedef double vd __attribute__((vector_size(VW * 8)));
typedef long long vi __attribute__((vector_size(VW * 8)));

/* m ? a : b, per lane (m is 0 or all ones). */
static inline __attribute__((always_inline)) vd vsel(vi m, vd a, vd b) { return (vd)(((vi)a & m) | ((vi)b & ~m)); }
static inline __attribute__((always_inline)) vi visel(vi m, vi a, vi b) { return (a & m) | (b & ~m); }

/*
 * in[l]: lane l's translated ref (0-4), query (0-4) and base qualities.
 * Scratch: f, b hold (l_query+1)*i_dim vd each; s holds l_query+2 vd;
 * R holds l_ref vi; Q holds l_query vi; QL holds l_query vd.
 * Output: state[l][0..l_query), q[l][0..l_query).
 */
static void KFN(const taps_lane_in *in, int l_ref, int l_query, int bw_in, double cd, double ce,
                vd *f, vd *b, vd *s, vi *R, vi *Q, vd *QL, int **state, uint8_t **q)
{
    double m[9], sI, sM, bI, bM;
    int bw, bw2, i, k, l;
    const vd one = (vd){0} + 1.;
    const vi three = (vi){0} + 3;

    for (k = 0; k < l_ref; ++k)
        for (l = 0; l < VW; ++l) R[k][l] = in[l].ref[k];
    for (i = 0; i < l_query; ++i)
        for (l = 0; l < VW; ++l) {
            Q[i][l] = in[l].query[i];
            QL[i][l] = (double)taps_qual2prob[in[l].iqual ? in[l].iqual[i] : 30];
        }

    /*** initialization ***/
    bw = l_ref > l_query ? l_ref : l_query;
    if (bw > bw_in) bw = bw_in;
    if (bw < abs(l_ref - l_query)) bw = abs(l_ref - l_query);
    bw2 = bw * 2 + 1;
    size_t i_dim = bw2 < l_ref ? (size_t)bw2 * 3 + 6 : (size_t)l_ref * 3 + 6;

    memset(f, 0, (size_t)(l_query + 1) * i_dim * sizeof(vd));
    memset(b, 0, (size_t)(l_query + 1) * i_dim * sizeof(vd));

    sM = sI = 1. / (2 * l_query + 2);
    m[0*3+0] = (1 - cd - cd) * (1 - sM);
    m[0*3+1] = m[0*3+2] = cd * (1 - sM);
    m[1*3+0] = (1 - ce) * (1 - sI);
    m[1*3+1] = ce * (1 - sI);
    m[1*3+2] = 0.;
    m[2*3+0] = 1 - ce;
    m[2*3+1] = 0.;
    m[2*3+2] = ce;
    bM = (1 - cd) / l_ref;
    bI = cd / l_ref;

    /*** forward ***/
    set_u(k, bw, 0, 0);
    f[0*i_dim+k] = s[0] = one;
    { // f[1]
        vd *fi = &f[1*i_dim], sum = (vd){0};
        int beg = 1, end = l_ref < bw + 1 ? l_ref : bw + 1;
        vi q0 = Q[0];
        vd E0 = QL[0] * EM, E1 = 1. - QL[0];
        for (k = beg; k <= end; ++k) {
            int u;
            vi r = R[k - 1];
            vd e = vsel((vi)(r > three) | (vi)(q0 > three), one, vsel((vi)(r == q0), E1, E0));
            set_u(u, bw, 1, k);
            fi[u+0] = e * bM; fi[u+1] = (vd){0} + EI * bI;
            sum += fi[u] + fi[u+1];
        }
        s[1] = sum;
    }
    for (i = 2; i <= l_query; ++i) {
        vd *fi = &f[i*i_dim], *fi1 = &f[(i-1)*i_dim], sum, qli = QL[i-1];
        int beg = 1, end = l_ref, x;
        vi qyi = Q[i - 1];
        x = i - bw; beg = beg > x ? beg : x;
        x = i + bw; end = end < x ? end : x;
        vd E0 = qli * EM, E1 = 1. - qli;
        vd M = 1. / s[i-1];
        vd xm0 = M * m[0], xm1 = M * m[3], xm2 = M * m[6];
        vd xm3 = EI * M * m[1], xm4 = EI * M * m[4];
        {
            int u, v11;
            set_u(u,   bw, i,   beg);
            set_u(v11, bw, i-1, beg-1);
            vd *xi = &fi[u];
            vd *yi = &fi1[v11];
            vd l_x0 = m[2] * xi[0];
            vd l_x2 = m[8] * xi[2];
            for (k = beg, sum = (vd){0}; k <= end; ++k, xi += 3, yi += 3) {
                vi r = R[k-1];
                vd e = vsel((vi)(r > three) | (vi)(qyi > three), one, vsel((vi)(r == qyi), E1, E0));
                vd z0 = xm0 * yi[0];
                vd z1 = xm1 * yi[1];
                vd z2 = xm2 * yi[2];
                vd z3 = xm3 * yi[3];
                vd z4 = xm4 * yi[4];
                xi[0] = e * (z0 + z1 + z2);
                xi[1] = z3 + z4;
                xi[2] = l_x0 + l_x2;
                sum  += xi[0] + xi[1] + xi[2];
                l_x0 = m[2] * xi[0];
                l_x2 = m[8] * xi[2];
            }
        }
        s[i] = sum;
    }
    { // f[l_query+1]
        vd sum = (vd){0};
        vd M = 1. / s[l_query];
        for (k = 1; k <= l_ref; ++k) {
            int u;
            set_u(u, bw, l_query, k);
            if (u < 3 || u >= (int)i_dim) continue;
            sum += M * f[l_query*i_dim + u+0] * sM + M * f[l_query*i_dim + u+1] * sI;
        }
        s[l_query+1] = sum;
    }
    /*** backward ***/
    for (k = 1; k <= l_ref; ++k) {
        int u;
        vd *bi = &b[l_query*i_dim];
        set_u(u, bw, l_query, k);
        if (u < 3 || u >= (int)i_dim) continue;
        bi[u+0] = sM / s[l_query] / s[l_query+1]; bi[u+1] = sI / s[l_query] / s[l_query+1];
    }
    for (i = l_query - 1; i >= 1; --i) {
        int beg = 1, end = l_ref, x;
        vd *bi = &b[i*i_dim], *bi1 = &b[(i+1)*i_dim], qli1 = QL[i];
        double y = (i > 1);
        vi qyi1 = Q[i];
        x = i - bw; beg = beg > x ? beg : x;
        x = i + bw; end = end < x ? end : x;
        vd E0 = qli1 * EM, E1 = 1. - qli1;
        {
            int u, v10;
            set_u(u,   bw, i,   end);
            set_u(v10, bw, i+1, end);
            vd *xi = &bi[u];
            vd *yi = &bi1[v10];
            vd xi_5 = xi[5];
            double e1 = EI * m[1];
            double e4 = EI * m[4];
            vd n = 1. / s[i];
            for (k = end; k >= beg; --k, xi -= 3, yi -= 3) {
                vd e;
                if (k >= l_ref) {
                    e = (vd){0};
                } else {
                    vi r = R[k];
                    e = vsel((vi)(r > three) | (vi)(qyi1 > three), one, vsel((vi)(r == qyi1), E1, E0)) * yi[3];
                }
                xi[1] =  e * m[3] +   e4 * yi[1];
                xi[0] =  e * m[0] +   e1 * yi[1]  + m[2] * xi_5;
                xi[2] = (e * m[6] + m[8] * xi_5) * y;
                xi_5 = xi[2];
                xi[1] *= n;
                xi[0] *= n;
                xi[2] *= n;
            }
        }
    }
    /*** MAP ***/
    for (i = 1; i <= l_query; ++i) {
        vd sum = (vd){0}, *fi = &f[i*i_dim], *bi = &b[i*i_dim], max = (vd){0};
        vi max_k = (vi){0} - 1;
        int beg = 1, end = l_ref, x;
        x = i - bw; beg = beg > x ? beg : x;
        x = i + bw; end = end < x ? end : x;
        vd M = 1. / s[i];
        int u;
        set_u(u, bw, i, beg);
        for (k = beg; k <= end; ++k, u += 3) {
            vd z1 = M * fi[u+0] * bi[u+0];
            vd z2 = M * fi[u+1] * bi[u+1];
            vi which = (vi)(z2 > z1);
            vd zm = vsel(which, z2, z1);
            vi upd = (vi)(zm > max);
            max = vsel(upd, zm, max);
            max_k = visel(upd, ((vi){0} + ((long long)(k-1) << 2)) | (which & 1), max_k);
            sum += z1 + z2;
        }
        max /= sum;
        for (l = 0; l < VW; ++l) {
            state[l][i-1] = (int)max_k[l];
            int kq = (int)(-4.343 * log(1. - max[l]) + .499);
            q[l][i-1] = kq > 100 ? 99 : kq;
        }
    }
}

#undef vd
#undef vi
#undef vsel
#undef visel
