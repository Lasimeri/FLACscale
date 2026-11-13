/*
 * fu_declip.c -- clipped-peak repair (see fu_declip.h).
 *
 * Math    : c = [1, lpc[0..P-1]] is the prediction-error filter. For unknown
 *           positions u_0 < ... < u_{G-1} the rows r = u_0 .. u_{G-1} + P are
 *           the ones an unknown enters. With the unknowns set to zero,
 *           d[r] = sum_k c[k] * x[r-k] is what the known samples leave in each
 *           row. Minimising sum_r (d[r] + sum_a c[r-u_a] * z_a)^2 gives
 *             M z = b,  M_ab = sum_r c[r-u_a] c[r-u_b],  b_a = -sum_r c[r-u_a] d[r].
 *           M is symmetric positive definite (c[0] = 1): Cholesky.
 *           When every row exists, M_ab depends only on |u_a - u_b| and is the
 *           autocorrelation of c; rows cut off by the end of the signal are
 *           summed explicitly.
 */
#include "fu_declip.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

enum { P = FU_LPC_ORDER, MAXU = FU_DECLIP_MAXU, SPAN = MAXU * (P + 1) + 2 * P + 2 };

struct fu_declip {
    fu_declip_cfg   cfg;
    fu_declip_stats st;
    size_t  u[MAXU];            /* absolute positions of the unknowns */
    double  level[MAXU];        /* the pinned value each unknown started from */
    double  z[MAXU];            /* current solution */
    uint8_t fixed[MAXU];        /* held at the ceiling by the constraint */
    int     freeidx[MAXU];
    double  rhs[MAXU];
    double  seg[SPAN];          /* local copy of the signal around the cluster */
    double  d[SPAN];
    double  m[MAXU * MAXU];
};

/* Opus's guards with one measured change: the noise floor on lag 0 is -80 dB
 * instead of Opus's -40 dB (+0.4 to +0.9 dB SDR on real music, declipcheck). */
fu_declip_cfg fu_declip_cfg_default(void)
{
    fu_declip_cfg cfg = { .lpc = fu_lpc_cfg_opus(), .passes = 1, .guard = 4.0 };
    cfg.lpc.noise_floor = 1e-8;
    return cfg;
}

fu_declip *fu_declip_new(const fu_declip_cfg *cfg)
{
    fu_declip *d = calloc(1, sizeof *d);
    if (d)
        d->cfg = cfg ? *cfg : fu_declip_cfg_default();
    return d;
}

void fu_declip_free(fu_declip *d)
{
    free(d);
}

const fu_declip_stats *fu_declip_stats_get(const fu_declip *d)
{
    return &d->st;
}

/* In-place Cholesky solve of the n x n system m * x = b (m row-major, lower
 * triangle used). Returns -1 when m is not positive definite. */
static int cholesky_solve(double *m, double *b, int n)
{
    for (int j = 0; j < n; j++) {
        double s = m[j * n + j];
        for (int k = 0; k < j; k++)
            s -= m[j * n + k] * m[j * n + k];
        if (!(s > 1e-300))
            return -1;
        const double piv = sqrt(s);
        m[j * n + j] = piv;
        for (int i = j + 1; i < n; i++) {
            double t = m[i * n + j];
            for (int k = 0; k < j; k++)
                t -= m[i * n + k] * m[j * n + k];
            m[i * n + j] = t / piv;
        }
    }
    for (int i = 0; i < n; i++) {
        double t = b[i];
        for (int k = 0; k < i; k++)
            t -= m[i * n + k] * b[k];
        b[i] = t / m[i * n + i];
    }
    for (int i = n - 1; i >= 0; i--) {
        double t = b[i];
        for (int k = i + 1; k < n; k++)
            t -= m[k * n + i] * b[k];
        b[i] = t / m[i * n + i];
    }
    return 0;
}

/* Solve one cluster of G unknowns with the filter c. base is the absolute
 * index of seg[0]; rows run to rowmax (absolute, inclusive). Fills d->z.
 * Returns -1 when the system is singular. */
static int solve_cluster(fu_declip *d, const double *c, const double *rc, int G, size_t base,
                         size_t rowmax)
{
    const size_t u0 = d->u[0];
    const int nrow = (int)(rowmax - u0) + 1;
    const int full = rowmax >= d->u[G - 1] + P;

    memset(d->fixed, 0, (size_t)G);
    for (int iter = 0; iter < 4; iter++) {
        int nf = 0;
        for (int a = 0; a < G; a++) {
            d->seg[d->u[a] - base] = d->fixed[a] ? d->level[a] : 0.0;
            if (!d->fixed[a])
                d->freeidx[nf++] = a;
        }
        if (nf == 0)
            return 0;
        for (int r = 0; r < nrow; r++) {
            const double *s = d->seg + (u0 - base) + r;
            double acc = 0.0;
            for (int k = 0; k <= P; k++)
                acc += c[k] * s[-k];
            d->d[r] = acc;
        }
        for (int i = 0; i < nf; i++) {
            const size_t ua = d->u[d->freeidx[i]];
            const int ra = (int)(ua - u0);
            const int top = ra + P < nrow - 1 ? ra + P : nrow - 1;
            double acc = 0.0;
            for (int r = ra; r <= top; r++)
                acc += c[r - ra] * d->d[r];
            d->rhs[i] = -acc;
            for (int j = 0; j <= i; j++) {
                const size_t ub = d->u[d->freeidx[j]];
                const int lag = (int)(ua - ub);
                double v = 0.0;
                if (lag <= P) {
                    if (full) {
                        v = rc[lag];
                    } else {
                        for (int r = ra; r <= (int)(ub - u0) + P && r < nrow; r++)
                            v += c[r - ra] * c[r - (int)(ub - u0)];
                    }
                }
                d->m[i * nf + j] = v;
            }
        }
        if (cholesky_solve(d->m, d->rhs, nf))
            return -1;
        int newly = 0;
        for (int i = 0; i < nf; i++) {
            const int a = d->freeidx[i];
            const double v = d->rhs[i];
            /* The floor of each repair is the sample's own pinned value. */
            const int short_of = d->level[a] > 0.0 ? v < d->level[a] : v > d->level[a];
            d->z[a] = short_of ? d->level[a] : v;
            d->fixed[a] |= (uint8_t)short_of;
            newly += short_of;
        }
        if (!newly)
            break;
    }
    return 0;
}

size_t fu_declip_span(fu_declip *d, double *x, size_t n, size_t from, size_t to,
                      double hi, double lo, int final)
{
    size_t i = from;

    while (i < to) {
        if (!(x[i] >= hi || x[i] <= lo)) {
            i++;
            continue;
        }
        /* Gather the cluster: pinned samples no further apart than P. */
        int G = 0;
        size_t last = i, j = i;
        for (; j < n && j - last <= P && G < MAXU; j++) {
            if (x[j] >= hi || x[j] <= lo) {
                d->u[G] = j;
                d->level[G] = x[j];
                G++;
                last = j;
            }
        }
        /* Right-hand context: P rows for the system, half a window for the LPC. */
        if (!final && last + P + FU_LPC_WINDOW / 2 >= n)
            return i;

        const size_t u0 = d->u[0];
        const size_t base = u0 >= P ? u0 - P : 0;
        const size_t rowmax = last + P < n ? last + P : n - 1;
        const size_t pad = P - (u0 - base);         /* zeros standing in for x[<0] */
        double ceiling = 0.0;                       /* largest pinned magnitude in the cluster */
        for (int a = 0; a < G; a++)
            ceiling = fmax(ceiling, fabs(d->level[a]));
        int ok = 1;

        d->st.clusters++;
        d->st.pinned += (uint64_t)G;
        for (int pass = 0; pass < d->cfg.passes && ok; pass++) {
            double lpc[P], c[P + 1], rc[P + 1];
            const size_t mid = (u0 + last) / 2;
            size_t w0 = mid >= FU_LPC_WINDOW / 2 ? mid - FU_LPC_WINDOW / 2 : 0;
            if (w0 + FU_LPC_WINDOW > n)
                w0 = n >= FU_LPC_WINDOW ? n - FU_LPC_WINDOW : 0;
            const size_t wn = n - w0 < FU_LPC_WINDOW ? n - w0 : FU_LPC_WINDOW;
            fu_lpc_analyze(x + w0, (int)wn, &d->cfg.lpc, lpc);

            c[0] = 1.0;
            memcpy(c + 1, lpc, sizeof lpc);
            for (int m = 0; m <= P; m++) {
                double acc = 0.0;
                for (int k = 0; k + m <= P; k++)
                    acc += c[k] * c[k + m];
                rc[m] = acc;
            }
            memset(d->seg, 0, pad * sizeof(double));
            memcpy(d->seg + pad, x + base, (rowmax - base + 1) * sizeof(double));
            if (solve_cluster(d, c, rc, G, base - pad, rowmax)) {
                ok = 0;
                break;
            }
            for (int a = 0; a < G && ok; a++)
                ok = fabs(d->z[a]) <= d->cfg.guard * ceiling;
            for (int a = 0; a < G && ok; a++)
                x[d->u[a]] = d->z[a];
        }
        if (!ok) {
            for (int a = 0; a < G; a++)
                x[d->u[a]] = d->level[a];
            d->st.refused++;
        } else {
            for (int a = 0; a < G; a++) {
                d->st.raised += !d->fixed[a];
                d->st.held += d->fixed[a];
                d->st.peak = fmax(d->st.peak, fabs(d->z[a]));
            }
        }
        i = last + 1;
    }
    return to;
}

/* ---- streaming form ---------------------------------------------------- */

enum { KEEP = FU_LPC_WINDOW };      /* settled frames kept as left-hand context */

struct fu_declip_stream {
    fu_declip *d;
    uint32_t   nch;
    double     hi[FU_DECLIP_MAXCH], lo[FU_DECLIP_MAXCH];
    double    *buf[FU_DECLIP_MAXCH];
    size_t     scan[FU_DECLIP_MAXCH];   /* per channel: first frame not settled yet */
    size_t     cap, fill;
    size_t     out;                     /* frames [0, out) were emitted and are history */
};

fu_declip_stream *fu_declip_stream_new(uint32_t nch, const double *hi, const double *lo,
                                       const fu_declip_cfg *cfg)
{
    fu_declip_stream *s = calloc(1, sizeof *s);
    if (!s || nch == 0 || nch > FU_DECLIP_MAXCH || !(s->d = fu_declip_new(cfg))) {
        free(s);
        return NULL;
    }
    s->nch = nch;
    memcpy(s->hi, hi, nch * sizeof *hi);
    memcpy(s->lo, lo, nch * sizeof *lo);
    return s;
}

void fu_declip_stream_free(fu_declip_stream *s)
{
    if (!s)
        return;
    for (uint32_t ch = 0; ch < FU_DECLIP_MAXCH; ch++)
        free(s->buf[ch]);
    fu_declip_free(s->d);
    free(s);
}

const fu_declip_stats *fu_declip_stream_stats(const fu_declip_stream *s)
{
    return fu_declip_stats_get(s->d);
}

size_t fu_declip_stream_bound(const fu_declip_stream *s, size_t n_in)
{
    return s->fill - s->out + n_in;
}

static size_t stream_step(fu_declip_stream *s, const double *const in[], size_t n_in,
                          double *const out[], size_t out_cap, int final)
{
    if (s->fill + n_in > s->cap) {
        const size_t cap = 2 * (s->fill + n_in) + KEEP;
        for (uint32_t ch = 0; ch < s->nch; ch++) {
            double *b = realloc(s->buf[ch], cap * sizeof *b);
            if (!b)
                return 0;
            s->buf[ch] = b;
        }
        s->cap = cap;
    }
    size_t settled = s->fill + n_in;
    for (uint32_t ch = 0; ch < s->nch; ch++) {
        if (n_in)
            memcpy(s->buf[ch] + s->fill, in[ch], n_in * sizeof(double));
        s->scan[ch] = fu_declip_span(s->d, s->buf[ch], s->fill + n_in, s->scan[ch],
                                     s->fill + n_in, s->hi[ch], s->lo[ch], final);
        settled = s->scan[ch] < settled ? s->scan[ch] : settled;
    }
    s->fill += n_in;

    size_t n = settled - s->out;
    n = n < out_cap ? n : out_cap;
    for (uint32_t ch = 0; ch < s->nch; ch++)
        memcpy(out[ch], s->buf[ch] + s->out, n * sizeof(double));
    s->out += n;

    const size_t drop = s->out > KEEP ? s->out - KEEP : 0;
    for (uint32_t ch = 0; ch < s->nch; ch++) {
        memmove(s->buf[ch], s->buf[ch] + drop, (s->fill - drop) * sizeof(double));
        s->scan[ch] -= drop;
    }
    s->fill -= drop;
    s->out -= drop;
    return n;
}

size_t fu_declip_stream_run(fu_declip_stream *s, const double *const in[], size_t n_in,
                            double *const out[], size_t out_cap)
{
    return stream_step(s, in, n_in, out, out_cap, 0);
}

size_t fu_declip_stream_flush(fu_declip_stream *s, double *const out[], size_t out_cap)
{
    return stream_step(s, NULL, 0, out, out_cap, 1);
}
