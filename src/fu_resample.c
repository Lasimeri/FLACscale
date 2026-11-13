/*
 * fu_resample.c -- rational polyphase resampler (see fu_resample.h).
 *
 * Math    : y[n] = sum_k h[n*M + c - k*L] * x[k],  c = L*N/2.
 *           With n*M + c = q*L + p:  y[n] = sum_j h[p + j*L] * x[q - j], j < N.
 *           c being a multiple of L makes output 0 use phase 0 and puts it at
 *           input time 0: the stage has no delay to compensate afterwards.
 * Layout  : row r of the table holds phase p_r = (r*M) mod L, taps reversed so
 *           the input window is read forwards. Outputs visit rows 0,1,2,...
 *           and wrap at L, so the table is walked sequentially; adv[r] says
 *           whether the window slides one extra frame after row r.
 * Buffer  : buf = [N/2 - 1 zeros][input ...]; output n reads buf[w .. w+N-1].
 *           After a run the consumed prefix is dropped and w returns to 0.
 * Design  : Kaiser window, A = 14.36 * (fstop - fpass) * N + 7.95 dB,
 *           beta = 0.1102 * (A - 8.7)  (Kaiser 1974; valid for A > 50 dB).
 */
#include "fu_resample.h"
#include "fu_dot.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

static double bessel_i0(double x)
{
    const double q = 0.25 * x * x;
    double sum = 1.0, term = 1.0;
    for (int k = 1; k < 1000; k++) {
        term *= q / ((double)k * (double)k);
        sum += term;
        if (term < sum * 1e-21)
            break;
    }
    return sum;
}

double fu_design_atten(uint32_t N, double fpass, double fstop)
{
    return 14.36 * (fstop - fpass) * (double)N + 7.95;
}

void fu_design(uint32_t L, uint32_t N, double fpass, double fstop, double *h)
{
    const size_t len  = (size_t)L * N;
    const double c    = 0.5 * (double)len;
    const double fc   = 0.5 * (fpass + fstop);
    const double beta = 0.1102 * (fu_design_atten(N, fpass, fstop) - 8.7);
    const double i0b  = bessel_i0(beta);

    h[0] = 0.0;
    for (size_t k = 1; k < len; k++) {
        const double t    = (double)k - c;
        const double u    = t / c;
        const double win  = bessel_i0(beta * sqrt(1.0 - u * u)) / i0b;
        const double arg  = M_PI * 2.0 * fc * t / (double)L;
        const double sinc = (t == 0.0) ? 1.0 : sin(arg) / arg;
        h[k] = 2.0 * fc * sinc * win;
    }
}

static void *alloc64(size_t bytes)
{
    return aligned_alloc(64, (bytes + 63) & ~(size_t)63);
}

int fu_stage_init(fu_stage *s, uint32_t L, uint32_t M, uint32_t N,
                  double fpass, double fstop, uint32_t nch)
{
    memset(s, 0, sizeof *s);
    if (L == 0 || M == 0 || M > L || N < 16 || N % 16 != 0 || nch == 0 || nch > FU_MAX_CH)
        return -1;

    s->L = L;
    s->M = M;
    s->N = N;
    s->nch  = nch;
    s->step = M / L;
    s->fill = N / 2 - 1;
    s->cap  = 4 * (size_t)N;

    double *h = malloc((size_t)L * N * sizeof *h);
    s->tab = alloc64((size_t)L * N * sizeof *s->tab);
    s->adv = malloc(L);
    if (!h || !s->tab || !s->adv) {
        free(h);
        fu_stage_free(s);
        return -1;
    }
    fu_design(L, N, fpass, fstop, h);
    for (uint32_t r = 0; r < L; r++) {
        const uint32_t p = (uint32_t)(((uint64_t)r * M) % L);
        for (uint32_t i = 0; i < N; i++)
            s->tab[(size_t)r * N + i] = h[p + (size_t)(N - 1 - i) * L];
        s->adv[r] = (uint8_t)(p + M % L >= L);
    }
    free(h);

    for (uint32_t ch = 0; ch < nch; ch++) {
        s->buf[ch] = calloc(s->cap, sizeof(double));
        if (!s->buf[ch]) {
            fu_stage_free(s);
            return -1;
        }
    }
    return 0;
}

void fu_stage_free(fu_stage *s)
{
    free(s->tab);
    free(s->adv);
    for (uint32_t ch = 0; ch < FU_MAX_CH; ch++)
        free(s->buf[ch]);
    memset(s, 0, sizeof *s);
}

size_t fu_stage_bound(const fu_stage *s, size_t n_in)
{
    return (size_t)(((uint64_t)n_in * s->L) / s->M) + 2;
}

/* Position of the next output: table row and window start. */
typedef struct cursor {
    uint32_t row;
    size_t   w;
} cursor;

static size_t kern2(const fu_stage *s, const double *xl, const double *xr,
                    double *yl, double *yr, size_t out_cap, cursor *cur)
{
    const uint32_t L = s->L, N = s->N, step = s->step;
    const double  *tab = s->tab;
    const uint8_t *adv = s->adv;
    uint32_t row = cur->row;
    size_t   w = cur->w, n = 0;

    if (s->fill < N)
        return 0;
    const size_t last = s->fill - N;
    while (n < out_cap && w <= last) {
        fu_dot2(tab + (size_t)row * N, xl + w, xr + w, N, yl + n, yr + n);
        w += step + adv[row];
        row = (row + 1 == L) ? 0 : row + 1;
        n++;
    }
    cur->row = row;
    cur->w = w;
    return n;
}

static size_t kern1(const fu_stage *s, const double *x, double *y, size_t out_cap, cursor *cur)
{
    const uint32_t L = s->L, N = s->N, step = s->step;
    const double  *tab = s->tab;
    const uint8_t *adv = s->adv;
    uint32_t row = cur->row;
    size_t   w = cur->w, n = 0;

    if (s->fill < N)
        return 0;
    const size_t last = s->fill - N;
    while (n < out_cap && w <= last) {
        y[n] = fu_dot1(tab + (size_t)row * N, x + w, N);
        w += step + adv[row];
        row = (row + 1 == L) ? 0 : row + 1;
        n++;
    }
    cur->row = row;
    cur->w = w;
    return n;
}

/* Append n_in frames (zeros when in == NULL), emit every output whose window is
 * complete, drop the consumed prefix. */
static size_t stage_push(fu_stage *s, const double *const in[], size_t n_in,
                         double *const out[], size_t out_cap)
{
    if (s->fill + n_in > s->cap) {
        const size_t cap = 2 * (s->fill + n_in);
        for (uint32_t ch = 0; ch < s->nch; ch++) {
            double *b = realloc(s->buf[ch], cap * sizeof *b);
            if (!b)
                return 0;
            s->buf[ch] = b;
        }
        s->cap = cap;
    }
    for (uint32_t ch = 0; ch < s->nch; ch++) {
        if (in)
            memcpy(s->buf[ch] + s->fill, in[ch], n_in * sizeof(double));
        else
            memset(s->buf[ch] + s->fill, 0, n_in * sizeof(double));
    }
    s->fill += n_in;

    const cursor start = { s->row, s->w };
    cursor end = start;
    size_t n = 0;
    uint32_t ch = 0;
    for (; ch + 1 < s->nch; ch += 2) {
        end = start;
        n = kern2(s, s->buf[ch], s->buf[ch + 1], out[ch], out[ch + 1], out_cap, &end);
    }
    if (ch < s->nch) {
        end = start;
        n = kern1(s, s->buf[ch], out[ch], out_cap, &end);
    }
    s->row = end.row;

    const size_t drop = end.w < s->fill ? end.w : s->fill;
    for (ch = 0; ch < s->nch; ch++)
        memmove(s->buf[ch], s->buf[ch] + drop, (s->fill - drop) * sizeof(double));
    s->fill -= drop;
    s->w = end.w - drop;
    return n;
}

size_t fu_stage_run(fu_stage *s, const double *const in[], size_t n_in,
                    double *const out[], size_t out_cap)
{
    return stage_push(s, in, n_in, out, out_cap);
}

size_t fu_stage_flush(fu_stage *s, double *const out[], size_t out_cap)
{
    return stage_push(s, NULL, s->N / 2, out, out_cap);
}

fu_quality fu_quality_default(void)
{
    return (fu_quality){ .taps_steep = 256, .taps_relaxed = 32, .pass = 0.907 };
}

static uint32_t gcd32(uint32_t a, uint32_t b)
{
    while (b) {
        const uint32_t t = a % b;
        a = b;
        b = t;
    }
    return a;
}

int fu_chain_init(fu_chain *c, uint32_t rate_in, uint32_t rate_out, uint32_t nch,
                  const fu_quality *q)
{
    const fu_quality d = q ? *q : fu_quality_default();
    memset(c, 0, sizeof *c);
    if (rate_in == 0 || rate_out < rate_in || nch == 0 || nch > FU_MAX_CH)
        return -1;
    c->nch = nch;
    c->rate_in = rate_in;
    c->rate_out = rate_out;
    if (rate_out == rate_in)
        return 0;

    if ((uint64_t)rate_out >= 2 * (uint64_t)rate_in) {
        if (fu_stage_init(&c->st[0], 2, 1, d.taps_steep, 0.5 * d.pass, 0.5, nch))
            return -1;
        c->nstage = 1;
        if ((uint64_t)rate_out > 2 * (uint64_t)rate_in) {
            const uint32_t g = gcd32(rate_out, 2 * rate_in);
            /* Stage input is 2 * rate_in: the source Nyquist sits at 0.25, its
             * first image at 0.75. Flat to 0.25 keeps stage 0's response intact. */
            if (fu_stage_init(&c->st[1], rate_out / g, 2 * rate_in / g,
                              d.taps_relaxed, 0.25, 0.75, nch)) {
                fu_chain_free(c);
                return -1;
            }
            c->nstage = 2;
        }
    } else {
        const uint32_t g = gcd32(rate_out, rate_in);
        if (fu_stage_init(&c->st[0], rate_out / g, rate_in / g, d.taps_steep,
                          0.5 * d.pass, 0.5, nch))
            return -1;
        c->nstage = 1;
    }
    return 0;
}

void fu_chain_free(fu_chain *c)
{
    fu_stage_free(&c->st[0]);
    fu_stage_free(&c->st[1]);
    for (uint32_t ch = 0; ch < FU_MAX_CH; ch++)
        free(c->mid[ch]);
    memset(c, 0, sizeof *c);
}

size_t fu_chain_bound(const fu_chain *c, size_t n_in)
{
    size_t n = n_in;
    for (uint32_t i = 0; i < c->nstage; i++)
        n = fu_stage_bound(&c->st[i], n + c->st[i].N / 2);
    return n;
}

static int mid_reserve(fu_chain *c, size_t frames)
{
    if (frames <= c->mid_cap)
        return 0;
    for (uint32_t ch = 0; ch < c->nch; ch++) {
        double *b = realloc(c->mid[ch], frames * sizeof *b);
        if (!b)
            return -1;
        c->mid[ch] = b;
    }
    c->mid_cap = frames;
    return 0;
}

size_t fu_chain_run(fu_chain *c, const double *const in[], size_t n_in,
                    double *const out[], size_t out_cap)
{
    if (c->nstage == 0) {
        const size_t n = n_in < out_cap ? n_in : out_cap;
        for (uint32_t ch = 0; ch < c->nch; ch++)
            memcpy(out[ch], in[ch], n * sizeof(double));
        return n;
    }
    if (c->nstage == 1)
        return fu_stage_run(&c->st[0], in, n_in, out, out_cap);

    const size_t need = fu_stage_bound(&c->st[0], n_in);
    if (mid_reserve(c, need))
        return 0;
    const size_t m = fu_stage_run(&c->st[0], in, n_in, c->mid, need);
    return fu_stage_run(&c->st[1], (const double *const *)c->mid, m, out, out_cap);
}

size_t fu_chain_flush(fu_chain *c, double *const out[], size_t out_cap)
{
    if (c->nstage == 0)
        return 0;
    if (c->nstage == 1)
        return fu_stage_flush(&c->st[0], out, out_cap);

    const size_t need = fu_stage_bound(&c->st[0], c->st[0].N / 2);
    if (mid_reserve(c, need))
        return 0;
    const size_t m = fu_stage_flush(&c->st[0], c->mid, need);
    size_t n = fu_stage_run(&c->st[1], (const double *const *)c->mid, m, out, out_cap);
    double *tail[FU_MAX_CH];
    for (uint32_t ch = 0; ch < c->nch; ch++)
        tail[ch] = out[ch] + n;
    n += fu_stage_flush(&c->st[1], tail, out_cap - n);
    return n;
}
