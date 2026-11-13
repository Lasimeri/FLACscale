/*
 * fu_cutoff.c -- source bandwidth detection (see fu_cutoff.h).
 */
#include "fu_cutoff.h"

#include <fftw3.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

enum { N = 2048, NBIN = N / 2 + 1, BELOW = 8 };
#define BAND_HZ   250.0
#define MIN_HZ    8000.0
#define WALL_DB   20.0

struct fu_cutoff {
    uint32_t      rate;
    size_t        fill, frames;
    double        win[N], buf[N], psd[NBIN];
    double       *frame;
    fftw_complex *spec;
    fftw_plan     plan;
};

fu_cutoff *fu_cutoff_new(uint32_t rate)
{
    fu_cutoff *c = calloc(1, sizeof *c);
    if (!c)
        return NULL;
    c->rate = rate;
    c->frame = fftw_malloc(N * sizeof *c->frame);
    c->spec = fftw_malloc(NBIN * sizeof *c->spec);
    if (!c->frame || !c->spec) {
        fu_cutoff_free(c);
        return NULL;
    }
    c->plan = fftw_plan_dft_r2c_1d(N, c->frame, c->spec, FFTW_ESTIMATE);
    for (int i = 0; i < N; i++)
        c->win[i] = 0.5 - 0.5 * cos(2.0 * M_PI * i / N);
    return c;
}

void fu_cutoff_free(fu_cutoff *c)
{
    if (!c)
        return;
    if (c->plan)
        fftw_destroy_plan(c->plan);
    fftw_free(c->frame);
    fftw_free(c->spec);
    free(c);
}

void fu_cutoff_feed(fu_cutoff *c, const double *x, size_t n)
{
    while (n) {
        const size_t take = n < N - c->fill ? n : N - c->fill;
        memcpy(c->buf + c->fill, x, take * sizeof *x);
        c->fill += take;
        x += take;
        n -= take;
        if (c->fill < N)
            continue;
        for (int i = 0; i < N; i++)
            c->frame[i] = c->buf[i] * c->win[i];
        fftw_execute(c->plan);
        for (int k = 0; k < NBIN; k++)
            c->psd[k] += c->spec[k][0] * c->spec[k][0] + c->spec[k][1] * c->spec[k][1];
        c->frames++;
        c->fill = 0;
    }
}

double fu_cutoff_decide(const fu_cutoff *c, double *drop_db)
{
    const double hz = (double)c->rate / N;
    const int nb = (int)(0.5 * c->rate / BAND_HZ);
    double level[512];

    if (drop_db)
        *drop_db = 0.0;
    if (c->frames == 0 || nb > 512 || nb < 40)
        return 0.0;
    for (int b = 0; b < nb; b++) {
        const int k0 = (int)ceil(b * BAND_HZ / hz), k1 = (int)ceil((b + 1) * BAND_HZ / hz);
        double p = 0.0;
        for (int k = k0; k < k1 && k < NBIN; k++)
            p += c->psd[k];
        level[b] = 10.0 * log10(p / ((k1 - k0) * (double)c->frames) + 1e-300);
    }
    /* loudest[b] = highest level among the bands from b to the top. */
    double loudest[513];
    loudest[nb] = -1e300;
    for (int b = nb - 1; b >= 0; b--)
        loudest[b] = level[b] > loudest[b + 1] ? level[b] : loudest[b + 1];

    for (int b = (int)(MIN_HZ / BAND_HZ); b + 2 <= nb; b++) {
        double below[BELOW];
        memcpy(below, level + b - BELOW, sizeof below);
        for (int i = 1; i < BELOW; i++) {               /* insertion sort, 8 values */
            const double v = below[i];
            int j = i;
            for (; j > 0 && below[j - 1] > v; j--)
                below[j] = below[j - 1];
            below[j] = v;
        }
        const double drop = 0.5 * (below[BELOW / 2 - 1] + below[BELOW / 2]) - loudest[b];
        if (drop >= WALL_DB) {
            if (drop_db)
                *drop_db = drop;
            return b * BAND_HZ;
        }
    }
    return 0.0;
}
