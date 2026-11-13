/*
 * fu_engine.c -- the signal path (see fu_engine.h).
 */
#include "fu_engine.h"
#include "fu_resample.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

/* A growable set of planar frames between two stages. */
typedef struct planes {
    double *p[FU_MAX_CH];
    size_t  cap;
} planes;

struct fu_engine {
    fu_engine_cfg     cfg;
    fu_declip_stream *declip;
    fu_chain          chain;
    fu_extend        *extend;       /* NULL = stage absent */
    planes            mid;          /* repaired frames at the source rate */
    planes            up;           /* resampled frames before the extension */
};

static const fu_extend_stats no_extension;

fu_engine_cfg fu_engine_cfg_plain(uint32_t rate_in, uint32_t rate_out, uint32_t channels)
{
    fu_engine_cfg cfg = { .rate_in = rate_in, .rate_out = rate_out, .channels = channels };
    for (int ch = 0; ch < FU_MAX_CH; ch++) {
        cfg.hi[ch] = INFINITY;
        cfg.lo[ch] = -INFINITY;
    }
    return cfg;
}

static int reserve(planes *b, uint32_t nch, size_t frames)
{
    if (frames <= b->cap)
        return 0;
    for (uint32_t ch = 0; ch < nch; ch++) {
        double *q = realloc(b->p[ch], frames * sizeof *q);
        if (!q)
            return -1;
        b->p[ch] = q;
    }
    b->cap = frames;
    return 0;
}

fu_engine *fu_engine_new(const fu_engine_cfg *cfg)
{
    fu_engine *e = calloc(1, sizeof *e);
    if (!e)
        return NULL;
    e->cfg = *cfg;
    if (fu_chain_init(&e->chain, cfg->rate_in, cfg->rate_out, cfg->channels, NULL) ||
        !(e->declip = fu_declip_stream_new(cfg->channels, cfg->hi, cfg->lo, NULL))) {
        fu_engine_free(e);
        return NULL;
    }
    if (cfg->src_top > 0.0) {
        const fu_extend_cfg x = fu_extend_cfg_default(cfg->rate_out, cfg->channels, cfg->src_top,
                                                      cfg->ref_top);
        e->extend = fu_extend_new(&x);
    }
    return e;
}

void fu_engine_free(fu_engine *e)
{
    if (!e)
        return;
    fu_declip_stream_free(e->declip);
    fu_chain_free(&e->chain);
    fu_extend_free(e->extend);
    for (int ch = 0; ch < FU_MAX_CH; ch++) {
        free(e->mid.p[ch]);
        free(e->up.p[ch]);
    }
    free(e);
}

int fu_engine_extending(const fu_engine *e)
{
    return e->extend != NULL;
}

size_t fu_engine_bound(const fu_engine *e, size_t n_in)
{
    const size_t n = fu_chain_bound(&e->chain, fu_declip_stream_bound(e->declip, n_in));
    return e->extend ? fu_extend_bound(e->extend, n) : n;
}

/* Repair + resample n_in frames (or flush both when final) into dst. */
static size_t front(fu_engine *e, const double *const in[], size_t n_in, int final,
                    double *const dst[], size_t dst_cap)
{
    const size_t need = fu_declip_stream_bound(e->declip, n_in);
    if (reserve(&e->mid, e->cfg.channels, need))
        return 0;
    const size_t m = final ? fu_declip_stream_flush(e->declip, e->mid.p, need)
                           : fu_declip_stream_run(e->declip, in, n_in, e->mid.p, need);
    size_t n = fu_chain_run(&e->chain, (const double *const *)e->mid.p, m, dst, dst_cap);
    if (final) {
        double *tail[FU_MAX_CH];
        for (uint32_t ch = 0; ch < e->cfg.channels; ch++)
            tail[ch] = dst[ch] + n;
        n += fu_chain_flush(&e->chain, tail, dst_cap - n);
    }
    return n;
}

static size_t step(fu_engine *e, const double *const in[], size_t n_in, int final,
                   double *const out[], size_t out_cap)
{
    if (!e->extend)
        return front(e, in, n_in, final, out, out_cap);

    const size_t need = fu_chain_bound(&e->chain, fu_declip_stream_bound(e->declip, n_in));
    if (reserve(&e->up, e->cfg.channels, need))
        return 0;
    const size_t k = front(e, in, n_in, final, e->up.p, need);
    size_t n = fu_extend_run(e->extend, (const double *const *)e->up.p, k, out, out_cap);
    if (final) {
        double *tail[FU_MAX_CH];
        for (uint32_t ch = 0; ch < e->cfg.channels; ch++)
            tail[ch] = out[ch] + n;
        n += fu_extend_flush(e->extend, tail, out_cap - n);
    }
    return n;
}

size_t fu_engine_run(fu_engine *e, const double *const in[], size_t n_in,
                     double *const out[], size_t out_cap)
{
    return step(e, in, n_in, 0, out, out_cap);
}

size_t fu_engine_flush(fu_engine *e, double *const out[], size_t out_cap)
{
    return step(e, NULL, 0, 1, out, out_cap);
}

const fu_declip_stats *fu_engine_declip_stats(const fu_engine *e)
{
    return fu_declip_stream_stats(e->declip);
}

const fu_extend_stats *fu_engine_extend_stats(const fu_engine *e)
{
    return e->extend ? fu_extend_stats_get(e->extend) : &no_extension;
}
