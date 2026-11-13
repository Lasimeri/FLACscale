/*
 * fu_engine.h -- the whole signal path as one streaming object.
 *
 * Path    : source-rate frames
 *             -> clipped-peak repair   (fu_declip, source rate)
 *             -> resampler chain       (fu_resample)
 *             -> band extension        (fu_extend, output rate)
 *             -> output-rate frames.
 *           The command line tool drives it for files; a real-time player
 *           drives the same object block by block.
 * Contract: planar doubles at +-1.0 in and out. fu_engine_run consumes all
 *           n_in frames; out must hold fu_engine_bound(n_in) frames.
 *           fu_engine_flush is called once after the last run and returns the
 *           tail. Total output = ceil(total input * rate_out / rate_in).
 *           hi / lo: per-channel ceilings from fu_cliplevel (+-infinity = that
 *           side is not clipped; the repair stage then passes samples through).
 *           src_top = 0 switches the band extension off; it is also off when
 *           the output rate has no room for it (fu_engine_extending tells).
 */
#ifndef FU_ENGINE_H
#define FU_ENGINE_H

#include "fu_audio.h"
#include "fu_declip.h"
#include "fu_extend.h"

#include <stddef.h>

typedef struct fu_engine_cfg {
    uint32_t rate_in, rate_out, channels;
    double   hi[FU_MAX_CH], lo[FU_MAX_CH];
    double   src_top, ref_top;      /* Hz, see fu_extend_cfg; src_top 0 = no extension */
} fu_engine_cfg;

typedef struct fu_engine fu_engine;

/* A configuration with every ceiling at infinity and no extension. */
fu_engine_cfg fu_engine_cfg_plain(uint32_t rate_in, uint32_t rate_out, uint32_t channels);

fu_engine *fu_engine_new(const fu_engine_cfg *cfg);
void       fu_engine_free(fu_engine *e);
size_t     fu_engine_bound(const fu_engine *e, size_t n_in);
size_t     fu_engine_run(fu_engine *e, const double *const in[], size_t n_in,
                         double *const out[], size_t out_cap);
size_t     fu_engine_flush(fu_engine *e, double *const out[], size_t out_cap);
int        fu_engine_extending(const fu_engine *e);
const fu_declip_stats *fu_engine_declip_stats(const fu_engine *e);
const fu_extend_stats *fu_engine_extend_stats(const fu_engine *e);

#endif
