/*
 * fu_extend.h -- band extension above the source's real bandwidth.
 *
 * Purpose : fill the spectrum the source never had, at the output rate, by
 *           two mechanisms taken from the Opus decoder and reworked:
 *
 *           1. Folding, src_top .. fold_top (25 kHz). CELT fills a band that
 *              received no bits with a copy of lower-band spectrum, normalised
 *              and rescaled to the band's energy (celt/bands.c). Here the band
 *              energies are not transmitted, so they are extrapolated: the
 *              level of the top real kilohertz, continued along the slope
 *              (dB/kHz) fitted to the top four real kilohertz.
 *
 *           2. Harmonics only, above fold_top. The pitch model of the Opus
 *              concealer says a tonal sound is a comb of partials. Each tonal
 *              peak found in the real band (15 dB above its local floor) gets
 *              its integer multiples h = 2 .. max_harm placed above fold_top
 *              with phase h * phi and amplitude reduced by the source's own
 *              spectral tilt (dB per octave) times log2(h). Noise-like content
 *              produces no peak and therefore nothing above fold_top.
 *
 * Contract: planar doubles in and out at `rate`, same frame count overall,
 *           same timing (the dry signal is delayed inside by the transform
 *           latency and the delay is removed again). The extension is added:
 *           the band below src_top is not processed. Everything written above
 *           src_top is synthesis, not recovery.
 *           out must hold fu_extend_bound(n_in) frames; flush once at the end.
 */
#ifndef FU_EXTEND_H
#define FU_EXTEND_H

#include <stdint.h>

#include <stddef.h>

typedef struct fu_extend_cfg {
    uint32_t rate, channels;
    double   src_top;       /* Hz: where the source band ends and folding begins */
    double   ref_top;       /* Hz: top of the band trusted for levels, slope, peaks and as fold
                               source (below any encoder or resampler roll-off) */
    double   fold_top;      /* Hz: folding stops, harmonics-only begins (25000) */
    double   fold_db;       /* trim on the folded band, dB (0) */
    double   harm_db;       /* trim on the harmonics, dB (0) */
    int      max_harm;      /* highest harmonic order (6) */
} fu_extend_cfg;

typedef struct fu_extend_stats {
    uint64_t frames;        /* transform frames processed (all channels) */
    uint64_t peaks;         /* tonal peaks that produced harmonics */
    uint64_t lobes;         /* harmonic partials written */
    uint64_t capped;        /* frames whose harmonics were scaled down by the -30 dB guard */
} fu_extend_stats;

typedef struct fu_extend fu_extend;

fu_extend_cfg fu_extend_cfg_default(uint32_t rate, uint32_t channels, double src_top,
                                    double ref_top);

/* NULL when the output rate leaves no room above src_top. */
fu_extend *fu_extend_new(const fu_extend_cfg *cfg);
void       fu_extend_free(fu_extend *x);
size_t     fu_extend_bound(const fu_extend *x, size_t n_in);
size_t     fu_extend_run(fu_extend *x, const double *const in[], size_t n_in,
                         double *const out[], size_t out_cap);
size_t     fu_extend_flush(fu_extend *x, double *const out[], size_t out_cap);
const fu_extend_stats *fu_extend_stats_get(const fu_extend *x);

#endif
