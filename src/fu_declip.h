/*
 * fu_declip.h -- clipped-peak repair: the Opus concealer, given both sides.
 *
 * Model   : the one celt_decode_lost() uses (libopus 1.6.1): a 24th-order LPC
 *           filter driven by an excitation, e[n] = x[n] + sum lpc[k-1]*x[n-k].
 *           Opus knows only the past of a lost frame, so it extends the
 *           excitation and runs the synthesis filter forward. A clipped run has
 *           intact samples on both sides, so here the pinned samples are the
 *           unknowns of one least-squares problem: minimise the excitation
 *           energy over every row the unknowns touch. With no right-hand side
 *           the solution is the Opus synthesis recursion itself.
 * Extra   : each pinned sample carries information Opus never has: the true
 *           value is at or beyond the ceiling, with known sign. A solution
 *           that falls short is held at the ceiling and the others are
 *           re-solved (active set).
 * Contract: a sample is pinned when x >= hi or x <= lo. Pinned samples closer
 *           than the LPC order are solved together (a cluster, at most
 *           FU_DECLIP_MAXU unknowns). fu_declip_span repairs, in place, the
 *           clusters that start in [from, to) using x[0..n) as context, and
 *           returns the index where it stopped: `to`, or earlier when a
 *           cluster needs samples beyond n and `final` is 0.
 */
#ifndef FU_DECLIP_H
#define FU_DECLIP_H

#include "fu_lpc.h"

#include <stddef.h>
#include <stdint.h>

#define FU_DECLIP_MAXU  256
#define FU_DECLIP_MAXCH 8

typedef struct fu_declip_cfg {
    fu_lpc_cfg lpc;
    int        passes;      /* LPC estimates per cluster; pass 2 re-fits on the repaired signal.
                               Default 1: measured +0.05 dB for 2x the time (declipcheck). */
    double     guard;       /* refuse a repair that exceeds guard * ceiling (4.0 = +12 dB) */
} fu_declip_cfg;

typedef struct fu_declip_stats {
    uint64_t pinned;        /* pinned samples met */
    uint64_t clusters;
    uint64_t raised;        /* samples moved beyond the ceiling */
    uint64_t held;          /* samples the constraint kept at the ceiling */
    uint64_t refused;       /* clusters left untouched (guard or singular system) */
    double   peak;          /* largest magnitude written */
} fu_declip_stats;

typedef struct fu_declip fu_declip;
typedef struct fu_declip_stream fu_declip_stream;

/* Streaming form: planar frames in, repaired frames out, same count overall.
 * hi[ch] / lo[ch] are the ceilings per channel (+-infinity = side not clipped).
 * A frame is emitted once every channel has settled it; a cluster near the end
 * of the buffered signal waits for its right-hand context (about 536 frames).
 * out must hold fu_declip_stream_bound(n_in) frames. */
fu_declip_stream *fu_declip_stream_new(uint32_t nch, const double *hi, const double *lo,
                                       const fu_declip_cfg *cfg);
void   fu_declip_stream_free(fu_declip_stream *s);
size_t fu_declip_stream_bound(const fu_declip_stream *s, size_t n_in);
size_t fu_declip_stream_run(fu_declip_stream *s, const double *const in[], size_t n_in,
                            double *const out[], size_t out_cap);
size_t fu_declip_stream_flush(fu_declip_stream *s, double *const out[], size_t out_cap);
const fu_declip_stats *fu_declip_stream_stats(const fu_declip_stream *s);

fu_declip_cfg fu_declip_cfg_default(void);
fu_declip    *fu_declip_new(const fu_declip_cfg *cfg);
void          fu_declip_free(fu_declip *d);
size_t        fu_declip_span(fu_declip *d, double *x, size_t n, size_t from, size_t to,
                             double hi, double lo, int final);
const fu_declip_stats *fu_declip_stats_get(const fu_declip *d);

#endif
