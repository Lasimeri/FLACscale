/*
 * fu_cutoff.h -- find where a source's real content ends.
 *
 * Purpose : lossy encoders and many masters low-pass the signal (MP3 at
 *           128-192 kb/s: about 16 kHz). The band extension must start there,
 *           not at the nominal Nyquist. This measures the long-term spectrum
 *           of the source and looks for that wall.
 * Method  : averaged power spectrum (FFT 2048, Hann, no overlap) in 250 Hz
 *           bands. A wall is the lowest boundary at or above 8 kHz where every
 *           band above it, up to Nyquist, lies at least 20 dB under the median
 *           of the 2 kHz below it. (Comparing only the neighbourhood of the
 *           boundary is not enough: one loud isolated partial fakes a drop.)
 * Contract: feed one channel (or a mono mix) at the source rate, then decide.
 *           Returns the wall in Hz, or 0 when the spectrum reaches Nyquist
 *           without one.
 */
#ifndef FU_CUTOFF_H
#define FU_CUTOFF_H

#include <stddef.h>
#include <stdint.h>

typedef struct fu_cutoff fu_cutoff;

fu_cutoff *fu_cutoff_new(uint32_t rate);
void       fu_cutoff_feed(fu_cutoff *c, const double *x, size_t n);
double     fu_cutoff_decide(const fu_cutoff *c, double *drop_db);
void       fu_cutoff_free(fu_cutoff *c);

#endif
