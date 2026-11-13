/*
 * flacspec.c -- measurement tool: where the power of a FLAC file sits.
 *
 * Purpose : Welch-averaged power spectrum of a whole file (all channels
 *           averaged in power), reported per frequency band: band power and
 *           the strongest spectral line inside the band.
 * Method  : frames of N samples, hop N/8 (the window is narrow in time; a
 *           coarser hop would skip signal), Kaiser window with beta 23.3
 *           (side lobes below -220 dB, so a full-scale tone does not leak
 *           into a band that is 170 dB down). One-sided power, normalised so
 *           the bands sum to the mean square of the signal.
 * Units   : dBFS = 10*log10(power / 0.5): 0 dBFS is a full-scale sine.
 *           A line's level sums the 21 bins around its peak (the main lobe).
 * Usage   : flacspec FILE [-n FFTSIZE] [-l LO HI]   (default FFTSIZE 65536)
 *           -l LO HI lists the 12 strongest spectral lines between LO and HI Hz
 *           and the power that remains in that range once they are removed.
 */
#include "fu_input.h"

#include <fftw3.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static double bessel_i0(double x)
{
    const double q = 0.25 * x * x;
    double sum = 1.0, term = 1.0;
    for (int k = 1; k < 1000 && term > sum * 1e-21; k++) {
        term *= q / ((double)k * (double)k);
        sum += term;
    }
    return sum;
}

int main(int argc, char **argv)
{
    const char *path = NULL;
    size_t N = 65536;
    double line_lo = 0.0, line_hi = 0.0;
    for (int i = 1; i < argc; i++) {
        if (!strcmp(argv[i], "-l") && i + 2 < argc) {
            line_lo = strtod(argv[++i], NULL);
            line_hi = strtod(argv[++i], NULL);
        } else if (!strcmp(argv[i], "-n") && i + 1 < argc)
            N = strtoul(argv[++i], NULL, 10);
        else
            path = argv[i];
    }
    if (!path || N < 1024 || (N & (N - 1))) {
        fprintf(stderr, "usage: flacspec FILE [-n FFTSIZE] [-l LO HI]   (FFTSIZE: power of two >= 1024)\n");
        return 2;
    }

    char err[512];
    fu_format info;
    fu_input *in = fu_input_open(path, &info, err, sizeof err);
    if (!in) {
        fprintf(stderr, "flacspec: %s\n", err);
        return 1;
    }
    const uint32_t nch = info.channels;
    const size_t nbin = N / 2 + 1;
    double *win = malloc(N * sizeof *win), *psd = calloc(nbin, sizeof *psd);
    double *frame = fftw_malloc(N * sizeof *frame);
    fftw_complex *spec = fftw_malloc(nbin * sizeof *spec);
    double *ring[FU_MAX_CH] = { 0 };
    for (uint32_t ch = 0; ch < nch; ch++)
        ring[ch] = calloc(N, sizeof(double));
    if (!win || !psd || !frame || !spec || !ring[nch - 1]) {
        fprintf(stderr, "flacspec: out of memory\n");
        return 1;
    }
    fftw_plan plan = fftw_plan_dft_r2c_1d((int)N, frame, spec, FFTW_ESTIMATE);

    double wsq = 0.0;
    const double i0b = bessel_i0(23.3);
    for (size_t i = 0; i < N; i++) {
        const double u = (2.0 * (double)i - (double)(N - 1)) / (double)(N - 1);
        win[i] = bessel_i0(23.3 * sqrt(1.0 - u * u)) / i0b;
        wsq += win[i] * win[i];
    }

    size_t fill = 0, nframes = 0;
    uint64_t total = 0;
    for (;;) {
        const double *planes[FU_MAX_CH];
        const long n = fu_input_read(in, planes);
        if (n < 0) {
            fprintf(stderr, "flacspec: %s: decode error\n", path);
            return 1;
        }
        if (n == 0)
            break;
        total += (uint64_t)n;
        for (long done = 0; done < n;) {
            const size_t take = (size_t)(n - done) < N - fill ? (size_t)(n - done) : N - fill;
            for (uint32_t ch = 0; ch < nch; ch++)
                memcpy(ring[ch] + fill, planes[ch] + done, take * sizeof(double));
            fill += take;
            done += (long)take;
            if (fill < N)
                continue;
            for (uint32_t ch = 0; ch < nch; ch++) {
                for (size_t i = 0; i < N; i++)
                    frame[i] = ring[ch][i] * win[i];
                fftw_execute(plan);
                for (size_t k = 0; k < nbin; k++)
                    psd[k] += spec[k][0] * spec[k][0] + spec[k][1] * spec[k][1];
                memmove(ring[ch], ring[ch] + N / 8, (N - N / 8) * sizeof(double));
            }
            fill = N - N / 8;
            nframes++;
        }
    }
    fu_input_close(in);
    if (nframes == 0) {
        fprintf(stderr, "flacspec: %s: shorter than one %zu-sample frame; use a smaller -n\n", path, N);
        return 1;
    }
    /* One-sided power per bin, averaged over frames and channels. */
    for (size_t k = 0; k < nbin; k++)
        psd[k] *= ((k == 0 || k == N / 2) ? 1.0 : 2.0) / ((double)N * wsq * (double)nframes * nch);

    static const double edges[] = { 0, 20, 60, 120, 250, 1000, 4000, 10000, 16000, 20000,
                                    22050, 25000, 30000, 40000, 60000, 96000, 192000, 384000 };
    const double hz = (double)info.rate / (double)N, nyq = 0.5 * info.rate;
    printf("flacspec: %s\n", path);
    printf("  %u Hz, %u bit, %u ch, %llu frames; FFT %zu (%.3f Hz per bin), %zu frames averaged\n",
           info.rate, info.bps, nch, (unsigned long long)total, N, hz, nframes);
    printf("  %9s %9s | %10s | %12s %10s\n", "from Hz", "to Hz", "power dBFS", "top line Hz", "line dBFS");
    double all = 0.0;
    for (size_t b = 0; b + 1 < sizeof edges / sizeof *edges && edges[b] < nyq; b++) {
        const double f0 = edges[b], f1 = edges[b + 1] < nyq ? edges[b + 1] : nyq;
        const size_t k0 = (size_t)ceil(f0 / hz), k1 = (size_t)ceil(f1 / hz) + (f1 >= nyq);
        double p = 0.0, top = -1.0;
        size_t ktop = k0;
        for (size_t k = k0; k < k1 && k < nbin; k++) {
            p += psd[k];
            if (psd[k] > top) {
                top = psd[k];
                ktop = k;
            }
        }
        double line = 0.0;
        for (size_t k = ktop > 10 ? ktop - 10 : 0; k <= ktop + 10 && k < nbin; k++)
            line += psd[k];
        all += p;
        printf("  %9.0f %9.0f | %10.2f | %12.1f %10.2f\n", f0, f1, 10.0 * log10(p / 0.5 + 1e-300),
               (double)ktop * hz, 10.0 * log10(line / 0.5 + 1e-300));
    }
    printf("  %9s %9s | %10.2f |\n", "all", "", 10.0 * log10(all / 0.5 + 1e-300));

    if (line_hi > line_lo) {
        /* Strongest lines first; each removes the 21 bins of its main lobe. */
        const size_t k0 = (size_t)ceil(line_lo / hz), k1 = (size_t)floor(line_hi / hz) < nbin ? (size_t)floor(line_hi / hz) : nbin - 1;
        double total = 0.0, lines = 0.0;
        for (size_t k = k0; k <= k1; k++)
            total += psd[k];
        printf("\n  lines between %.0f and %.0f Hz (band power %.2f dBFS)\n  %12s %10s\n", line_lo, line_hi,
               10.0 * log10(total / 0.5 + 1e-300), "line Hz", "dBFS");
        for (int i = 0; i < 12; i++) {
            size_t kt = k0;
            for (size_t k = k0; k <= k1; k++)
                if (psd[k] > psd[kt])
                    kt = k;
            double p = 0.0;
            for (size_t k = kt > k0 + 10 ? kt - 10 : k0; k <= kt + 10 && k <= k1; k++) {
                p += psd[k];
                psd[k] = 0.0;
            }
            if (p <= 0.0)
                break;
            lines += p;
            printf("  %12.1f %10.2f\n", (double)kt * hz, 10.0 * log10(p / 0.5 + 1e-300));
        }
        printf("  %12s %10.2f   (everything else in the range)\n", "rest", 10.0 * log10((total - lines) / 0.5 + 1e-300));
    }

    fftw_destroy_plan(plan);
    fftw_free(frame);
    fftw_free(spec);
    free(win);
    free(psd);
    for (uint32_t ch = 0; ch < nch; ch++)
        free(ring[ch]);
    return 0;
}
