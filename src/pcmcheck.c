/*
 * pcmcheck.c -- gate for fu_pcm: the AVX2 forms must equal the plain forms and
 *               must not be slower than them.
 *
 * Purpose : (1) equality at 32, 24 and 16 bits on synthetic data that straddles
 *               every boundary of the clip rule (exact ties, both saturation
 *               edges, tiny values) plus 2^22 pseudo-random samples at gains
 *               that force clipping;
 *           (2) timing, 5 interleaved rounds over 2^24 samples.
 * Usage   : pcmcheck   (no arguments; exit 0 = pass)
 */
#include "fu_pcm.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

static double now_ms(void)
{
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return 1e3 * (double)ts.tv_sec + 1e-6 * (double)ts.tv_nsec;
}

int main(void)
{
    enum { N = 1 << 24, NCHECK = 1 << 22, ROUNDS = 5 };
    double  *x = malloc(N * sizeof *x);
    int32_t *a = malloc(N * sizeof *a), *b = malloc(N * sizeof *b);
    if (!x || !a || !b)
        return 2;

    /* Boundary offsets around +-full scale, in units of one output step. */
    static const double off[] = { -1.5, -1.0, -0.75, -0.5, -0.25, 0.0, 0.25, 0.5, 1.0, 1e3 };
    static const unsigned widths[] = { 32, 24, 16 };
    size_t nedge = 0;
    for (size_t w = 0; w < 3; w++) {
        const double fs = (double)(1ull << (widths[w] - 1));
        for (size_t i = 0; i < sizeof off / sizeof *off; i++) {
            x[nedge++] = (fs + off[i]) / fs;
            x[nedge++] = (-fs - off[i]) / fs;
        }
        x[nedge++] = 0.5 / fs;
        x[nedge++] = -0.5 / fs;
        x[nedge++] = 1.5 / fs;
        x[nedge++] = 2.5 / fs;
    }
    uint64_t rng = 0xD1B54A32D192ED03ull;
    for (size_t i = nedge; i < N; i++) {
        rng ^= rng << 13;
        rng ^= rng >> 7;
        rng ^= rng << 17;
        x[i] = (double)(rng >> 11) * 0x1p-52 - 1.0;
    }

    int fail = 0;
    static const double gains[] = { 1.0, 0.5, 1.0000001, 1.7, 3.0 };
    for (size_t w = 0; w < 3; w++) {
        for (size_t g = 0; g < sizeof gains / sizeof *gains; g++) {
            for (size_t n = NCHECK - 5; n <= NCHECK; n++) {
                const size_t ca = fu_pcm_to_int(x, a, n, gains[g], widths[w]);
                const size_t cb = fu_pcm_to_int_ref(x, b, n, gains[g], widths[w]);
                if (ca != cb || memcmp(a, b, n * sizeof *a) != 0) {
                    printf("FAIL to_int: %u bit gain %g n %zu clipped %zu vs %zu\n", widths[w],
                           gains[g], n, ca, cb);
                    fail = 1;
                }
            }
        }
        printf("to_int %2u bit: %zu samples, clipped %zu at gain 1, %zu at gain 3, identical to "
               "the plain form\n", widths[w], (size_t)NCHECK,
               fu_pcm_to_int(x, a, NCHECK, 1.0, widths[w]),
               fu_pcm_to_int(x, a, NCHECK, 3.0, widths[w]));
    }
    for (size_t n = NCHECK - 20; n <= NCHECK; n++) {
        if (fu_pcm_peak(x, n) != fu_pcm_peak_ref(x, n)) {
            printf("FAIL peak: n %zu\n", n);
            fail = 1;
        }
    }
    printf("peak: %.17g, identical to the plain form\n", fu_pcm_peak(x, NCHECK));

    double best[4] = { 1e30, 1e30, 1e30, 1e30 }, worst[4] = { 0 };
    double sink = 0.0;
    for (int r = 0; r < ROUNDS; r++) {
        for (int v = 0; v < 4; v++) {
            const double t0 = now_ms();
            switch (v) {
            case 0: sink += (double)fu_pcm_to_int_ref(x, a, N, 1.0, 32); break;
            case 1: sink += (double)fu_pcm_to_int(x, a, N, 1.0, 32); break;
            case 2: sink += fu_pcm_peak_ref(x, N); break;
            default: sink += fu_pcm_peak(x, N); break;
            }
            const double ms = now_ms() - t0;
            best[v] = ms < best[v] ? ms : best[v];
            worst[v] = ms > worst[v] ? ms : worst[v];
        }
    }
    static const char *const name[4] = { "to_int plain", "to_int AVX2", "peak plain", "peak AVX2" };
    printf("\n%-14s %10s %10s %14s   (%d samples, %d interleaved rounds)\n", "kernel", "best ms",
           "worst ms", "Msamples/s", (int)N, (int)ROUNDS);
    for (int v = 0; v < 4; v++)
        printf("%-14s %10.2f %10.2f %14.0f\n", name[v], best[v], worst[v], N / best[v] / 1e3);
    printf("timing checksum %.17g\n", sink);
    if (best[1] > best[0] || best[3] > best[2]) {
        printf("FAIL: an AVX2 form is slower than its plain form\n");
        fail = 1;
    }
    free(x);
    free(a);
    free(b);
    printf("%s\n", fail ? "pcmcheck: FAIL" : "pcmcheck: PASS");
    return fail;
}
