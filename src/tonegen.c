/*
 * tonegen.c -- SYNTHETIC test signal for the band-extension gate (not music).
 *
 * Purpose : a signal whose correct extension is known in advance:
 *             sines at 9000 Hz and 14000 Hz, -20 dBFS each
 *           + white noise low-passed at 15500 Hz, -50 dBFS RMS
 *           10 s, 44100 Hz, 16 bit, stereo (both channels identical), WAV.
 *           The converter must report a wall at 15750 Hz, and above 25 kHz the
 *           output must hold integer harmonics of 9000 and 14000 Hz only:
 *           27000, 28000, 36000, 42000, 45000, 54000, 56000, 70000, 84000 Hz.
 * Method  : xorshift64 noise (fixed seed, so the file is the same on every
 *           machine up to libm rounding), 511-tap Blackman-windowed sinc.
 * Usage   : tonegen OUT.wav
 */
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

enum { RATE = 44100, SECONDS = 10, FRAMES = RATE * SECONDS, TAPS = 511, CHANNELS = 2, BITS = 16 };

static const double TONE_HZ[2] = { 9000.0, 14000.0 };
static const double TONE_DBFS  = -20.0;
static const double NOISE_DBFS = -50.0;
static const double NOISE_TOP  = 15500.0;

/* Uniform white noise in [-1, 1): n values from xorshift64. */
static void fill_noise(double *x, size_t n)
{
    uint64_t r = 0x2545F4914F6CDD1Dull;
    for (size_t i = 0; i < n; i++) {
        r ^= r << 13;
        r ^= r >> 7;
        r ^= r << 17;
        x[i] = (double)(r >> 11) * 0x1p-52 - 1.0;
    }
}

/* Low-pass prototype: Blackman-windowed sinc, cutoff fc in cycles per sample. */
static void design_lowpass(double *h, int taps, double fc)
{
    for (int i = 0; i < taps; i++) {
        const double t = i - (taps - 1) / 2.0;
        const double w = 0.42 - 0.5 * cos(2 * M_PI * i / (taps - 1))
                       + 0.08 * cos(4 * M_PI * i / (taps - 1));
        h[i] = (t == 0.0 ? 2 * fc : sin(2 * M_PI * fc * t) / (M_PI * t)) * w;
    }
}

/* 16-bit PCM WAV header for `frames` interleaved stereo frames. */
static int write_header(FILE *f, uint32_t frames)
{
    const uint32_t bytes = frames * CHANNELS * (BITS / 8), riff = 36 + bytes, fmtlen = 16;
    const uint32_t rate = RATE, byterate = RATE * CHANNELS * (BITS / 8);
    const uint16_t pcmtag = 1, ch = CHANNELS, align = CHANNELS * (BITS / 8), bits = BITS;
    size_t ok = fwrite("RIFF", 1, 4, f);
    ok += fwrite(&riff, 4, 1, f);
    ok += fwrite("WAVEfmt ", 1, 8, f);
    ok += fwrite(&fmtlen, 4, 1, f);
    ok += fwrite(&pcmtag, 2, 1, f);
    ok += fwrite(&ch, 2, 1, f);
    ok += fwrite(&rate, 4, 1, f);
    ok += fwrite(&byterate, 4, 1, f);
    ok += fwrite(&align, 2, 1, f);
    ok += fwrite(&bits, 2, 1, f);
    ok += fwrite("data", 1, 4, f);
    ok += fwrite(&bytes, 4, 1, f);
    return ok == 4 + 1 + 8 + 7 + 4 + 1;
}

int main(int argc, char **argv)
{
    if (argc != 2) {
        fprintf(stderr, "usage: tonegen OUT.wav\n");
        return 2;
    }
    double  *noise = calloc(FRAMES + TAPS, sizeof *noise);
    double  *h = malloc(TAPS * sizeof *h);
    int16_t *pcm = malloc((size_t)CHANNELS * FRAMES * sizeof *pcm);
    if (!noise || !h || !pcm) {
        fprintf(stderr, "tonegen: out of memory\n");
        return 1;
    }
    fill_noise(noise, FRAMES + TAPS);
    design_lowpass(h, TAPS, NOISE_TOP / RATE);

    /* Uniform noise in [-1, 1) has RMS 0.577; the low-pass keeps 15.5 / 22.05 of
     * the band, which leaves RMS 0.484. The sqrt(0.5) puts the noise on the same
     * dBFS scale as the tones (0 dBFS = a full-scale sine). */
    const double noise_gain = pow(10.0, NOISE_DBFS / 20.0) / 0.484 * sqrt(0.5);
    const double tone_gain = pow(10.0, TONE_DBFS / 20.0);
    for (int n = 0; n < FRAMES; n++) {
        double v = 0.0;
        for (int i = 0; i < TAPS; i++)
            v += h[i] * noise[n + i];
        const double x = noise_gain * v + tone_gain * sin(2 * M_PI * TONE_HZ[0] * n / RATE)
                       + tone_gain * sin(2 * M_PI * TONE_HZ[1] * n / RATE);
        pcm[2 * n] = pcm[2 * n + 1] = (int16_t)lrint(x * 32767.0);
    }

    FILE *f = fopen(argv[1], "wb");
    if (!f) {
        fprintf(stderr, "tonegen: cannot write %s\n", argv[1]);
        return 1;
    }
    const int ok = write_header(f, FRAMES)
                && fwrite(pcm, sizeof *pcm, (size_t)CHANNELS * FRAMES, f) == (size_t)CHANNELS * FRAMES;
    if (fclose(f) != 0 || !ok) {
        fprintf(stderr, "tonegen: write failed: %s\n", argv[1]);
        return 1;
    }
    free(noise);
    free(h);
    free(pcm);
    return 0;
}
