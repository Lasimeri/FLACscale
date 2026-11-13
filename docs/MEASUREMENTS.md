# Measurements

Every number in the README and in DESIGN.md comes from here.

- **Machine:** AMD Ryzen 7 5800X, Linux 7.3.1 (CachyOS), gcc 16.2.1, flags as in
  the `Makefile`.
- **Libraries:** libFLAC 1.5.0, mpg123 1.33.7, libsndfile 1.2.2, FFTW 3.3.11,
  libsoxr 0.1.3.
- **Date:** 2026-11-13, unless a row says otherwise.
- **Reproduce:** `make check` runs sections 1 to 5 on synthetic signals.
  Section 6 needs an unclipped audio file of your own.

Speed figures move by about 10 % between runs on this machine; best and worst
of the interleaved rounds are both given where the tool reports them.

## 1. PCM kernels (`out/pcmcheck`)

AVX2 forms against the plain scalar forms: identical results at 32, 24 and 16
bits (4,194,304 samples each, including every boundary of the clip rule), then
timing over 16,777,216 samples, 5 interleaved rounds.

| Kernel | Best ms | Worst ms | Msamples/s |
|---|---|---|---|
| `to_int` plain | 105.06 | 131.98 | 160 |
| `to_int` AVX2 | 7.51 | 8.40 | 2233 |
| peak plain | 39.94 | 47.19 | 420 |
| peak AVX2 | 4.04 | 17.13 | 4154 |

Best against best: 14.0x for the conversion, 9.9x for the peak search.

## 2. Resampler filters (`out/rsbench`, table 1)

Response computed from the prototype taps (DTFT), frequencies at the stage
input.

| Design | L | Taps | Pass Hz | Stop Hz | Passband ripple dB p-p | Stopband peak dB |
|---|---|---|---|---|---|---|
| single 640/147 | 640 | 256 | 19999 | 22050 | 3.653e-08 | -170.7 |
| stage A, 2/1 steep | 2 | 256 | 19999 | 22050 | 3.655e-08 | -170.6 |
| stage B, 320/147 relaxed | 320 | 32 | 22050 | 66150 | 9.153e-11 | -221.0 |

The Kaiser design formula predicts 178.9 dB for the steep filter. The measured
stopband peak is -170.7 dB. The measured number is the one quoted.

## 3. Resampler correctness and accuracy (`out/rsbench`, tables 2 and 3)

On 10 s of synthetic stereo noise (441,000 frames in, 1,920,000 out):

- AVX2 kernels against the scalar formula: largest difference 1.665e-15.
- Ragged chunk sizes against 4096-frame chunks: 0 (bit-identical).

Synthetic sine, amplitude 0.5, 1 s, fitted over 0.1 s to 0.9 s at 192 kHz:

| Tone Hz | FLACscale gain dB | FLACscale residual dB | libsoxr VHQ gain dB | libsoxr VHQ residual dB |
|---|---|---|---|---|
| 30 | -8.57e-10 | -209.6 | -6.34e-09 | -187.6 |
| 1000 | -8.36e-10 | -209.3 | -3.64e-09 | -187.1 |
| 10000 | -9.80e-10 | -206.7 | -1.11e-09 | -191.3 |
| 19000 | -5.59e-09 | -190.2 | -3.00e-03 | -184.5 |
| 20000 | -2.95e-08 | -188.3 | -6.22e-03 | -191.8 |

Time offset: 0.000 ns for both at every tone.

## 4. Resampler speed (`out/rsbench`, table 4)

10 s stereo, 44100 -> 192000, 1 thread, 3 interleaved rounds.

| Variant | Best ms | Worst ms | x real time |
|---|---|---|---|
| scalar formula, single stage 640/147 x 256 | 1816.8 | 1901.0 | 5.5 |
| AVX2 phase-order, single stage | 99.4 | 111.0 | 100.6 |
| AVX2 visit-order, single stage | 73.7 | 82.0 | 135.7 |
| AVX2 two-stage chain (the one used) | 42.2 | 43.2 | 236.8 |
| libsoxr VHQ (28-bit recipe) | 30.1 | 53.7 | 331.8 |
| libsoxr 32-bit recipe | 46.4 | 58.4 | 215.5 |

**Not overcome:** libsoxr VHQ is 1.40x faster than FLACscale's chain in this
run (1.37x in a run on 2026-11-12: 27.5 ms against 37.7 ms). FLACscale is the
more accurate of the two at 19 kHz and 20 kHz (table above) and 1.10x faster
than libsoxr's 32-bit recipe.

## 5. Band extension on a synthetic signal (`make check`, last three lines)

Input from `out/tonegen`: sines at 9000 Hz and 14000 Hz (-20 dBFS each) plus
white noise low-passed at 15500 Hz (-50 dBFS RMS), 10 s, 44100 Hz, 16 bit.

- Detected wall: 15750 Hz (a 49.5 dB drop).
- Folded band: 15750 to 25000 Hz, from the band below 15500 Hz.
- Total power between 25 kHz and 96 kHz: -46.99 dBFS, in a file whose total
  power is -16.98 dBFS (30.01 dB below it: the guard at work).

Strongest lines between 25 kHz and 96 kHz (`flacspec -l 25000 96000`):

| Line Hz | dBFS | Is |
|---|---|---|
| 27999.0 | -51.65 | 2 x 14000 |
| 42000.0 | -55.17 | 3 x 14000 |
| 27000.0 | -55.22 | 3 x 9000 |
| 36000.0 | -57.71 | 4 x 9000 |
| 56001.0 | -57.66 | 4 x 14000 |
| 45000.0 | -59.64 | 5 x 9000 |
| 69999.0 | -59.60 | 5 x 14000 |
| 84000.0 | -61.19 | 6 x 14000 |
| 54000.0 | -61.22 | 6 x 9000 |
| everything else in the range | -95.30 | residue |

The nine strongest lines are exactly the nine integer multiples of the two
tones that fall between 25 kHz and 96 kHz. Everything else in that range sums
to -95.30 dBFS. (Line frequencies are bin centres, 2.93 Hz apart.)

## 6. Clipped-peak repair against a known truth (`out/declipcheck FILE`)

Method: take an unclipped recording (248.041 s, 44100 Hz, stereo, checked with
`clipscan`), clip it artificially so a chosen share of samples is pinned,
repair it, and score the pinned samples only against the original:
`SDR = 10 * log10(sum(truth^2) / sum((truth - estimate)^2))`.

Default settings (noise floor -80 dB, 1 pass):

| Pinned share | Ceiling dBFS | SDR clipped dB | SDR repaired dB | Gain dB | Raised | Held | Refused |
|---|---|---|---|---|---|---|---|
| 0.3 % | -3.70 | 22.57 | 31.30 | +8.73 | 61,290 | 4,322 | 0 |
| 0.9 % | -4.61 | 19.40 | 26.00 | +6.60 | 175,053 | 21,773 | 0 |
| 3.0 % | -6.21 | 15.40 | 19.53 | +4.12 | 529,867 | 126,312 | 0 |

The same input with other settings, gain in dB at 0.3 / 0.9 / 3.0 %:

| Setting | 0.3 % | 0.9 % | 3.0 % | Highest repaired peak dBFS |
|---|---|---|---|---|
| Opus guards (floor -40 dB), 1 pass | +7.90 | +6.02 | +3.78 | +0.07 |
| Opus guards, 2 passes | +7.88 | +6.06 | +3.88 | +0.20 |
| floor -80 dB, 1 pass (default) | +8.73 | +6.60 | +4.12 | +0.07 |
| floor -80 dB, 2 passes | +8.63 | +6.63 | +4.24 | +1.76 |
| floor -80 dB, lag window 0.002, 1 pass | +8.70 | +6.56 | +4.09 | +0.07 |

A second pass changes the result by -0.10 to +0.12 dB and produces overshoot
(+1.76 dBFS), so the default is one pass. The -80 dB floor is worth +0.34 to
+0.83 dB over Opus's -40 dB. The repair gets weaker as clipping gets heavier: at 3 % of
samples pinned it recovers 4.12 dB, not 8.73 dB.

This is one recording. It is not a claim about all music.

## 7. Input paths (2026-11-12)

- **MP3:** decoded output against ffmpeg's decode of the same file, three
  files: residual -119 to -121 dB at lag 0 (`out/flacdiff`).
- **WAV:** against the same audio as FLAC: residual exactly zero.
- **Resample-only output** (`-D -X`): power above the source's Nyquist sits at
  the 32-bit quantisation floor, -195.5 dBFS (`out/flacspec`).

## 8. Whole path on real files

Thirteen music files from a private library (not distributed), full default
path, 2026-11-12 and 2026-11-13. Titles are omitted except for the worked
example in the README.

| Source | Wall Hz | Repair: raised / held | Gain dB |
|---|---|---|---|
| FLAC 44.1 kHz / 16 bit | 20250 | 99,193 / 108,889 | -3.1708 |
| FLAC 44.1 kHz / 16 bit | none | 97,344 / 364,668 | -6.2882 |
| FLAC 44.1 kHz / 16 bit | none | 26,030 / 1,172 | -1.7294 |
| FLAC 44.1 kHz / 16 bit | 17250 | 347 / 62 | -0.6419 |
| FLAC 44.1 kHz / 16 bit | none | not clipped by the detector's rule | -0.0362 |
| FLAC 44.1 kHz / 16 bit | none | not clipped | 0 |
| FLAC 48 kHz / 24 bit | 22500 | 1,457 / 5,155 | -0.7678 |
| WAV 44.1 kHz / 16 bit | none | not clipped | -0.0873 |
| MP3 | 16000 | lossy: skipped | -2.4547 |
| MP3 | 15750 | lossy: skipped | -1.8466 |
| MP3 | 14500 | lossy: skipped | -0.2644 |
| MP3 | 15250 | lossy: skipped | -0.2801 |
| MP3 | 15250 | lossy: skipped | 0 |

- Every output passed `flac -t` and reported 0 clipped samples.
- Speed of the whole path (scan, measure, write): 23x to 33x real time.
- Heavily clipped masters lose the most level, because restored peaks need
  headroom: the -6.2882 dB file had a highest repaired peak of +6.287 dBFS.

## 9. A second machine: AMD EPYC 7742 (64 cores, 128 threads)

Built from a fresh clone of this repository, same compiler and library
versions, Linux 7.3.3. One FLAC file, 588.516 s, 44.1 kHz / 16 bit, default
path. The machine was also running an idle inference server.

**Same bits.** The output's audio MD5 equals the one produced on the Ryzen
5800X, for this file and for the `make check` tones.

**One conversion.** The signal path is single-threaded; `-t` sets the FLAC
encoder's threads only.

| `-t` | Total s | x real time | Encode s |
|---|---|---|---|
| 1 | 26.76 | 22.0 | 6.54 |
| 4 (default) | 20.60 | 28.6 | 0.36 |
| 16 | 20.97 | 28.1 | 0.69 |
| 64 | 20.85 | 28.2 | 0.63 |

For comparison, the 5800X converts the same file in 19.69 s (29.9x). libFLAC
1.5.0 refused 96 and 128 threads and encoded on one (6.53 s); `flacscale` now
steps down to the largest count the library accepts.

**Many conversions at once** (`-t 1` each, output discarded):

| Simultaneous | Wall s | Each instance s | Audio converted | Throughput |
|---|---|---|---|---|
| 1 | 26.76 | 26.76 | 588.5 s | 22.0x real time |
| 64 | 38.37 | 36.12 to 38.34 | 37,665 s | 981.7x real time |
| 128 | 68.83 | 67.32 to 68.80 | 75,330 s | 1094.4x real time |

Filling the second hardware thread of every core adds 11.5 % throughput. Each
instance slows from 26.76 s alone to about 37 s with all 64 cores busy.

Gates on that machine: `pcmcheck` PASS (conversion 10.0x, peak 7.0x over the
plain forms); resampler chain 55.7 ms per 10 s stereo (179.6x), libsoxr VHQ
34.0 ms (1.64x faster there).
