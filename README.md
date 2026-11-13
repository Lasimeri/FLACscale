# FLACscale

Upscales audio files to 192 kHz / 32-bit FLAC. FLAC, WAV, AIFF and MP3 in, FLAC
out. Written in C with hand-written AVX2 + FMA3 kernels.

It does three things beyond a plain sample-rate conversion, and the first two
are a native rework of the loss-concealment machinery of the Opus decoder:

1. **Repairs clipped peaks** with Opus's 24th-order linear predictor, solved
   from both sides of the clipped run.
2. **Fills the band the source never had** (above an MP3's 16 kHz low-pass, for
   example) by Opus-style band folding up to 25 kHz, and above 25 kHz by
   integer harmonics of the tones it finds, and nothing else.
3. **Never clips the result:** the true peak is measured first, and the file is
   attenuated only as much as that peak requires.

## Read this before you expect miracles

- **It cannot recover information.** A 44.1 kHz file holds 0 to 22.05 kHz
  completely. Everything FLACscale writes above the source's bandwidth is
  synthesis: built from the content below it, bounded in level, and labelled
  as synthesis in the report.
- **More bits do not add resolution to the source.** The 32-bit output holds
  the result of double-precision processing without a coarser rounding step.
- **What is real:** clipped peaks reconstructed (scored against a known truth:
  +8.73 dB at light clipping, +4.12 dB at heavy), a resampling filter measured
  at -170.6 dB in the stopband, and exact time alignment with the source.
- **Opus forward error correction (FEC) is not used.** It needs redundancy made
  by an encoder, which a finished file does not contain. The part of Opus that
  is reused is the decoder-side concealment. Details:
  [docs/DESIGN.md](docs/DESIGN.md).

## Requirements

- x86-64 CPU with AVX2, FMA3 and BMI2 (Intel Haswell or AMD Zen and newer).
- gcc with C23 support (`-std=gnu2x`) and make.
- libFLAC 1.5.0 or newer (its multithreaded encoder is used), libmpg123,
  libsndfile, FFTW 3, libsoxr (only for the benchmark `rsbench`).

Built and tested on one system only: CachyOS (Arch Linux), gcc 16.2.1,
libFLAC 1.5.0, mpg123 1.33.7, libsndfile 1.2.2, FFTW 3.3.11, libsoxr 0.1.3.

    sudo pacman -S --needed gcc make flac mpg123 libsndfile fftw libsoxr

## Build and check

    make          # every tool lands in out/
    make check    # gates on synthetic signals, about 30 s

`make check` exits non-zero when a kernel differs from its plain form or is
slower than it. Its tables are explained in
[docs/MEASUREMENTS.md](docs/MEASUREMENTS.md).

## Use

One file, with a read-only checksummed copy of the input and a verified output:

    scripts/upscale "/path/to/song.flac"
    scripts/upscale "/path/to/song.flac" myfolder -Q     # options after the folder name go to flacscale

The result is `out/<song>/<song>.192k32.flac`, with the report in
`out/<song>/run.log`. An MP3 input also gets `<song>.flac`, a plain 24-bit
conversion at the source rate.

The converter directly:

    out/flacscale IN [-o OUT.flac] [-r RATE] [-g auto|DB] [-t THREADS] [-n] [-D] [-X] [-Q]

| Option | Effect | Default |
|---|---|---|
| `-o OUT.flac` | output file | input's name + `.192k32.flac`, in the current directory |
| `-r RATE` | output sample rate in Hz | 192000 |
| `-g auto` or `-g DB` | output gain; a number is a fixed gain, e.g. `-g -3` | auto: only as much attenuation as the true peak needs |
| `-t THREADS` | FLAC encoder threads | 4 |
| `-n` | lossy input: skip the extra 24-bit FLAC copy | off |
| `-D` | skip the clipped-peak repair | off |
| `-X` | skip the band extension (plain resample only) | off |
| `-Q` | restoration mode: first back-calculate the precision that rounding to the source's bit depth removed (integer PCM sources only; see below) | off |

## What a run reports

A 310 s MP3 (Sytricka, "Running In The 90s (Vaporwave Remix)"), with the
paths and the file name shortened:

```
flacscale: Sytricka - Running In The 90s (Vaporwave Remix).mp3
  input      MP3   44100 Hz  lossy   2 ch  13676544 frames  310.126 s  8.6 MiB
  converted  FLAC  44100 Hz  24 bit  2 ch  13676544 frames  61.7 MiB  peak -0.123 dBFS  gain +0.0000 dB  clipped 0
  bandwidth  the source's content ends at 15250 Hz (a 37.4 dB wall)
  extension  folded 15250 .. 25000 Hz from the band below 15000 Hz; above 25000 Hz harmonics only: 4577 tonal peaks -> 17506 partials (56.4 per second)
  output     FLAC 192000 Hz  32 bit  2 ch  59544138 frames  253.8 MiB
  true peak  -0.010 dBFS at 192000 Hz   gain +0.0000 dB (auto)   clipped samples 0
  time       scan+measure 4.91 s | decode 0.26 s  dsp 4.14 s  encode 0.37 s | total 9.72 s  (31.9x real time)
```

Reading it: the MP3 encoder cut everything above 15250 Hz; the band from there
to 25 kHz was synthesised by folding; above 25 kHz only harmonics were added;
the level was left alone because the true peak fits. The audio itself is not
distributed here.

## Restoration mode (`-Q`)

Off by default; without it the converter behaves exactly as before. With it,
the 16-bit samples are treated as what they are: statements that the true
value lay within half a step of each of them. A Wiener filter removes rounding
noise where it stands alone, and every sample is then forced back inside its
interval, so no sample ever moves more than half a step from the file.

Measured against a 96 kHz / 32-bit master cut down to 44.1 kHz / 16 bit: the
result is 0.42 dB closer to the truth as mastered, 2.82 dB closer in a passage
60 dB down, and 0.22 dB closer after the whole upscale. It is a small, real
gain that grows as the music gets quieter. A 24-bit source has nothing to
recover (samples moved 0.007 steps RMS on a real one). Lossy and float sources
are skipped.

## Signal path

1. **Scan:** per channel, is it clipped and at what level; where does the
   source's real bandwidth end.
2. **Repair** (source rate): clusters of clipped samples are re-solved under
   the Opus LPC model with both sides known. Every repaired sample stays at or
   beyond its clipped value.
3. **Resample:** two polyphase stages in double precision, time-aligned to the
   input.
4. **Extend** (output rate): folding from the source's bandwidth to 25 kHz;
   above 25 kHz only integer harmonics of tonal peaks, capped per frame at
   -30 dB of the real signal's power.
5. **Level and encode:** smallest attenuation that avoids clipping, 32-bit
   FLAC, tags and cover art carried over.

The full description, with a diagram, is in [docs/DESIGN.md](docs/DESIGN.md).

## Measured

On a Ryzen 7 5800X, from [docs/MEASUREMENTS.md](docs/MEASUREMENTS.md):

| What | Result |
|---|---|
| Whole path on music files | 23x to 33x real time |
| Resampler alone, 44.1 -> 192 kHz stereo | 236.8x real time |
| Resampler stopband | -170.6 dB |
| Resampler error on a 20 kHz tone | gain -2.95e-08 dB, residual -188.3 dB, offset 0.000 ns |
| Clip repair, 0.3 % / 0.9 % / 3 % of samples clipped | +8.73 / +6.60 / +4.12 dB SDR |
| Synthetic tones: lines above 25 kHz | exactly the 9 expected harmonics; all else -95.30 dBFS |

Not overcome: libsoxr's VHQ recipe is 1.40x faster than FLACscale's resampler.

## Limits

- Repaired files are quieter, because restored peaks need headroom. A heavily
  clipped master can lose 6 dB.
- Clipping that went through a lossy encoder is not repaired.
- A channel with fewer than 32 samples at its extreme counts as not clipped
  and is left alone.
- The bandwidth detector needs a drop of at least 20 dB. A gentler roll-off is
  treated as content that reaches Nyquist.
- Files only. There is no real-time player yet; see
  [docs/ROADMAP.md](docs/ROADMAP.md).

## Tools

All in `out/` after `make`.

| Tool | Purpose |
|---|---|
| `flacscale` | the converter |
| `flacspec FILE [-n FFTSIZE] [-l LO HI]` | band powers and the strongest spectral lines of a file |
| `flacdiff A B` | null test between two files: offset, gain, residual |
| `clipscan FILE [-t TOL]` | how many samples sit pinned at the extremes, and in what runs |
| `declipcheck FILE` | scores the repair against a known truth (needs an unclipped file) |
| `dequantcheck FILE` | scores the dequantizer against a known truth (needs a file finer than 16 bits) |
| `rsbench`, `pcmcheck` | gates for the resampler and the PCM kernels |
| `tonegen OUT.wav` | the synthetic test signal used by `make check` |

## Layout

    src/        C sources (fu_* = engine modules, the rest = tools)
    scripts/    upscale: copy, convert, verify
    docs/       DESIGN.md, MEASUREMENTS.md, ROADMAP.md
    input/      local: read-only copies of your files + SHA256SUMS (not in git)
    out/        local: built tools, converted files, logs (not in git)

## Licence

MIT, see [LICENSE](LICENSE). The LPC routines ported from libopus keep its
BSD-3-Clause notice, and the linked libraries have their own licences (FFTW is
GPL, which matters if you redistribute binaries): see [NOTICE.md](NOTICE.md).

Contributions: see [CONTRIBUTING.md](CONTRIBUTING.md).
