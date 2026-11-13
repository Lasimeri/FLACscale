# Design

How FLACscale is built, stage by stage, and which part of each stage comes from
the Opus decoder. Every number here is measured; the commands and tables are in
[MEASUREMENTS.md](MEASUREMENTS.md).

```
╔══════════════════════════════════════════════════════════════════════════════════════╗
║ desktop: Ryzen 7 5800X (AVX2 + FMA3)   out/flacscale        23-33x real time         ║
║                                                                                      ║
║ ┌──────────────────────────────────────────────────────────────────────────────────┐ ║
║ │ fu_input   FLAC (libFLAC) | WAV, AIFF (libsndfile) | MP3 (libmpg123, float out)  │ ║
║ └────────────────────────────────────────┬─────────────────────────────────────────┘ ║
║                                          ▼ planar f64 at the source rate             ║
║ ┌──────────────────────────────────────────────────────────────────────────────────┐ ║
║ │ pass 0 scan    fu_cliplevel: pinned levels per channel | fu_cutoff: wall in Hz   │ ║
║ └────────────────────────────────────────┬─────────────────────────────────────────┘ ║
║                                          ▼ ceilings, src_top, ref_top                ║
║ ┌──────────────────────────────────────────────────────────────────────────────────┐ ║
║ │ fu_engine (passes 1 and 2)                                                       │ ║
║ │ 1 fu_declip    Opus LPC-24, two-sided least squares    +8.73 dB SDR at 0.3 %     │ ║
║ │ 2 fu_resample  2/1 x 256 taps, then 320/147 x 32 taps  236.8x alone, -170.6 dB   │ ║
║ │ 3 fu_extend    fold wall .. 25 kHz; harmonics above    guard -30 dB per frame    │ ║
║ └────────────────────────────────────────┬─────────────────────────────────────────┘ ║
║                                          ▼ planar f64 at 192000 Hz                   ║
║ ┌──────────────────────────────────────────────────────────────────────────────────┐ ║
║ │ pass 1: true peak -> gain    pass 2: fu_pcm f64 -> i32 (14.0x plain) -> libFLAC  │ ║
║ └────────────────────────────────────────┬─────────────────────────────────────────┘ ║
║                                          ▼ 192000 Hz, 32 bit, tags + cover carried   ║
║ ┌──────────────────────────────────────────────────────────────────────────────────┐ ║
║ │ out/<dir>/<name>.192k32.flac     (+ <name>.flac, 24 bit, for a lossy input)      │ ║
║ └──────────────────────────────────────────────────────────────────────────────────┘ ║
╚══════════════════════════════════════════════════════════════════════════════════════╝
```

## Ground rules

- **Double precision from decoder to encoder.** Samples are planar `double` at
  +-1.0 through every stage. The only rounding is the final conversion to
  32-bit integers.
- **Nothing is clipped.** The file is run twice: once to measure the true peak
  at the output rate, once to write it with the smallest uniform attenuation
  that fits. A file whose peak already fits is written at unity gain.
- **The base band is not processed by the extension.** Synthesis is added above
  the source's bandwidth; the band below it passes through the resampler only.
- **One streaming object.** `fu_engine` holds the whole path and is driven
  block by block. The command line tool drives it for files; a real-time
  player can drive the same object.

## What comes from Opus

Opus protects against lost packets in two ways (forward error correction and
concealment) and uses band folding to save bits. FLACscale reworks the
concealment's linear predictor and the band folding.

| Opus mechanism | What it is | In FLACscale |
|---|---|---|
| In-band FEC (LBRR) | The encoder puts a low-resolution copy of each frame into the next packet. | Not used. It needs redundancy made by an encoder; a finished file contains none. |
| Packet loss concealment, LPC part (`celt_decode_lost`) | The decoder fits a 24th-order linear predictor to the last 1024 samples and continues the signal from it. | Ported (`fu_lpc`) and reworked (`fu_declip`): a clipped run has intact samples on both sides, so the lost samples are solved from both. |
| Packet loss concealment, pitch part | The decoder repeats the last pitch period of the excitation. | Not ported. Its premise (a tonal sound is a comb of partials) is what the harmonic stage of `fu_extend` is built on. |
| Band folding (`celt/bands.c`) | A band that received no bits is filled with a copy of lower-band spectrum, rescaled to the band's energy. | Reworked (`fu_extend`): the missing band above the source's bandwidth is filled the same way, with the energy extrapolated instead of transmitted. |

## Pass 0: scan

The source is decoded once before any processing.

**Clip levels (`fu_cliplevel`).** Per channel, a 65537-bin histogram on a
16-bit grid plus the exact extremes. Hard clipping piles samples onto the
extreme level. A side counts as clipped when its extreme bin holds at least 32
samples and at least 8 times the mean of the 7 bins just inside it. The ceiling
handed to the repair is the inner edge of that bin, because a 24-bit limiter
ceiling wobbles inside one 16-bit step. Lossy sources skip this: the encoder
smeared the flat tops, so they cannot be located.

**Bandwidth (`fu_cutoff`).** Averaged power spectrum (FFT 2048, Hann) in 250 Hz
bands. The "wall" is the lowest boundary at or above 8 kHz where every band
above it, up to Nyquist, lies at least 20 dB under the median of the 2 kHz
below it. No wall means the content reaches Nyquist.

## Stage 1: clipped-peak repair (`fu_declip`, source rate)

Model: `e[n] = x[n] + sum(lpc[k-1] * x[n-k])`, order 24, analysed on up to
1024 samples with the CELT window taper, exactly as `celt_decode_lost` does.

- Pinned samples no further apart than the LPC order form a **cluster** (at
  most 256 unknowns).
- The unknowns of a cluster are one least-squares problem: minimise the
  excitation energy over every row the unknowns touch. It is solved by
  Cholesky factorisation. With no known samples on the right, the solution is
  the Opus synthesis recursion itself.
- **Constraint Opus never has:** the true value of a pinned sample lies at or
  beyond its pinned value, with known sign. A solution that falls short is held
  at its pinned value and the others are re-solved (active set).
- **Explosion guard:** a cluster whose solution exceeds 4 times the largest
  pinned magnitude in it is refused and left as it was.
- One measured change to Opus's numeric guards: the noise floor on lag 0 is
  -80 dB instead of -40 dB (+0.34 to +0.83 dB SDR). A second LPC pass was
  measured and dropped: it changes the SDR by -0.10 to +0.12 dB and overshoots
  (highest repaired peak +1.76 dBFS against +0.07 dBFS with one pass).

## Stage 2: resampling (`fu_resample`)

Rational polyphase, Kaiser-windowed sinc, zero delay: output `n` sits at input
time `n * M / L`, and `n` input frames give exactly `ceil(n * L / M)` output
frames.

- **Two stages for 44100 -> 192000 (640/147):** a steep 2/1 stage (256 taps per
  output, pass 19999 Hz, stop 22050 Hz), then a relaxed 320/147 stage (32
  taps), which only has to remove images that start 44 kHz away.
- **Visit-order coefficient table:** the polyphase rows are stored in the order
  the outputs visit them, with a byte table for how far the input window
  advances after each row. The kernel walks memory linearly.
- **Stereo-fused FMA kernel:** both channels share one pass over the
  coefficients, 4 accumulator chains per channel (`fu_dot.h`).

Rates below 2x the source rate use one steep stage.

## Stage 3: band extension (`fu_extend`, output rate)

Short-time Fourier transform at the output rate: frames of 4096, hop 1024,
Hann window on analysis and synthesis. Bins are phase-referenced to the frame
centre, so a copy shifted by a multiple of 4 bins overlap-adds coherently.

**From the wall to 25 kHz: folding.** Spectrum from just below the wall is
copied upward, normalised by its own local level and rescaled to an
extrapolated envelope: the level of the top real kilohertz, continued along the
slope fitted to the top four real kilohertz (clamped to -12 .. -1 dB/kHz).

**Above 25 kHz: integer harmonics only.** A bin that stands 15 dB above its
local floor is a tonal peak. Each peak gets its multiples h = 2 .. 6 that land
above 25 kHz, with phase `h * phi` and amplitude reduced by the source's own
spectral tilt (clamped to 6 .. 24 dB per octave) times `log2(h)`. Noise-like
content produces no peak and therefore nothing above 25 kHz.

**Guard.** Per frame, the synthesised power above 25 kHz is capped at -30 dB of
the real signal's power in that frame. A wrong guess cannot load a tweeter.

## Level and output (`fu_pcm`, `fu_flac`)

`fu_pcm_to_int` scales, rounds (to nearest, ties to even), saturates and counts
clipped samples in one branchless AVX2 pass; the count comes from a population
count of the comparison masks. `fu_flac` writes 32-bit FLAC (compression level
5, a seek point every 10 s, threaded encoder), carries the source's tags and
cover art, and adds a `FLACSCALE=` comment that records the conversion.

A lossy input is also written as a plain 24-bit FLAC at the source rate. The
upscaler is fed from the decoder's floating-point output, not from that file.

## Modules

| Module | Role |
|---|---|
| `flacscale.c` | the command: options, the three passes, the report |
| `fu_engine` | the whole signal path as one streaming object |
| `fu_input`, `fu_flac`, `fu_mp3`, `fu_snd`, `fu_tags` | readers (format recognised by content), FLAC writer, tag carrying |
| `fu_cliplevel`, `fu_cutoff` | pass 0 detectors |
| `fu_lpc`, `fu_declip` | Opus LPC analysis and the two-sided repair |
| `fu_resample`, `fu_dot.h` | resampler chain and its FMA kernels |
| `fu_resample_ref` | plain reference forms, used only by `rsbench` |
| `fu_extend` | folding and harmonic synthesis |
| `fu_pcm` | double to integer conversion, peak measurement |
| `flacspec`, `flacdiff`, `clipscan` | measurement tools |
| `rsbench`, `pcmcheck`, `declipcheck`, `tonegen` | gates |

## What is hand-vectorised, and what is not

Hand-written AVX2/FMA: the resampler kernels (`fu_dot.h`) and the PCM
conversion and peak search (`fu_pcm.c`). Each is checked against a plain form
for equal results and must not be slower than it (`rsbench`, `pcmcheck`).

Plain C, vectorised only as far as the compiler manages: the LPC
autocorrelation, the least-squares solver and the extension stage. They are the
next targets; see [ROADMAP.md](ROADMAP.md).
