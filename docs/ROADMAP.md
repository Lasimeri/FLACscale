# Roadmap

What exists, what does not, and what is known to be weak. Nothing below the
first section is built.

## Done

- Resampler: two-stage polyphase chain, gated against a scalar formula and
  against libsoxr.
- File path: FLAC, WAV, AIFF, MP3 in; 32-bit FLAC out; true-peak gain; tags and
  cover art carried over.
- Clipped-peak repair from the Opus LPC model, gated against a known truth.
- Band extension: folding to 25 kHz, harmonics only above, with a power guard.
- `fu_engine`: the whole path as one streaming object.

## Planned

1. **Real-time player.** Drives `fu_engine` block by block into the sound
   device at the output rate. It needs a device path that really runs at
   192 kHz (a sound server locked to 48 kHz would resample the result back
   down), and a hard volume ceiling of its own.
2. **Graphical front end**, in Rust, on top of the player. It comes last.

## Known weaknesses, in the order they hurt

1. **The pitch half of Opus's concealment is not ported.** The repair uses the
   short-term predictor only. Long clipped runs in strongly periodic material
   are where a pitch predictor would help most.
2. **Speed against libsoxr.** Its VHQ recipe is 1.40x faster than the
   resampler here.
3. **Scalar stages.** The LPC autocorrelation, the least-squares solver and the
   extension stage are plain C. The extension stage plans its transforms with
   `FFTW_ESTIMATE`.
4. **One outlier sets the level.** A single repaired peak far above the rest
   decides the attenuation of the whole file.
5. **Detector thresholds are fixed constants.** 32 samples at the extreme for
   "clipped", a 20 dB drop for a bandwidth wall. Lightly clipped masters and
   gentle roll-offs fall under them.
6. **No dequantisation.** The 16-bit rounding of the source is carried into the
   output unchanged.
7. **Tested on one machine and one distribution.**
