# Contributing

Patches are welcome. These are the rules the existing code was held to.

## Build rules

- **Zero diagnostics.** The build uses `-Wall -Wextra -Werror`. Each changed
  file must also be clean under clangd:

      cd src && clangd --check=FILE.c --tweaks=NoSuchTweak

- **`make check` must pass** before and after the change.
- C23 as gcc implements it (`-std=gnu2x`). No new dependencies without a
  reason stated in the pull request.

## Measurement rules

- **A claim carries its number.** "Faster" and "more accurate" come with the
  measured figure, its unit, the command that produced it and the machine.
- **A hand-vectorised or otherwise unusual kernel needs a gate:** the plain
  form stays in the tree, the gate proves equal results and proves the unusual
  form is not slower. `pcmcheck` and `rsbench` are the pattern.
- **When the data contradicts the expectation, the data wins.** Say so in the
  change. `docs/MEASUREMENTS.md` already records one such case (the Kaiser
  design formula against the measured stopband).
- **Anything that touches the signal path** reports the audio MD5 of a
  converted file before and after (`metaflac --show-md5sum`), or explains why
  it changed.

## Test data

- **No audio files in the repository.** `input/` and `out/` are ignored.
- **Synthetic signals are labelled synthetic** and come from a generator in
  the tree (`tonegen`, or the ones inside `rsbench` and `pcmcheck`).

## Style

- One module, one job; the header states purpose and contract at the top.
- Comments say why, and cite the libopus file and line when a routine is
  ported from it. New ported code is listed in `NOTICE.md`.
