# FLACscale -- build. Targets land in out/. Flags match src/compile_flags.txt.
#
#   make          build every tool into out/
#   make check    run the gates that need no audio file of yours (about 30 s)
CC     := gcc
CFLAGS := -std=gnu2x -O3 -mavx2 -mfma -mbmi2 -mtune=znver3 -Wall -Wextra -Werror -Isrc

HEADERS := $(wildcard src/*.h)
INPUT   := src/fu_input.c src/fu_flac.c src/fu_mp3.c src/fu_snd.c src/fu_tags.c
INLIBS  := -lFLAC -lmpg123 -lsndfile
DECLIP  := src/fu_declip.c src/fu_lpc.c
ENGINE  := src/fu_engine.c src/fu_resample.c src/fu_cliplevel.c src/fu_cutoff.c src/fu_extend.c $(DECLIP)

TOOLS := flacscale flacspec flacdiff clipscan declipcheck dequantcheck rsbench pcmcheck tonegen
BINS  := $(TOOLS:%=out/%)
all: $(BINS)

$(BINS): | out
out:
	mkdir -p out

out/flacscale: src/flacscale.c src/fu_pcm.c $(ENGINE) $(INPUT) $(HEADERS)
	$(CC) $(CFLAGS) -o $@ src/flacscale.c src/fu_pcm.c $(ENGINE) $(INPUT) $(INLIBS) -lfftw3 -lm

out/flacspec: src/flacspec.c $(INPUT) $(HEADERS)
	$(CC) $(CFLAGS) -o $@ src/flacspec.c $(INPUT) $(INLIBS) -lfftw3 -lm

out/flacdiff: src/flacdiff.c $(INPUT) $(HEADERS)
	$(CC) $(CFLAGS) -o $@ src/flacdiff.c $(INPUT) $(INLIBS) -lm

out/clipscan: src/clipscan.c $(INPUT) $(HEADERS)
	$(CC) $(CFLAGS) -o $@ src/clipscan.c $(INPUT) $(INLIBS) -lm

out/declipcheck: src/declipcheck.c $(DECLIP) $(INPUT) $(HEADERS)
	$(CC) $(CFLAGS) -o $@ src/declipcheck.c $(DECLIP) $(INPUT) $(INLIBS) -lm

out/dequantcheck: src/dequantcheck.c src/fu_dequant.c $(INPUT) $(HEADERS)
	$(CC) $(CFLAGS) -fopenmp -o $@ src/dequantcheck.c src/fu_dequant.c $(INPUT) $(INLIBS) -lfftw3 -lm

out/rsbench: src/rsbench.c src/fu_resample_ref.c src/fu_resample.c $(HEADERS)
	$(CC) $(CFLAGS) -o $@ src/rsbench.c src/fu_resample_ref.c src/fu_resample.c -lsoxr -lm

out/pcmcheck: src/pcmcheck.c src/fu_pcm.c $(HEADERS)
	$(CC) $(CFLAGS) -o $@ src/pcmcheck.c src/fu_pcm.c -lm

out/tonegen: src/tonegen.c
	$(CC) $(CFLAGS) -o $@ src/tonegen.c -lm

# Gates on synthetic signals only. pcmcheck and rsbench exit non-zero on a
# failure; the last three lines print the band-extension result for the
# synthetic tones, to be read against docs/MEASUREMENTS.md.
check: all
	out/pcmcheck
	out/rsbench
	out/tonegen out/check-tones.wav
	out/flacscale out/check-tones.wav -o out/check-tones.192k32.flac
	out/flacspec out/check-tones.192k32.flac -l 25000 96000

.PHONY: all check
