# Third-party notices

FLACscale's own code is under the MIT licence (see `LICENSE`). This file lists
what comes from elsewhere.

## Code ported from libopus

`src/fu_lpc.c` and `src/fu_lpc.h` contain routines ported from libopus 1.6.1
and rewritten in double precision:

| FLACscale | libopus 1.6.1 | Copyright in the original file |
|---|---|---|
| `celt_window` | `celt/modes.c:392` | (c) 2007-2008 CSIRO, (c) 2007-2009 Xiph.Org Foundation, (c) 2008 Gregory Maxwell; written by Jean-Marc Valin and Gregory Maxwell |
| `fu_autocorr` | `celt/celt_lpc.c:281` (`_celt_autocorr`) | (c) 2009-2010 Xiph.Org Foundation; written by Jean-Marc Valin |
| `fu_levinson` | `celt/celt_lpc.c:37` (`_celt_lpc`) | (c) 2009-2010 Xiph.Org Foundation; written by Jean-Marc Valin |
| `fu_lpc_analyze` | the analysis sequence of `celt_decode_lost`, `celt/celt_decoder.c:856-877` | (c) 2007-2008 CSIRO, (c) 2007-2010 Xiph.Org Foundation, (c) 2008 Gregory Maxwell; written by Jean-Marc Valin and Gregory Maxwell |

Two other modules take a method from libopus and no code: `src/fu_declip.c`
(the excitation model of `celt_decode_lost`) and `src/fu_extend.c` (band
folding as in `celt/bands.c`).

The ported routines remain under the libopus licence, reproduced here in full
from the `COPYING` file of libopus 1.6.1:

```
Copyright 2001-2023 Xiph.Org, Skype Limited, Octasic,
                    Jean-Marc Valin, Timothy B. Terriberry,
                    CSIRO, Gregory Maxwell, Mark Borgerding,
                    Erik de Castro Lopo, Mozilla, Amazon

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions
are met:

- Redistributions of source code must retain the above copyright
notice, this list of conditions and the following disclaimer.

- Redistributions in binary form must reproduce the above copyright
notice, this list of conditions and the following disclaimer in the
documentation and/or other materials provided with the distribution.

- Neither the name of Internet Society, IETF or IETF Trust, nor the
names of specific contributors, may be used to endorse or promote
products derived from this software without specific prior written
permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

Opus is subject to the royalty-free patent licenses which are
specified at:

Xiph.Org Foundation:
https://datatracker.ietf.org/ipr/1524/

Microsoft Corporation:
https://datatracker.ietf.org/ipr/1914/

Broadcom Corporation:
https://datatracker.ietf.org/ipr/1526/
```

FLACscale is not Opus, is not endorsed by its authors, and neither reads nor
writes Opus streams.

## Libraries linked at build time

No library source is included in this repository. The programs link against
the system's shared libraries:

| Library | Used by | Licence |
|---|---|---|
| libFLAC | every tool that reads or writes FLAC | BSD-3-Clause (Xiph.Org) |
| libmpg123 | MP3 input | LGPL-2.1 |
| libsndfile | WAV, RF64, AIFF, CAF input | LGPL-2.1-or-later |
| FFTW 3 | `flacscale` (band extension, bandwidth scan), `flacspec` | GPL-2.0-or-later |
| libsoxr | `rsbench` only, as the speed and accuracy baseline | LGPL-2.1-or-later |

Consequence for anyone who redistributes binaries: `flacscale` and `flacspec`
link FFTW, so a binary of either is a combined work that must be distributed
under the terms of the GPL. The source in this repository stays MIT.
