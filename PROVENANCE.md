# Provenance

This repository contains no MATLAB source code. Every file was written for
this project and is released under the MIT licence in `LICENSE`.

The Python package is an implementation of the published synchrosqueezing
algorithm. MATLAB R2025b with the Wavelet Toolbox was used as a reference
oracle: its functions were run on test signals, and the outputs were compared
with the Python implementation. MATLAB's `.m` files were read to establish
conventions that a port has to reproduce if the numbers are to agree, such as
where the frequency grid starts, how the filters are normalised, and which
threshold discards small coefficients.

## Files consulted for conventions

| MATLAB file | Convention taken from it |
| --- | --- |
| `wavelet/wsst.m` | order of operations, `sqrt(eps)` coefficient threshold, log-frequency row mapping, signal mean removal, reflection padding rule |
| `wavelet/iwsst.m` | reconstruction as `2/R_psi * real(sum(sst))`, ridge-restricted reconstruction, `NumFrequencyBins` default of 16 |
| `wavelet/wsstridge.m` | argument order, ridge frequency vector, defaults for penalty and bins |
| `wavelet/cwt.m` | defaults for the CWT entry point (Morse, 10 voices, reflection, `[3 60]`) |
| `wavelet/icwt.m` | scale ratio from the frequency vector, admissibility constant, band selection |
| `@cwtfilterbank/cwtfilterbank.m` | filter layout, frequency grid, centre frequency definition |
| `@cwtfilterbank/wt.m` | analysis by multiplication in the frequency domain, padding |
| `+cwt/wavbpfilters.m` | Morlet and bump responses, peak magnitude 2 |
| `+cwt/morsebpfilters.m` | Morse response and its factor |
| `+cwt/getCWTScales.m`, `+cwt/cwtfreqlimits.m` | scale grid and the coarsest scale implied by the signal length |
| `+cwt/wavCFandSD.m` | centre frequency and time-domain width per wavelet |
| `+cwt/findFreqIndices.m` | which rows a frequency band keeps |
| `+cwt/morsenormconstant.m`, `+cwt/morsepeakfreq.m`, `+cwt/morseproperties.m` | Morse normalisation and peak frequency |
| `+cwt/admConstant.m` | admissibility constants used by the inverse transforms |
| `+signalwavelet/+internal/+tfridge/extractRidges.m` | cost `-log|sst|`, squaring of the penalty, blanking of rows between ridges |

Nothing in that list was transliterated. The Python is structured differently:
vectorised over scales and columns, with the filter bank as an object, explicit
band handling, and ridge extraction by a dynamic programme written from the
recurrence rather than from the compiled helper.

## Automated check

`scripts/check_provenance.py` compares this repository against an installed
MATLAB tree. It runs three passes:

1. Marker scan for strings that only occur in toolbox source: `coder.`
   pragmas, `validateattributes`, `bsxfun`, `accumarray`, `quadgk`,
   `underlyingType`, `@cwtfilterbank`, MathWorks copyright lines.
2. Eight-token window matching against every `.m` file in the Wavelet
   Toolbox, the shared signalwavelet package and the signal toolbox, keeping
   only windows that contain at least four identifier tokens.
3. Longest common token run per repository file against the eighteen files
   listed above, reporting only runs that contain five or more identifiers.

The script and this file are excluded from the scan, because both quote the
marker strings and the example snippets on purpose.

Result at the commit that added this file, MATLAB R2025b on Debian 13:

```
repository text files: 34
pass 1, marker scan: 1 hit(s)
  README.md:201  [MathWorks]  This is not affiliated with MathWorks and ...
pass 2, 6037 .m files, 1746663 distinct windows
  shared 8-token windows: 3
    probe/probe6.m:  ; end for k = 1 : numel
    probe/probe6.m:  end for k = 1 : numel (
    probe/probe6.m:  sst , f ] = wsst ( x
pass 3, longest shared run: 0 tokens
pass 4, 958 .m files, 382343 distinct word windows
  files with a run of 8 or more shared words: 0
```

The single marker hit is the word "MathWorks" in the README's disclaimer. The
three window hits are in `probe/probe6.m`, a script written here that calls
MATLAB's public API: a `for` loop header and the call `[sst, f] = wsst(x)`.
Neither is toolbox source. No file shares a run of five or more consecutive
identifier tokens with any consulted MATLAB file.

Reproduce with:

```bash
python scripts/check_provenance.py --matlab /usr/local/MATLAB/R2025b
```

The check needs a MATLAB installation. It is not part of CI for that reason.

## What the shared formulas are

Some expressions necessarily look alike because they are the same mathematics:
the analytic Morlet response `2 exp(-(s omega - 6)^2 / 2)`, the Morse response
`2 a w^beta exp(-w^gamma)`, the scale grid `s0 2^(k/nv)`, and the reassignment
weight `log 2 / nv`. These are definitions from the literature (Daubechies,
Lu and Wu 2011; Thakur et al. 2013; Olhede and Walden 2002 for the Morse
wavelet), and the normalisation to peak magnitude 2 is the L1 convention that
makes MATLAB's coefficients come out at the scale they do. A port that changed
them would return different numbers.

## Test fixtures

`tests/data/reference.npz` holds numeric output: synchrosqueezed matrices,
coefficient arrays and ridge indices produced by running MATLAB on the test
signals in `probe/probe7.m`. It contains no MATLAB code. The fixtures exist so
that `pytest` runs without a MATLAB licence. They can be regenerated from
scratch with `probe/probe7.m` and `scripts/make_fixtures.py`.

## Trademarks

MATLAB and Wavelet Toolbox are trademarks of The MathWorks, Inc. This project
is not affiliated with, endorsed by, or sponsored by MathWorks.
