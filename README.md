# wsstpy

Wavelet synchrosqueezing in Python, with the same numbers as MATLAB's
Wavelet Toolbox.

| MATLAB | wsstpy |
| --- | --- |
| `wsst(x)` | `wsstpy.wsst(x)` |
| `iwsst(sst)` | `wsstpy.iwsst(sst)` |
| `wsstridge(sst)` | `wsstpy.wsstridge(sst)` |
| `cwt(x)` | `wsstpy.cwt(x)` |
| `icwt(cfs)` | `wsstpy.icwt(cfs)` |
| `cwtfilterbank` | `wsstpy.CWTFilterBank` |

## Accuracy

Checked against MATLAB R2025b, signal lengths 256 to 4096, Morlet, Morse and
bump:

| quantity | maximum difference |
| --- | --- |
| `wsst` matrix | 2.6e-15 |
| `cwt` coefficients | 2.2e-14 |
| `iwsst` / `icwt`, Morse | 6.3e-15 |
| `iwsst` / `icwt`, Morlet | 7.2e-08 relative |
| `wsstridge` row indices | equal in every case tested |

Ridge indices are integers and match exactly. Continuous arrays agree to a few
units in the last place, not bit for bit, because the two sides sum and FFT in
a different order.

The Morlet inverse differs at 1e-8 relative because MATLAB evaluates the
admissibility constant with `integral` at its default tolerance of 1e-6.
wsstpy computes 0.8610275430057707 where MATLAB's output implies
0.861027604938177. Tables for every case are in [VALIDATION.md](VALIDATION.md).

![MATLAB and wsstpy output, and the difference](figures/sst_comparison.png)

Left and middle are the same matrix from MATLAB and from wsstpy. Right is the
difference on a log scale, peaking at 1.6e-15.

## Install

```bash
pip install git+https://github.com/Halhousni/wsstpy.git
```

Python 3.9 or newer, NumPy, SciPy. Not on PyPI yet.

## Usage

```python
import numpy as np
import wsstpy

fs = 1000.0
k = np.arange(2000)
x = np.cos(2 * np.pi * 50 * k / fs) + 0.7 * np.cos(2 * np.pi * 180 * k / fs)

sst, f = wsstpy.wsst(x, fs)                    # complex (n_scales, n), f in Hz
fridge, iridge = wsstpy.wsstridge(sst, penalty=0.5, num_ridges=2, f=f)
tones = wsstpy.iwsst(sst, ridge_indices=iridge)  # one signal per ridge
xrec = wsstpy.iwsst(sst)                         # whole signal back
```

Plot it with `plt.pcolormesh(np.arange(x.size), f, np.abs(sst))`. Runnable
examples are in `examples/`.

## API

| function | does |
| --- | --- |
| `wsst(x, fs=None, wavelet="amor", voices_per_octave=32, extend_signal=False)` | forward transform, rows low to high frequency |
| `iwsst(sst, f=None, freq_range=None, ridge_indices=None)` | inverse, full, band or per ridge |
| `wsstridge(sst, penalty=0, num_ridges=1, f=None)` | ridge rows, penalty on the squared row change |
| `cwt(x, fs=None, wavelet="morse", voices_per_octave=10)` | continuous wavelet transform |
| `icwt(cfs, f=None, freq_range=None)` | inverse CWT |
| `CWTFilterBank(signal_length, ...)` | filters, scales, centre frequencies |

`wavelet` is `"amor"`, `"morse"` or `"bump"`; Morse takes
`wavelet_parameters=(gamma, time_bandwidth)`.

## How it works

Filters are the wavelet's Fourier transform on a geometric scale grid, peaking
at magnitude 2. The CWT comes from multiplying the signal's FFT by each filter,
and its derivative from the same product with `1j * omega`. The phase transform
`imag(dCWT / CWT) / (2 pi)` gives each coefficient an instantaneous frequency,
and the coefficient moves to the row nearest that frequency on a log scale,
weighted by `log(2) / voices_per_octave`. Coefficients below `sqrt(eps)` are
dropped. Ridges minimise `-log|SST|` plus `penalty * (i - j)^2` along a path.

## Notes

- `voices_per_octave` accepts any even integer from 2 up; MATLAB allows 10 to
  48. Same value, same result.
- `ScalingCoefficients` and the low-pass term of `cwt` are not implemented.
  Use `extend_signal=True`, which is MATLAB's default.
- MATLAB's `wsst` plots when called without outputs. Here it always returns
  arrays.
- `icwt` does not accept time tables.

## Licence and provenance

MIT, see [LICENSE](LICENSE). No MATLAB source code is used; the audit and the
list of files read for conventions are in [PROVENANCE.md](PROVENANCE.md).

## References

Daubechies, Lu and Wu, *Appl. Comput. Harmon. Anal.* 30(2), 2011. Thakur,
Brevdo, Fuckar and Wu, *Signal Processing* 93(5), 2013.
