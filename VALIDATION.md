# Validation against MATLAB

This file records how `wsstpy` was checked against MATLAB, and the exact
numbers the check produced. The numbers are not rounded or summarised: they
are the output of the scripts in `probe/`, printed by `probe/compare5.py` and
`probe/compare6.py`.

## Environment

| | |
| --- | --- |
| MATLAB | 25.2.0.3177638 (R2025b) Update 5, Wavelet Toolbox |
| Operating system | Debian GNU/Linux 13, x86_64, kernel 6.12 |
| CPU | Intel Core i7-4870HQ |
| Python | 3.13.5 |
| NumPy | 2.5.4 |
| SciPy | 1.18.1 |

## Method

`probe/probe5.m` and `probe/probe6.m` run the MATLAB functions and save every
output in double precision. `probe/compare5.py` and `probe/compare6.py` run
the Python functions on the same inputs and print the element-wise difference.
No post-processing, averaging or rounding is applied before the comparison.

Test signals:

| name | definition |
| --- | --- |
| tone | `cos(2*pi*0.05*k) + 0.5*cos(2*pi*0.20*k)` |
| chirp | `cos(2*pi*(0.02*k + 0.4*k.^2/N))` |
| noise | `randn(N,1)` with `rng(0)` |
| short | `cos(2*pi*0.13*k) + 0.3*sin(2*pi*0.31*k)` |
| impulse | one sample equal to 1 at `floor(N/3)` |

Lengths 256, 512, 1000, 1024, 1025 and 4096. `k = 0..N-1`.

## Forward transform

`wsst` with the analytic Morlet wavelet, 32 voices per octave, periodic
boundary (the MATLAB default), and with 10 voices and with `ExtendSignal`.

| comparison | max abs diff | max rel diff | RMS diff |
| --- | --- | --- | --- |
| wsst tone N=256 | 2.168e-15 | 1.008e-14 | 7.713e-17 |
| wsst vpo=10 N=256 | 2.959e-16 | 1.372e-15 | 2.455e-17 |
| wsst extend N=256 | 2.445e-15 | 1.137e-14 | 7.268e-17 |
| wsst tone N=512 | 2.478e-15 | 1.152e-14 | 6.406e-17 |
| wsst vpo=10 N=512 | 4.282e-16 | 1.981e-15 | 2.507e-17 |
| wsst extend N=512 | 2.162e-15 | 1.005e-14 | 4.764e-17 |
| wsst tone N=1000 | 7.005e-16 | 3.256e-15 | 3.716e-17 |
| wsst vpo=10 N=1000 | 2.243e-16 | 1.042e-15 | 1.438e-17 |
| wsst extend N=1000 | 1.670e-15 | 7.761e-15 | 3.907e-17 |
| wsst tone N=1024 | 2.205e-15 | 1.025e-14 | 4.598e-17 |
| wsst vpo=10 N=1024 | 4.475e-16 | 2.064e-15 | 1.663e-17 |
| wsst extend N=1024 | 2.071e-15 | 9.628e-15 | 4.159e-17 |
| wsst tone N=1025 | 2.134e-15 | 9.919e-15 | 4.911e-17 |
| wsst vpo=10 N=1025 | 4.679e-16 | 2.160e-15 | 2.481e-17 |
| wsst extend N=1025 | 2.103e-15 | 9.776e-15 | 4.600e-17 |
| wsst tone N=4096 | 2.168e-15 | 1.008e-14 | 3.399e-17 |
| wsst vpo=10 N=4096 | 3.310e-16 | 1.535e-15 | 1.402e-17 |
| wsst extend N=4096 | 2.448e-15 | 1.138e-14 | 3.464e-17 |
| wsst morse (3,60) | 2.178e-15 | 1.330e-14 | 5.159e-17 |
| wsst morse (3,30) | 2.185e-15 | 9.321e-15 | 5.098e-17 |
| wsst morse (2,40) | 2.567e-15 | 1.285e-14 | 8.927e-17 |
| wsst bump | 1.939e-15 | 2.662e-14 | 7.721e-17 |
| wsst fs=100 | 2.205e-15 | 1.025e-14 | 4.598e-17 |

Max relative difference is below 3e-14 in every case. For comparison, double
precision carries about 1e-16 relative resolution, so these are differences of
a few units in the last place.

## Waveforms other than steady tones

These are the cases where the phase transform moves between rows, so they
exercise the reassignment rule rather than a single fixed row.

| case | max abs diff | max rel diff |
| --- | --- | --- |
| wsst chirp N=256 | 8.951e-16 | 4.403e-15 |
| wsst noise N=256 | 8.365e-16 | 5.470e-15 |
| wsst short N=256 | 1.330e-15 | 6.272e-15 |
| wsst impulse N=256 | 1.982e-17 | 2.273e-15 |
| wsst chirp N=512 | 1.183e-15 | 5.672e-15 |
| wsst noise N=512 | 1.444e-15 | 9.355e-15 |
| wsst short N=512 | 1.430e-15 | 6.742e-15 |
| wsst impulse N=512 | 1.528e-17 | 1.499e-15 |
| wsst chirp N=1000 | 1.585e-15 | 7.424e-15 |
| wsst noise N=1000 | 1.065e-15 | 6.992e-15 |
| wsst short N=1000 | 7.920e-16 | 3.734e-15 |
| wsst impulse N=1000 | 1.508e-17 | 1.485e-15 |
| wsst chirp N=1024 | 1.590e-15 | 7.431e-15 |
| wsst noise N=1024 | 1.579e-15 | 8.438e-15 |
| wsst short N=1024 | 1.105e-15 | 5.211e-15 |
| wsst impulse N=1024 | 1.541e-17 | 1.517e-15 |
| wsst chirp N=1025 | 1.660e-15 | 7.765e-15 |
| wsst noise N=1025 | 1.269e-15 | 5.541e-15 |
| wsst short N=1025 | 1.463e-15 | 6.899e-15 |
| wsst impulse N=1025 | 1.778e-17 | 2.026e-15 |
| wsst chirp N=4096 | 2.365e-15 | 1.099e-14 |
| wsst noise N=4096 | 1.356e-15 | 5.943e-15 |
| wsst short N=4096 | 1.438e-15 | 6.778e-15 |
| wsst impulse N=4096 | 1.538e-17 | 1.752e-15 |

The impulse case is smaller because its synchrosqueezed matrix is sparse: few
coefficients survive the threshold, so there are fewer values to differ in.

## Continuous wavelet transform

| comparison | max abs diff | max rel diff | RMS diff |
| --- | --- | --- | --- |
| cwt Morse default | 8.111e-15 | 8.136e-15 | 1.544e-15 |
| cwt Morse limits [0.02 0.25] | 5.997e-15 | 6.035e-15 | 1.341e-15 |
| cwt amor 32 voices | 2.143e-14 | 2.085e-14 | 4.735e-15 |

## Inverse transforms

`iwsst` and `icwt` differ from MATLAB by about 1e-8 relative in two cases.
That is not a difference in the transform. It comes from the admissibility
constant `C_psi`, a single number that both implementations evaluate by
numerical integration. The integrand behaves like `1/w` near zero, so the
value is set by where the integration stops, and MATLAB stops at `sqrt(eps)`
with `integral`'s default relative tolerance of 1e-6.

The constants implied by the MATLAB outputs can be recovered from the
reconstructions:

| constant | MATLAB output implies | wsstpy | relative difference |
| --- | --- | --- | --- |
| `R_psi`, amor | 0.4305137819610273 | 0.43051377150288545 | 2.43e-08 |
| `C_psi`, amor | 0.861027604938177 | 0.8610275430057707 | 7.19e-08 |

The wsstpy values were checked two ways (adaptive quadrature and an analytic
split of the singular part) and agree with each other to 5.6e-17. MATLAB's
values carry the error of its own quadrature. The Morse constant has a closed
form, and there the two implementations agree exactly:

| comparison | max abs diff | max rel diff | RMS diff |
| --- | --- | --- | --- |
| icwt full | 6.217e-15 | 4.129e-15 | 2.823e-15 |
| icwt band [0.03 0.08] | 4.663e-15 | 4.313e-15 | 2.599e-15 |
| icwt vpo=10 | 6.217e-15 | 4.129e-15 | 2.823e-15 |
| icwt Morse limits [0.03 0.1] | 4.663e-15 | 4.279e-15 | 2.525e-15 |
| icwt amor | 3.461e-07 | 7.193e-08 | 1.819e-07 |
| iwsst full | 3.651e-08 | 2.429e-08 | 1.922e-08 |
| iwsst band [0.03 0.08] | 2.583e-08 | 2.429e-08 | 1.716e-08 |
| iwsst by ridge | 2.545e-08 | 2.429e-08 | 1.356e-08 |
| iwsst Morse | 4.441e-15 | 2.952e-15 | 2.083e-15 |

## Ridge extraction

`wsstridge` returns row indices, so the comparison is exact equality, not a
tolerance. The reference is MATLAB's own output on its own synchrosqueezed
matrix.

| case | identical rows |
| --- | --- |
| penalty 0 | 1024/1024 |
| penalty 0.5 | 1024/1024 |
| penalty 1 | 1024/1024 |
| penalty 5 | 1024/1024 |
| penalty 10 | 1024/1024 |
| penalty 100 | 1024/1024 |
| 2 ridges | 2048/2048 |
| 3 ridges, 6 frequency bins | 3072/3072 |
| chirp, penalty 2 | 1024/1024 |
| noise, penalty 3 | 1024/1024 |

The ridge penalty is applied to the squared row difference. With a linear
penalty 1015 of 1024 rows matched on the chirp and 769 of 1024 on the noise,
which is how the squared form was identified.

## Reproducing this

On a machine with MATLAB and the Wavelet Toolbox:

```bash
cd probe
matlab -batch probe5     # writes data5/ref5.mat
matlab -batch probe6     # writes data6/ref6.mat
matlab -batch probe7     # writes data7/ref7.mat (small fixtures)
cd ..
python probe/compare5.py     # prints the first table
python probe/compare6.py     # prints the waveform table
python scripts/make_fixtures.py
```

Without MATLAB, `pytest` runs the same comparisons against the committed
fixtures in `tests/data/reference.npz`, which are the arrays from
`probe7.m`.
