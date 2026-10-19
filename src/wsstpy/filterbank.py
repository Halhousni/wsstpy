"""Continuous-wavelet filter bank with the same layout as MATLAB's
``cwtfilterbank``.

The filter bank stores one frequency-domain filter per scale, laid out on the
length-``n_fft`` frequency grid.  Analysis is a plain multiplication in the
frequency domain, which is what both ``cwt`` and ``wsst`` use.
"""

from __future__ import annotations

import numpy as np

from ._wavelet import SUPPORTED_WAVELETS, frequency_response
from ._wavelet import wavelet_parameters as _wav_params

__all__ = ["CWTFilterBank", "omega_grid"]


def omega_grid(n: int) -> np.ndarray:
    """Angular-frequency grid of length ``n``, in radians/sample.

    Bin 0 is DC, then positive frequencies up to Nyquist, then the negative
    frequencies in increasing order, matching ``fft``'s output ordering.
    """
    n = int(n)
    half = n // 2
    pos = 2.0 * np.pi * np.arange(1, half + 1) / n
    n_neg = (n - 1) // 2
    neg = -pos[n_neg - 1 :: -1] if n_neg > 0 else np.empty(0)
    return np.concatenate(([0.0], pos, neg))


class CWTFilterBank:
    """Frequency-domain filter bank for the analytic CWT and SST.

    Parameters
    ----------
    signal_length : int
        Number of samples of the signals to be analysed (>= 4).
    voices_per_octave : int, optional
        Number of scales per octave.  Must be even; MATLAB restricts it to
        ``10 <= voices_per_octave <= 48``, this port allows any even value
        ``>= 2``.
    wavelet : {"amor", "morse", "bump"}, optional
        Analysis wavelet.  ``"amor"`` is MATLAB's default.
    sampling_frequency : float, optional
        Sampling frequency in Hz (or any unit).  Frequencies are reported in
        the same unit.
    frequency_limits : (float, float), optional
        Band of interest in Hz.  A lower limit of ``0`` is replaced by the
        lowest frequency the filter bank can resolve, which is what MATLAB's
        ``wsst`` does.  ``None`` uses the cutoff-based default band.
    boundary : {"periodic", "reflection"}, optional
        Boundary handling.  ``"periodic"`` means no padding, ``"reflection"``
        pads by half the signal length on each side with mirrored samples.
    wavelet_parameters : (float, float), optional
        ``(gamma, time_bandwidth)`` for the Morse wavelet, default ``(3, 60)``.

    Attributes
    ----------
    scales : ndarray, shape (n_scales,)
        Dilation factors, finest scale first (highest frequency first).
    omega : ndarray, shape (n_fft,)
        Angular-frequency grid in radians/sample.
    psi_dft : ndarray, shape (n_scales, n_fft)
        Wavelet Fourier transforms, peak magnitude 2.
    center_frequencies : ndarray, shape (n_scales,)
        Peak frequency of each filter, in Hz, ordered like ``scales``.
    n_fft : int
        Length of the FFT grid, ``signal_length + 2 * pad``.
    pad : int
        Number of mirrored samples added on each side.
    """

    def __init__(
        self,
        signal_length: int,
        voices_per_octave: int = 32,
        wavelet: str = "amor",
        sampling_frequency: float = 1.0,
        frequency_limits=None,
        boundary: str = "periodic",
        wavelet_parameters=None,
    ):
        signal_length = int(signal_length)
        if signal_length < 4:
            raise ValueError("signal_length must be >= 4")
        if int(voices_per_octave) != voices_per_octave or voices_per_octave < 2:
            raise ValueError("voices_per_octave must be an integer >= 2")
        if voices_per_octave % 2:
            raise ValueError("voices_per_octave must be even (MATLAB requires 10..48)")
        name = str(wavelet).lower()
        if name not in SUPPORTED_WAVELETS:
            raise ValueError(
                f"unknown wavelet {wavelet!r}; expected one of {SUPPORTED_WAVELETS}"
            )
        if boundary not in ("periodic", "reflection"):
            raise ValueError("boundary must be 'periodic' or 'reflection'")
        fs = float(sampling_frequency)
        if not fs > 0:
            raise ValueError("sampling_frequency must be positive")

        self.signal_length = signal_length
        self.voices_per_octave = int(voices_per_octave)
        self.wavelet = name
        self.sampling_frequency = fs
        self.boundary = boundary
        self.wavelet_parameters = (
            None if wavelet_parameters is None else tuple(wavelet_parameters)
        )
        if name == "morse" and self.wavelet_parameters is None:
            self.wavelet_parameters = (3.0, 60.0)

        cf, sigma_t, gamma, beta = _wav_params(name, self.wavelet_parameters)
        self._cf = cf
        self._sigma_t = sigma_t
        self._gamma = gamma
        self._beta = beta

        self.frequency_limits = self._resolve_frequency_limits(frequency_limits, fs)

        if boundary == "reflection":
            self.pad = signal_length // 2
        else:
            self.pad = 0
        self.n_fft = signal_length + 2 * self.pad
        self.omega = omega_grid(self.n_fft)

        self.scales = self._scales()
        self.center_frequencies = (cf / (2.0 * np.pi)) / self.scales * fs
        self.psi_dft = frequency_response(
            name, self.scales[:, None] * self.omega[None, :], self.wavelet_parameters
        )

    # -- construction helpers ------------------------------------------------
    def _cutoff(self) -> float:
        """Percentage of the peak response that defines the band edge."""
        return 0.50 if self.wavelet == "morse" else 0.10

    def _min_frequency(self, num_std: int = 2) -> float:
        """Lowest frequency the filter bank can resolve, in Hz."""
        max_scale = self.signal_length / (self._sigma_t * num_std)
        fourier_factor = 2.0 * np.pi / self._cf
        return 1.0 / (max_scale * fourier_factor) * self.sampling_frequency

    def _band_edges(self, num_std: int = 2):
        """Return ``(min_freq, max_freq)`` in Hz for the default band."""
        fs = self.sampling_frequency
        max_scale = self.signal_length / (self._sigma_t * num_std)
        min_freq = 1.0 / (max_scale * (2.0 * np.pi / self._cf)) * fs

        cutoff = self._cutoff()
        if self.wavelet == "amor":
            # 2*exp(-(w-cf)^2/2) = 2*cutoff  ->  w = cf + sqrt(-2*ln(cutoff))
            omega_c = self._cf + np.sqrt(-2.0 * np.log(cutoff))
        elif self.wavelet == "bump":
            alpha = cutoff
            eps = np.finfo(float).eps
            eps_lo, eps_hi = eps, 1.0 - eps
            f = lambda om: 1.0 / (1.0 - om**2) + np.log(alpha) - np.log(2.0) - 1.0
            lo, hi = eps_lo, eps_hi
            for _ in range(200):
                mid = 0.5 * (lo + hi)
                if f(lo) * f(mid) <= 0:
                    hi = mid
                else:
                    lo = mid
            omega_c = 0.6 * 0.5 * (lo + hi) + self._cf
        else:  # morse
            from scipy.optimize import brentq

            ga, beta = self._gamma, self._beta
            anorm = 2.0 * np.exp(beta / ga * (1.0 + np.log(ga) - np.log(beta)))
            f = lambda om: 2.0 * cutoff - anorm * om**beta * np.exp(-om**ga)
            omax = 750.0 ** (1.0 / ga)
            if f(self._cf) >= 0:
                omega_c = omax if f(omax) == f(self._cf) else self._cf
            else:
                omega_c = brentq(f, self._cf, omax)

        min_scale = omega_c / np.pi
        max_freq = 1.0 / (min_scale * (2.0 * np.pi / self._cf))
        if max_freq > fs / 2.0:
            max_freq = fs / 2.0
        return min_freq, max_freq

    def _resolve_frequency_limits(self, frequency_limits, fs: float):
        if frequency_limits is None:
            lo, hi = self._band_edges()
            return (lo, hi)
        lo, hi = (float(v) for v in frequency_limits)
        if not hi > lo:
            raise ValueError("frequency_limits must be increasing")
        if hi <= 0 or lo >= fs / 2.0:
            raise ValueError(f"frequency_limits must lie inside [0, {fs / 2.0}]")
        if lo < self._min_frequency():
            lo = self._min_frequency()
        if hi > fs / 2.0:
            hi = fs / 2.0
        if np.log2(hi) - np.log2(lo) < 1.0 / self.voices_per_octave:
            raise ValueError(
                "frequency_limits are too close for the requested voices per octave"
            )
        return (lo, hi)

    def _scales(self) -> np.ndarray:
        w_lo = 2.0 * np.pi * self.frequency_limits[0] / self.sampling_frequency
        w_hi = 2.0 * np.pi * self.frequency_limits[1] / self.sampling_frequency
        s0 = self._cf / w_hi
        smax = self._cf / w_lo
        n_octaves = np.log2(smax / s0)
        n = int(np.floor(self.voices_per_octave * n_octaves)) + 1
        return s0 * 2.0 ** (np.arange(n) / self.voices_per_octave)

    # -- use -----------------------------------------------------------------
    def pad_signal(self, x: np.ndarray) -> np.ndarray:
        """Apply the boundary extension used before the FFT."""
        if self.pad == 0:
            return np.asarray(x, dtype=float).ravel()
        x = np.asarray(x, dtype=float).ravel()
        if x.size != self.signal_length:
            raise ValueError(
                f"expected a signal of length {self.signal_length}, got {x.size}"
            )
        return np.concatenate((x[: self.pad][::-1], x, x[::-1][: self.pad]))

    def __repr__(self) -> str:  # pragma: no cover - cosmetic
        return (
            f"CWTFilterBank(signal_length={self.signal_length}, "
            f"wavelet={self.wavelet!r}, voices_per_octave={self.voices_per_octave}, "
            f"sampling_frequency={self.sampling_frequency:g}, "
            f"frequency_limits=({self.frequency_limits[0]:.6g}, "
            f"{self.frequency_limits[1]:.6g}), boundary={self.boundary!r}, "
            f"n_scales={self.scales.size})"
        )
