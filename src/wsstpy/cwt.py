"""Continuous wavelet transform and its inverse.

``cwt`` matches MATLAB's ``cwt``: the same filter bank, the same default band
and the same reconstruction constant, so coefficients agree to floating point
precision.
"""

from __future__ import annotations

import numpy as np

from ._wavelet import morse_norm_constant, wavelet_parameters
from .filterbank import CWTFilterBank

__all__ = ["cwt", "icwt", "adm_constant"]


def cwt(
    x,
    fs: float | None = None,
    wavelet: str = "morse",
    voices_per_octave: int = 10,
    extend_signal: bool = True,
    frequency_limits=None,
    wavelet_parameters=(3.0, 60.0),
    return_all: bool = False,
):
    """Continuous wavelet transform using an analytic wavelet.

    Parameters
    ----------
    x : array_like, shape (n,)
        Real or complex signal, 4 samples or longer.
    fs : float, optional
        Sampling frequency.  Frequencies are reported in cycles/sample when
        omitted.
    wavelet : {"morse", "amor", "bump"}, optional
        Analysis wavelet, the Morse wavelet by default as in MATLAB.
    voices_per_octave : int, optional
        Number of scales per octave, 10 by default.
    extend_signal : bool, optional
        Mirror the signal at both ends (default true) before the transform.
    frequency_limits : (float, float), optional
        Band to analyse in Hz.  The default depends on the wavelet and the
        signal length.
    wavelet_parameters : (float, float), optional
        ``(gamma, time_bandwidth)`` for the Morse wavelet, ``(3, 60)`` by
        default.
    return_all : bool, optional
        Also return the filter bank.

    Returns
    -------
    cfs : ndarray, shape (n_scales, n)
        Complex wavelet coefficients, finest scale first.
    f : ndarray, shape (n_scales,)
        Frequency of each row.
    fb : CWTFilterBank
        Only when ``return_all`` is true.
    """
    x = np.asarray(x)
    if x.ndim != 1:
        x = x.ravel()
    n = x.size
    if n < 4:
        raise ValueError("x needs 4 samples or more")
    if not np.all(np.isfinite(x)):
        raise ValueError("x must be finite")
    fs = 1.0 if fs is None else float(fs)

    fb = CWTFilterBank(
        n,
        voices_per_octave=voices_per_octave,
        wavelet=wavelet,
        sampling_frequency=fs,
        frequency_limits=frequency_limits,
        boundary="reflection" if extend_signal else "periodic",
        wavelet_parameters=wavelet_parameters,
    )

    if np.iscomplexobj(x):
        xh = x / 2.0
        xv = fb.pad_signal(xh)
        pos = np.fft.ifft(fb.psi_dft * np.fft.fft(xv)[None, :], axis=1)
        neg = np.conj(
            np.fft.ifft(fb.psi_dft * np.fft.fft(np.conj(xv))[None, :], axis=1)
        )
        cfs = np.stack((pos, neg), axis=2)
    else:
        cfs = np.fft.ifft(fb.psi_dft * np.fft.fft(fb.pad_signal(x))[None, :], axis=1)
    if fb.pad:
        cfs = cfs[:, fb.pad : fb.pad + n]

    f = fb.center_frequencies.copy()
    if return_all:
        return cfs, f, fb
    return cfs, f


def adm_constant(wavelet: str = "morse", wavelet_parameters=(3.0, 60.0)) -> float:
    """Admissibility constant ``int_0^inf |psi(w)| / w dw``."""
    import warnings

    from scipy.integrate import IntegrationWarning, quad
    from scipy.special import gammaln

    name = str(wavelet).lower()
    if name in ("amor", "morlet", "amorlet"):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", IntegrationWarning)
            val, _ = quad(lambda w: abs(2.0 * np.exp(-0.5 * (w - 6.0) ** 2)) / w,
                          np.sqrt(np.finfo(float).eps), np.inf, limit=800,
                          epsabs=1e-15, epsrel=1e-15)
        return float(val)
    if name == "bump":
        mu, sigma = 5.0, 0.6
        lo = mu - sigma
        hi = mu + sigma
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", IntegrationWarning)
            val, _ = quad(
                lambda w: abs(2.0 * np.e * np.exp(-1.0 / (1.0 - ((w - mu) / sigma) ** 2))) / w,
                lo, hi, limit=800, epsabs=1e-15, epsrel=1e-15)
        return float(val)
    if name == "morse":
        ga, tb = wavelet_parameters
        beta = float(tb) / float(ga)
        return float(morse_norm_constant(ga, beta) / ga * np.exp(gammaln(beta / ga)))
    raise ValueError(f"unknown wavelet {wavelet!r}")


def icwt(
    cfs,
    fs: float | None = None,
    wavelet: str = "morse",
    voices_per_octave: int | None = None,
    f=None,
    freq_range=None,
    wavelet_parameters=(3.0, 60.0),
    signal_mean=None,
):
    """Inverse continuous wavelet transform.

    Parameters
    ----------
    cfs : ndarray, shape (n_scales, n)
        Wavelet coefficients from :func:`cwt` (real or complex input).
    fs : float, optional
        Sampling frequency, only used when ``f`` must be rebuilt.
    wavelet : {"morse", "amor", "bump"}, optional
        Wavelet used for the analysis.
    voices_per_octave : int, optional
        Voices per octave used for the analysis.  When omitted the scale ratio
        is read from ``f`` (or assumed to be 10).
    f : ndarray, optional
        Frequencies of the rows of ``cfs``; used to obtain the scale ratio and
        to select ``freq_range``.
    freq_range : (float, float), optional
        Band to reconstruct in the unit of ``f``.
    wavelet_parameters : (float, float), optional
        Morse wavelet parameters.
    signal_mean : float, optional
        Value added back after reconstruction.

    Returns
    -------
    ndarray, shape (n,)
        Reconstructed signal.
    """
    cfs = np.asarray(cfs)
    if cfs.ndim not in (2, 3):
        raise ValueError("cfs must be 2-D (or 3-D for complex input)")
    if not np.all(np.isfinite(cfs)):
        raise ValueError("cfs must be finite")
    n_scales, n = cfs.shape[:2]

    if voices_per_octave is not None:
        if f is not None:
            raise ValueError("pass either voices_per_octave or f, not both")
        ds = 1.0 / float(voices_per_octave)
    elif f is not None:
        f = np.asarray(f, dtype=float).ravel()
        if f.size != n_scales:
            raise ValueError(f"f must have {n_scales} entries, got {f.size}")
        ds = float(np.mean(np.diff(np.log2(1.0 / f))))
    else:
        ds = 0.1

    wt = cfs
    if f is not None and freq_range is not None:
        f = np.asarray(f, dtype=float).ravel()
        # ``cwt`` orders rows from the highest frequency down, so the band is
        # the contiguous run between the last row above the band and the first
        # row below it (same rule as MATLAB's findFreqIndices).
        above = np.nonzero(f >= freq_range[1])[0]
        below = np.nonzero(f <= freq_range[0])[0]
        keep = np.zeros(n_scales, dtype=bool)
        if above.size and below.size:
            keep[above.max() : below.min() + 1] = True
        wt = cfs.copy()
        wt[~keep] = 0

    cpsi = adm_constant(wavelet, wavelet_parameters)
    a0 = 2.0 ** ds
    if wt.ndim == 3:
        xrec = 2.0 * np.log(a0) / cpsi * np.real(wt[:, :, 0] + wt[:, :, 1]).sum(axis=0)
    else:
        xrec = 2.0 * np.log(a0) / cpsi * np.real(wt).sum(axis=0)
    if signal_mean is not None:
        xrec = xrec + float(signal_mean)
    return np.asarray(xrec).ravel()
