"""Wavelet synchrosqueezing transform and its inverse.

``wsst`` reallocates the wavelet coefficients along the frequency axis using
the phase transform, which turns a smeared scalogram into sharp ridges.  The
implementation follows the algorithm of Daubechies, Lu and Wu (2011) with the
discretisation used by MATLAB's ``wsst``, so both give the same numbers.
"""

from __future__ import annotations

import numpy as np

from .filterbank import CWTFilterBank

__all__ = ["wsst", "iwsst"]


def wsst(
    x,
    fs: float | None = None,
    wavelet: str = "amor",
    voices_per_octave: int = 32,
    extend_signal: bool = False,
    wavelet_parameters=None,
    return_all: bool = False,
):
    """Wavelet synchrosqueezed transform of a real signal.

    Parameters
    ----------
    x : array_like, shape (n,)
        Real-valued signal, at least 4 samples.
    fs : float, optional
        Sampling frequency.  When omitted the frequencies are returned in
        cycles/sample.
    wavelet : {"amor", "morse", "bump"}, optional
        Analytic wavelet, ``"amor"`` (analytic Morlet) by default.
    voices_per_octave : int, optional
        Even integer, 32 by default.
    extend_signal : bool, optional
        Mirror the signal at both ends before the transform (MATLAB's
        ``ExtendSignal``).  This reduces edge artefacts at the cost of
        assuming the signal is symmetric outside its support.
    wavelet_parameters : (float, float), optional
        ``(gamma, time_bandwidth)`` for the Morse wavelet.
    return_all : bool, optional
        Also return the filter bank, which gives access to the scales, the
        wavelet filters and the frequencies of the rows.

    Returns
    -------
    sst : ndarray, shape (n_scales, n)
        Complex synchrosqueezed transform, rows ordered from the lowest to
        the highest frequency.
    f : ndarray, shape (n_scales,)
        Frequency of each row, in the unit of ``fs``.
    fb : CWTFilterBank
        Only returned when ``return_all`` is true.

    Examples
    --------
    >>> import numpy as np
    >>> t = np.arange(1024) / 1024
    >>> x = np.cos(2 * np.pi * 0.1 * np.arange(1024))
    >>> sst, f = wsst(x)
    >>> bool(abs(f[np.argmax(abs(sst).sum(axis=1))] - 0.1) < 0.005)
    True
    """
    x = np.asarray(x, dtype=float).ravel()
    n = x.size
    if n < 4:
        raise ValueError("x must have at least 4 samples")
    if not np.all(np.isfinite(x)):
        raise ValueError("x must be finite")
    fs = 1.0 if fs is None else float(fs)
    if not fs > 0:
        raise ValueError("fs must be positive")

    fb = CWTFilterBank(
        n,
        voices_per_octave=voices_per_octave,
        wavelet=wavelet,
        sampling_frequency=fs,
        frequency_limits=(0.0, fs / 2.0),
        boundary="reflection" if extend_signal else "periodic",
        wavelet_parameters=wavelet_parameters,
    )

    threshold = np.sqrt(np.finfo(float).eps)
    vpo = fb.voices_per_octave

    xin = x - x.mean()
    xdft = np.fft.fft(fb.pad_signal(xin))

    cfs = np.fft.ifft(fb.psi_dft * xdft[None, :], axis=1)
    dcfs = np.fft.ifft(1j * fb.omega[None, :] * fb.psi_dft * xdft[None, :], axis=1)
    if fb.pad:
        sl = slice(fb.pad, fb.pad + n)
        cfs = cfs[:, sl]
        dcfs = dcfs[:, sl]

    # Rows go from the coarse to the fine scales so that frequencies increase.
    cfs = cfs[::-1, :] / 2.0
    dcfs = dcfs[::-1, :] / 2.0
    f = fb.center_frequencies[::-1].copy()

    with np.errstate(divide="ignore", invalid="ignore"):
        phase = np.imag(dcfs / cfs) / (2.0 * np.pi)
    phase[~np.isfinite(phase) | (phase < 0.0)] = 0.0

    # the phase transform is in cycles/sample, so the row grid is too
    sst = _squeeze(cfs, phase, threshold, np.log2(f[0] / fs), vpo)
    if return_all:
        return sst, f, fb
    return sst, f


def _round_half_away(x: np.ndarray) -> np.ndarray:
    """MATLAB's ``round``: halves go away from zero."""
    return np.where(x >= 0.0, np.floor(x + 0.5), np.ceil(x - 0.5))


def _squeeze(cfs, phase, threshold, log2_fmin, vpo):
    """Accumulate the wavelet coefficients into the synchrosqueezed rows."""
    n_scales, n = cfs.shape
    with np.errstate(divide="ignore", invalid="ignore"):
        row = 1.0 + vpo * (np.log2(phase) - log2_fmin)
    row = _round_half_away(row)
    row = np.clip(np.nan_to_num(row, nan=0.0, posinf=n_scales + 1.0,
                                neginf=-1.0), 0, n_scales)

    keep = (np.abs(cfs) > threshold) & (row >= 1) & (row <= n_scales)
    if not keep.any():
        return np.zeros((n_scales, n), dtype=complex)

    rows, cols = np.nonzero(keep)
    sst = np.zeros((n_scales, n), dtype=complex)
    np.add.at(
        sst,
        (row[keep].astype(np.intp) - 1, np.nonzero(keep)[1]),
        cfs[keep] * np.log(2.0) / vpo,
    )
    return sst


def _wavelet_reconstruction_constant(wavelet: str, wavelet_parameters=None) -> float:
    """``R_psi`` used by the inverse transform, ``int_0^inf psi(w)/w dw``."""
    import warnings

    from scipy.integrate import IntegrationWarning, quad

    name = str(wavelet).lower()
    if name in ("amor", "morlet", "amorlet"):
        f = lambda w: np.exp(-0.5 * (w - 6.0) ** 2) / w
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", IntegrationWarning)
            val, _ = quad(f, np.sqrt(np.finfo(float).eps), np.inf, limit=800,
                          epsabs=1e-15, epsrel=1e-15)
        return float(val)
    if name == "bump":
        mu, sigma = 5.0, 1.0
        f = lambda w: np.e * np.exp(-1.0 / (1.0 - ((w - mu) / sigma) ** 2)) / w
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", IntegrationWarning)
            val, _ = quad(f, mu - sigma, mu + sigma, limit=800,
                          epsabs=1e-15, epsrel=1e-15)
        return float(val)
    if name == "morse":
        from scipy.special import gammaln

        from ._wavelet import morse_norm_constant

        ga, tb = (3.0, 60.0) if wavelet_parameters is None else wavelet_parameters
        beta = float(tb) / float(ga)
        anorm = morse_norm_constant(ga, beta)
        return float(anorm / (2.0 * ga) * np.exp(gammaln(beta / ga)))
    raise ValueError(f"unknown wavelet {wavelet!r}")


def iwsst(
    sst,
    f=None,
    freq_range=None,
    wavelet: str = "amor",
    wavelet_parameters=None,
    ridge_indices=None,
    num_frequency_bins: int = 16,
):
    """Inverse wavelet synchrosqueezed transform.

    Parameters
    ----------
    sst : ndarray, shape (n_scales, n)
        Synchrosqueezed transform, as returned by :func:`wsst`.
    f : ndarray, optional
        Frequencies of the rows.  Required together with ``freq_range`` to
        reconstruct a band-limited signal.
    freq_range : (float, float), optional
        Band to reconstruct, in the same unit as ``f``.
    wavelet : {"amor", "morse", "bump"}, optional
        Wavelet used to build ``sst``.
    wavelet_parameters : (float, float), optional
        Morse wavelet parameters.
    ridge_indices : ndarray, shape (n, n_ridges), optional
        Row indices of ridges to reconstruct separately, as returned by
        :func:`wsstridge`.  The reconstruction keeps ``num_frequency_bins``
        bins on either side of each ridge.
    num_frequency_bins : int, optional
        Half-width, in rows, of the band kept around each ridge.

    Returns
    -------
    ndarray, shape (n,) or (n, n_ridges)
        Reconstructed signal, or one reconstruction per ridge.
    """
    sst = np.asarray(sst)
    if sst.ndim != 2:
        raise ValueError("sst must be a 2-D array")
    if not np.all(np.isfinite(sst)):
        raise ValueError("sst must be finite")
    n_scales, n = sst.shape
    rpsi = _wavelet_reconstruction_constant(wavelet, wavelet_parameters)

    if ridge_indices is None:
        if freq_range is None:
            return (2.0 / rpsi * np.real(sst.sum(axis=0))).ravel()
        if f is None:
            raise ValueError("f is required when freq_range is given")
        f = np.asarray(f, dtype=float).ravel()
        low = np.searchsorted(f, freq_range[0], side="left")
        high = np.searchsorted(f, freq_range[1], side="right")
        tx = np.zeros_like(sst)
        tx[low:high, :] = sst[low:high, :]
        return (2.0 / rpsi * np.real(tx.sum(axis=0))).ravel()

    idx = np.asarray(ridge_indices)
    if idx.ndim != 2 or idx.shape[0] != n:
        raise ValueError("ridge_indices must have shape (n, n_ridges)")
    nbins = int(num_frequency_bins)
    if nbins < 2:
        raise ValueError("num_frequency_bins must be >= 2")
    out = np.empty((n, idx.shape[1]))
    for k in range(idx.shape[1]):
        c = idx[:, k].astype(np.intp) - 1
        mask = np.zeros((n_scales, n), dtype=bool)
        for j in range(n):
            lo = max(0, c[j] - nbins)
            hi = min(n_scales - 1, c[j] + nbins)
            mask[lo : hi + 1, j] = True
        out[:, k] = 2.0 / rpsi * np.real(np.where(mask, sst, 0).sum(axis=0))
    if out.shape[1] == 1:
        return out[:, 0].copy()
    return out
