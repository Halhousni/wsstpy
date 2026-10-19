"""Wavelet definitions used by the filter bank.

The conventions here mirror those of the MATLAB Wavelet Toolbox so that the
numbers produced by :func:`wsstpy.wsst` match ``wsst`` to floating point
precision.  All wavelets are analytic and are defined through their Fourier
transform, which peaks at a value of 2 (unit L1 norm).

Notation
--------
``s_omega`` is the product of scale and angular frequency in radians/sample.
"""

from __future__ import annotations

import numpy as np
from scipy.special import gammaln

__all__ = ["SUPPORTED_WAVELETS", "wavelet_parameters", "frequency_response"]

SUPPORTED_WAVELETS = ("amor", "morse", "bump")

#: Measured time-domain standard deviation of the bump wavelet (MATLAB value).
_BUMP_SIGMA_T = 5.847705
_BUMP_CF = 5.0
_BUMP_SIGMA = 0.6
_AMOR_CF = 6.0
_AMOR_SIGMA_T = np.sqrt(2.0)


def _morse_sigma_t(gamma: float, beta: float) -> float:
    """Time-domain standard deviation of the Morse wavelet.

    Evaluated in closed form (all integrals are Gamma functions) using the
    same normalisation as MATLAB.  ``beta`` is the time-bandwidth product
    divided by ``gamma``.
    """
    ga, be = float(gamma), float(beta)

    def _int_ratio(p):
        # integral of om**p * exp(-2*om**ga) from 0 to inf
        return 2.0 ** (-(p + 1.0) / ga) * np.exp(gammaln((p + 1.0) / ga)) / ga

    num = (
        be**2 * _int_ratio(2 * be - 2)
        - 2 * be * ga * _int_ratio(2 * be + ga - 2)
        + ga**2 * _int_ratio(2 * be + 2 * ga - 2)
    )
    den = _int_ratio(2 * be)
    cf = (be / ga) ** (1.0 / ga)
    return float(np.sqrt(cf**2 * num / den))


def morse_peak_frequency(gamma: float, beta: float) -> float:
    """Peak angular frequency of the Morse wavelet."""
    return (float(beta) / float(gamma)) ** (1.0 / float(gamma))


def morse_norm_constant(gamma: float, beta: float) -> float:
    """Peak value of ``om**beta * exp(-om**gamma)`` scaled to peak magnitude 2."""
    ga, be = float(gamma), float(beta)
    return 2.0 * np.exp(be / ga * (1.0 + np.log(ga) - np.log(be)))


def wavelet_parameters(wavelet: str, wavelet_parameters=None):
    """Return ``(cf, sigma_t, gamma, beta)`` for the requested wavelet.

    ``cf`` is the peak angular frequency in radians/sample, ``sigma_t`` the
    time-domain standard deviation (used to pick the coarsest scale).
    ``gamma``/``beta`` are only meaningful for the Morse wavelet.
    """
    name = str(wavelet).lower()
    if name in ("amor", "morlet", "amorlet"):
        return _AMOR_CF, _AMOR_SIGMA_T, None, None
    if name == "bump":
        return _BUMP_CF, _BUMP_SIGMA_T, None, None
    if name == "morse":
        if wavelet_parameters is None:
            ga, tb = 3.0, 60.0
        else:
            ga, tb = (float(v) for v in wavelet_parameters)
            if ga < 1:
                raise ValueError("Morse symmetry parameter (gamma) must be >= 1")
            if tb < ga:
                raise ValueError("Morse time-bandwidth parameter must be >= gamma")
            if tb / ga > 40:
                raise ValueError("Morse beta = time-bandwidth/gamma must be <= 40")
        beta = tb / ga
        return morse_peak_frequency(ga, beta), _morse_sigma_t(ga, beta), ga, beta
    raise ValueError(
        f"unknown wavelet {wavelet!r}; expected one of {SUPPORTED_WAVELETS}"
    )


_params = wavelet_parameters  # module-level alias, the argument below shadows the name


def frequency_response(wavelet: str, s_omega: np.ndarray, wavelet_parameters=None):
    """Evaluate the wavelet's Fourier transform at ``s_omega``.

    Parameters
    ----------
    wavelet : {"amor", "morse", "bump"}
    s_omega : ndarray
        Scale times angular frequency, shape ``(n_scales, n_freq)``.
    wavelet_parameters : sequence of 2 floats, optional
        ``(gamma, time_bandwidth)`` for the Morse wavelet.

    Returns
    -------
    ndarray of the same shape as ``s_omega``, real valued.
    """
    name = str(wavelet).lower()
    so = np.asarray(s_omega, dtype=float)
    positive = so > 0.0

    if name in ("amor", "morlet", "amorlet"):
        with np.errstate(over="ignore"):
            psi = 2.0 * np.exp(-0.5 * (so - _AMOR_CF) ** 2)
        return np.where(positive, psi, 0.0)

    if name == "morse":
        _, _, ga, beta = _params(name, wavelet_parameters)
        fo = morse_peak_frequency(ga, beta)
        factor = np.exp(-beta * np.log(fo) + fo**ga)
        with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
            psi = 2.0 * factor * np.exp(beta * np.log(np.abs(so)) - np.abs(so) ** ga)
        return np.where(positive, psi, 0.0)

    if name == "bump":
        w = (so - _BUMP_CF) / _BUMP_SIGMA
        inside = np.abs(w) < 1.0 - np.finfo(float).eps
        with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
            psi = 2.0 * np.e * np.exp(-1.0 / (1.0 - w**2))
        psi = np.where(np.isfinite(psi), psi, 0.0)
        return np.where(inside, psi, 0.0)

    raise ValueError(
        f"unknown wavelet {wavelet!r}; expected one of {SUPPORTED_WAVELETS}"
    )
