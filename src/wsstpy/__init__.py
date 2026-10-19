"""wsstpy -- wavelet synchrosqueezing in Python.

A self-contained NumPy implementation of the wavelet synchrosqueezed
transform: forward and inverse transforms, the underlying continuous wavelet
transform and filter bank, and ridge extraction.
"""

from __future__ import annotations

from ._wavelet import SUPPORTED_WAVELETS
from .cwt import adm_constant, cwt, icwt
from .filterbank import CWTFilterBank, omega_grid
from .ridge import wsstridge
from .wsst import iwsst, wsst

__version__ = "0.1.0"

__all__ = [
    "wsst",
    "iwsst",
    "wsstridge",
    "cwt",
    "icwt",
    "adm_constant",
    "CWTFilterBank",
    "omega_grid",
    "SUPPORTED_WAVELETS",
    "__version__",
]
