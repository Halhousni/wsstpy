"""Ridge extraction from a synchrosqueezed transform.

Ridges are the curves of locally maximal energy in the time-frequency plane.
They are found with dynamic programming: the total energy along a ridge is
maximised while moving between rows costs ``penalty`` times the square of the
row difference, which is the penalisation MATLAB's ``wsstridge`` uses.

The recurrence, for column ``n`` and row ``i``::

    D[n, i] = E[n, i] + min_j ( D[n-1, j] + penalty * (i - j)^2 )

where ``E = -log|sst|``.  Solving it exactly costs ``O(n_rows**2)`` per column,
which is small next to the transform itself.
"""

from __future__ import annotations

import numpy as np

__all__ = ["wsstridge"]

_REALMIN_LOG = -np.log(np.finfo(float).tiny)


def _ridge_once(cost: np.ndarray, penalty: float) -> np.ndarray:
    """Cheapest path through ``cost`` (n_rows x n_cols), as row indices."""
    n_rows, n_cols = cost.shape
    back = np.zeros((n_cols, n_rows), dtype=np.intp)
    total = cost[:, 0].copy()
    rows = np.arange(n_rows)
    dist = penalty * (rows[:, None] - rows[None, :]) ** 2
    for n in range(1, n_cols):
        cand = total[None, :] + dist
        best = np.argmin(cand, axis=1)
        back[n] = best
        total = cost[:, n] + cand[rows, best]
    path = np.empty(n_cols, dtype=np.intp)
    path[-1] = int(np.argmin(total))
    for n in range(n_cols - 1, 0, -1):
        path[n - 1] = back[n, path[n]]
    return path


def wsstridge(
    sst,
    penalty: float = 0.0,
    num_ridges: int = 1,
    num_frequency_bins: int = 4,
    f=None,
    filterbank=None,
):
    """Extract ridges from a synchrosqueezed transform.

    Parameters
    ----------
    sst : ndarray, shape (n_scales, n)
        Synchrosqueezed transform, as returned by :func:`wsstpy.wsst`.
    penalty : float, optional
        Cost of moving between rows, applied to the squared row difference.
        ``0`` lets the ridge jump freely, larger values give smoother curves.
    num_ridges : int, optional
        Number of ridges to extract.  Each ridge is removed from the cost
        before the next one is searched for.
    num_frequency_bins : int, optional
        Half-width, in rows, of the band blanked out around an extracted
        ridge.  Only used when ``num_ridges > 1``.
    f : ndarray, optional
        Frequencies of the rows of ``sst``.  When omitted they are reconstructed
        from ``filterbank`` (or from the default filter bank parameters).
    filterbank : CWTFilterBank, optional
        Filter bank used to produce ``sst``, needed to label the ridges when
        ``f`` is not given.

    Returns
    -------
    fridge : ndarray, shape (n, num_ridges)
        Frequency of each ridge at each time sample.
    iridge : ndarray, shape (n, num_ridges)
        Row index (1-based) of each ridge, ready to pass to
        :func:`wsstpy.iwsst`.

    Examples
    --------
    >>> import numpy as np, wsstpy
    >>> n = 256
    >>> x = np.cos(2 * np.pi * 0.1 * np.arange(n))
    >>> sst, f = wsstpy.wsst(x)
    >>> fridge, iridge = wsstpy.wsstridge(sst, f=f)
    >>> bool(np.allclose(fridge[:, 0], 0.1, atol=0.01))
    True
    """
    sst = np.asarray(sst)
    if sst.ndim != 2:
        raise ValueError("sst must be a 2-D array")
    n_rows, n_cols = sst.shape
    if num_ridges < 1:
        raise ValueError("num_ridges must be >= 1")
    if num_frequency_bins < 0:
        raise ValueError("num_frequency_bins must be >= 0")

    if f is None:
        if filterbank is None:
            from .filterbank import CWTFilterBank

            filterbank = CWTFilterBank(n_cols, 32, "amor", 1.0, (0.0, 0.5), "periodic")
        f = filterbank.center_frequencies[::-1]
    f = np.asarray(f, dtype=float).ravel()
    if f.size != n_rows:
        raise ValueError(f"f must have {n_rows} entries, got {f.size}")

    magnitude = np.abs(sst)
    scale = magnitude.sum()
    cost = -np.log(magnitude + 1e-8 * scale)

    iridge = np.empty((n_cols, num_ridges), dtype=np.intp)
    for k in range(num_ridges):
        path = _ridge_once(cost, float(penalty))
        iridge[:, k] = path + 1
        if k + 1 < num_ridges:
            for offset in range(num_frequency_bins + 1):
                lower = np.maximum(path - offset, 0)
                upper = np.minimum(path + offset, n_rows - 1)
                cost[lower, np.arange(n_cols)] = _REALMIN_LOG
                if offset > 0:
                    cost[upper, np.arange(n_cols)] = _REALMIN_LOG

    fridge = f[iridge - 1]
    return fridge, iridge
