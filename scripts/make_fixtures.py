"""Build the committed test fixtures from a MATLAB reference run.

The arrays come from ``probe/probe7.m`` (MATLAB R2025b, Wavelet Toolbox) and
are stored in double precision so the tests can assert agreement at 1e-12.

    python scripts/make_fixtures.py
"""
from __future__ import annotations

import os

import numpy as np
import scipy.io as sio

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
SRC = os.path.join(ROOT, "probe", "data7", "ref7.mat")
OUT = os.path.join(ROOT, "tests", "data", "reference.npz")

ref = sio.loadmat(SRC, squeeze_me=True, struct_as_record=False)
S, W, C, I, G = ref["S"], ref["W"], ref["C"], ref["I"], ref["G"]

out = {}
for name in ("tone", "chirp", "noise", "short", "impulse", "odd"):
    out[f"sig_{name}"] = np.asarray(getattr(S, name), dtype=np.float64).ravel()
    out[f"wsst_{name}"] = np.asarray(getattr(W, f"sst_{name}"), dtype=np.complex128)
    out[f"f_{name}"] = np.asarray(getattr(W, f"f_{name}"), dtype=np.float64).ravel()

for tag in ("vpo10", "extend", "morse", "morse30", "bump", "fs"):
    out[f"wsst_{tag}"] = np.asarray(getattr(W, f"sst_{tag}"), dtype=np.complex128)
    out[f"f_{tag}"] = np.asarray(getattr(W, f"f_{tag}"), dtype=np.float64).ravel()

out["cwt_morse"] = np.asarray(C.cfs_morse, dtype=np.complex128)
out["f_cwt_morse"] = np.asarray(C.f_morse, dtype=np.float64).ravel()
out["cwt_amor"] = np.asarray(C.cfs_amor, dtype=np.complex128)
out["cwt_limits"] = np.asarray(C.cfs_lim, dtype=np.complex128)
out["f_cwt_limits"] = np.asarray(C.f_lim, dtype=np.float64).ravel()
out["icwt_full"] = np.asarray(C.icwt_full, dtype=np.float64).ravel()
out["icwt_band"] = np.asarray(C.icwt_band, dtype=np.float64).ravel()
out["icwt_vpo10"] = np.asarray(C.icwt_vpo, dtype=np.float64).ravel()
out["icwt_amor"] = np.asarray(C.icwt_amor, dtype=np.float64).ravel()

out["iwsst_full"] = np.asarray(I.full, dtype=np.float64).ravel()
out["iwsst_band"] = np.asarray(I.band, dtype=np.float64).ravel()
out["iwsst_morse"] = np.asarray(I.morse, dtype=np.float64).ravel()
out["iwsst_ridge"] = np.asarray(I.ridge, dtype=np.float64)

for pen in (0, 1, 5, 100):
    g = getattr(G, f"pen{pen}")
    out[f"iridge_pen{pen}"] = np.asarray(g.iridge, dtype=np.int64)
    out[f"fridge_pen{pen}"] = np.asarray(g.fridge, dtype=np.float64).ravel()
out["iridge_two"] = np.asarray(G.r2.iridge, dtype=np.int64)
out["iridge_three_bins6"] = np.asarray(G.r3.iridge, dtype=np.int64)
out["iridge_chirp"] = np.asarray(G.chirp.iridge, dtype=np.int64)
out["iridge_noise"] = np.asarray(G.noise.iridge, dtype=np.int64)

os.makedirs(os.path.dirname(OUT), exist_ok=True)
np.savez_compressed(OUT, **out)
print(f"wrote {OUT} with {len(out)} arrays, {os.path.getsize(OUT) / 1e6:.2f} MB")
