"""Separate two tones with the synchrosqueezed transform.

Run with:  python examples/separate_tones.py
"""
import numpy as np

import wsstpy

fs = 250.0
t = np.arange(0, 8, 1 / fs)
slow = np.cos(2 * np.pi * 10 * t)
fast = 0.6 * np.cos(2 * np.pi * 45 * t + 0.3 * np.sin(2 * np.pi * 0.5 * t))
x = slow + fast

sst, f = wsstpy.wsst(x, fs)
fridge, iridge = wsstpy.wsstridge(sst, penalty=1.0, num_ridges=2, f=f)
parts = wsstpy.iwsst(sst, ridge_indices=iridge)

order = np.argsort(np.median(fridge, axis=0))
for k, idx in enumerate(order):
    err = np.abs(parts[:, idx] - (slow if k == 0 else fast)).max()
    print(f"ridge at {np.median(fridge[:, idx]):6.2f} Hz, "
          f"max error against the true component: {err:.3f}")
