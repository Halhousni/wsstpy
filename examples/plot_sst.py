"""Plot a chirp's synchrosqueezed transform with its ridge.

Run with:  python examples/plot_sst.py
"""
import matplotlib.pyplot as plt
import numpy as np

import wsstpy

fs = 500.0
t = np.arange(0, 4, 1 / fs)
f0, f1 = 20.0, 120.0
phase = 2 * np.pi * (f0 * t + 0.5 * (f1 - f0) / t[-1] * t**2)
x = np.cos(phase) + 0.05 * np.random.default_rng(0).standard_normal(t.size)

sst, f = wsstpy.wsst(x, fs)
fridge, iridge = wsstpy.wsstridge(sst, penalty=0.5, f=f)
tone = wsstpy.iwsst(sst, ridge_indices=iridge)

fig, axes = plt.subplots(3, 1, figsize=(9, 8), sharex=True)

axes[0].plot(t, x, lw=0.6)
axes[0].set_ylabel("signal")

axes[1].pcolormesh(t, f, np.abs(sst), shading="auto", cmap="magma")
axes[1].plot(t, fridge[:, 0], lw=0.8, color="cyan")
axes[1].set_ylim(0, 200)
axes[1].set_ylabel("frequency (Hz)")

axes[2].plot(t, tone, lw=0.8)
axes[2].set_ylabel("ridge component")
axes[2].set_xlabel("time (s)")

fig.tight_layout()
fig.savefig("sst_chirp.png", dpi=150)
print("wrote sst_chirp.png")
