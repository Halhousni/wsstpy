"""Make the validation figures: MATLAB output next to wsstpy output.

Reads the double-precision references in probe/data5 and probe/data6, runs
wsstpy on the same signals, and writes figures/sst_comparison.{png,pdf} and
figures/residuals_and_ridges.{png,pdf}.

    python scripts/make_validation_figure.py
"""
from __future__ import annotations

import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import scipy.io as sio
from matplotlib.colors import LogNorm

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "src"))
import wsstpy

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
FIG = os.path.join(ROOT, "figures")
os.makedirs(FIG, exist_ok=True)

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Nimbus Sans", "Helvetica", "Liberation Sans", "DejaVu Sans"],
    "font.size": 8,
    "axes.linewidth": 0.6,
    "xtick.major.width": 0.6,
    "ytick.major.width": 0.6,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "legend.frameon": False,
    "savefig.dpi": 200,
})

ref6 = sio.loadmat(os.path.join(ROOT, "probe", "data6", "ref6.mat"),
                   squeeze_me=True, struct_as_record=False)
ref5 = sio.loadmat(os.path.join(ROOT, "probe", "data5", "ref5.mat"),
                   squeeze_me=True, struct_as_record=False)
S6, W6 = ref6["S"], ref6["W"]
S5, W5 = ref5["S"], ref5["W"]

FS_EPS = np.finfo(float).eps

cases = []
for tag, name in [("chirp", "chirp"), ("noise", "random noise")]:
    x = np.asarray(getattr(S6, f"{tag}_1024"), dtype=float).ravel()
    matlab = np.asarray(getattr(W6, f"sst_{tag}_1024"))
    f = np.asarray(getattr(W6, f"f_{tag}_1024")).ravel()
    python = wsstpy.wsst(x)[0]
    cases.append(dict(tag=tag, label=name, x=x, matlab=matlab, python=python, f=f))

tone_x = np.asarray(S5.tone_1024).ravel()
tone_m = np.asarray(W5.sst_1024)
tone_p = wsstpy.wsst(tone_x)[0]

# ---------------------------------------------------------------- figure 1
fig, axes = plt.subplots(2, 3, figsize=(7.2, 4.4), constrained_layout=True)

for row, case in enumerate(cases):
    m, p, f = case["matlab"], case["python"], case["f"]
    diff = np.abs(p - m)
    vmax = np.abs(m).max() / 3.0
    for col, (data, title) in enumerate([(np.abs(m), "MATLAB"),
                                         (np.abs(p), "wsstpy")]):
        ax = axes[row, col]
        ax.pcolormesh(np.arange(m.shape[1]), f, data, shading="auto",
                      cmap="magma", vmin=0, vmax=vmax, rasterized=True)
        ax.set_title(title, fontsize=8)
        ax.set_ylabel("frequency (cycles/sample)")
        if row == 1:
            ax.set_xlabel("sample")
    ax = axes[row, 2]
    im = ax.pcolormesh(np.arange(m.shape[1]), f, np.maximum(diff, 1e-18),
                       shading="auto", cmap="inferno",
                       norm=LogNorm(vmin=1e-17, vmax=1e-14), rasterized=True)
    ax.set_title("absolute difference", fontsize=8)
    ax.set_xlabel("sample")
    ax.text(0.03, 0.03, f"max |Δ| = {diff.max():.1e}",
            transform=ax.transAxes, color="white", fontsize=7)
    cb = fig.colorbar(im, ax=ax, pad=0.02, fraction=0.046)
    cb.set_label("|wsstpy − MATLAB|", fontsize=7)
    cb.ax.tick_params(labelsize=6)
    axes[row, 0].text(-0.30, 0.5, case["label"], transform=axes[row, 0].transAxes,
                      rotation=90, va="center", ha="center", fontsize=8)

fig.savefig(os.path.join(FIG, "sst_comparison.png"), bbox_inches="tight")
fig.savefig(os.path.join(FIG, "sst_comparison.pdf"), bbox_inches="tight")
plt.close(fig)

# ---------------------------------------------------------------- figure 2
fig, axes = plt.subplots(1, 3, figsize=(7.2, 2.5), constrained_layout=True)

# (a) ridges: MATLAB and wsstpy on top of each other
ax = axes[0]
xc = np.asarray(S6.chirp_1024).ravel()
sst_p, f_c = wsstpy.wsst(xc)
fr_p, ir_p = wsstpy.wsstridge(sst_p, penalty=2.0, f=f_c)
ir_m = np.asarray(ref6["G"].chirp2.iridge).ravel()
t = np.arange(ir_m.size)
ax.plot(t, ir_m, lw=3.0, color="0.75", label="MATLAB")
ax.plot(t, ir_p[:, 0], lw=0.8, ls="--", color="#b2182b", label="wsstpy")
ax.set_xlabel("sample")
ax.set_ylabel("ridge row index")
ax.set_title("ridge extraction, chirp", fontsize=8)
ax.legend(loc="upper left", fontsize=7)
ax.text(0.45, 0.06, f"{int(np.sum(ir_p[:, 0] == ir_m))}/{ir_m.size} rows identical, "
                    f"penalty = 2",
        transform=ax.transAxes, fontsize=7)

# (b) cumulative distribution of the residual, in units of eps of the peak
ax = axes[1]
for case, colour in zip(cases, ["#2166ac", "#762a83"]):
    m = np.abs(case["matlab"])
    d = np.abs(case["python"] - case["matlab"])
    keep = m > 0
    scale = np.abs(case["matlab"]).max() * FS_EPS
    vals = np.sort(d[keep] / scale)
    frac = np.arange(1, vals.size + 1) / vals.size
    ax.semilogx(vals, frac, lw=1.0, color=colour, label=case["label"])
    median = vals[vals.size // 2]
    ax.plot([median], [0.5], marker="o", ms=3, color=colour)
ax.set_xlabel("|Δ| / (max|SST| · ε)")
ax.set_ylabel("fraction of coefficients")
ax.set_title("cumulative residual", fontsize=8)
ax.set_ylim(0, 1)
ax.set_xlim(1e-4, 2e2)
ax.axvline(1.0, color="0.4", lw=0.8, ls=":")
ax.annotate("1 ε", xy=(1.0, 0.08), xycoords=("data", "axes fraction"),
            xytext=(4, 0), textcoords="offset points", fontsize=7, color="0.35")
ax.legend(loc="lower left", fontsize=7)

# (c) relative difference against coefficient magnitude
ax = axes[2]
for case, colour in zip(cases, ["#2166ac", "#762a83"]):
    m = np.abs(case["matlab"]).ravel()
    d = np.abs(case["python"] - case["matlab"]).ravel()
    keep = m > 1e-6
    ax.scatter(m[keep], (d[keep] / m[keep]), s=1.5, alpha=0.35, color=colour,
               linewidths=0, label=case["label"], rasterized=True)
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlabel("|SST| (MATLAB)")
ax.set_ylabel("|Δ| / |SST|")
ax.set_title("relative difference", fontsize=8)
ax.axhline(FS_EPS, color="0.4", lw=0.8, ls=":")
ax.legend(loc="lower right", fontsize=7, markerscale=4)

fig.savefig(os.path.join(FIG, "residuals_and_ridges.png"), bbox_inches="tight")
fig.savefig(os.path.join(FIG, "residuals_and_ridges.pdf"), bbox_inches="tight")
plt.close(fig)

# ---------------------------------------------------------------- numbers
print("tone  N=1024: max|Δ| =", np.abs(tone_p - tone_m).max())
for case in cases:
    d = np.abs(case["python"] - case["matlab"])
    rel = d[np.abs(case["matlab"]) > 1e-6].max() / np.abs(case["matlab"]).max()
    print(f"{case['tag']:6s} N=1024: max|Δ| = {d.max():.3e}, peak-rel = {rel:.2e}")
print("ridge identical:", int(np.sum(ir_p[:, 0] == ir_m)), "/", ir_m.size)
print("wrote", FIG)
