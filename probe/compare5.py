"""Compare wsstpy against double-precision MATLAB R2025b references.

Writes a markdown table to stdout and JSON to probe/results5.json.
"""
import json
import os
import sys

import numpy as np
import scipy.io as sio

sys.path.insert(0, "/home/hilal/dev/wsstpy/src")
import wsstpy

DATA = "/home/hilal/dev/wsstpy/probe/data5"
ref = sio.loadmat(os.path.join(DATA, "ref5.mat"), squeeze_me=True, struct_as_record=False)

S = ref["S"]
W = ref["W"]
C = ref["C"]
I = ref["I"]
G = ref["G"]

rows = []


def cmp_matrix(name, mine, theirs):
    mine = np.asarray(mine)
    theirs = np.asarray(theirs)
    if mine.shape != theirs.shape:
        rows.append((name, f"SHAPE {mine.shape} vs {theirs.shape}", "", ""))
        return
    d = np.abs(mine - theirs)
    scale = np.abs(theirs).max()
    rows.append((name, f"{d.max():.3e}", f"{(d.max() / scale if scale else 0):.3e}",
                 f"{np.sqrt(np.mean(d**2)):.3e}"))


def cmp_vector(name, mine, theirs):
    cmp_matrix(name, np.asarray(mine).ravel(), np.asarray(theirs).ravel())


def cmp_exact(name, mine, theirs):
    mine = np.asarray(mine).ravel()
    theirs = np.asarray(theirs).ravel()
    n_same = int(np.sum(mine == theirs))
    rows.append((name, f"{n_same}/{mine.size} identical",
                 "0" if n_same == mine.size else f"{np.abs(mine - theirs).max()} rows off", ""))


Ns = [256, 512, 1000, 1024, 1025, 4096]

print("running wsst comparisons ...")
for N in Ns:
    x = getattr(S, f"tone_{N}")
    sst, f = wsstpy.wsst(x)
    cmp_matrix(f"wsst tone N={N}", sst, getattr(W, f"sst_{N}"))
    cmp_vector(f"  f tone N={N}", f, getattr(W, f"f_{N}"))
    sst10, f10 = wsstpy.wsst(x, voices_per_octave=10)
    cmp_matrix(f"wsst vpo=10 N={N}", sst10, getattr(W, f"sst_vpo10_{N}"))
    sste, fe = wsstpy.wsst(x, extend_signal=True)
    cmp_matrix(f"wsst extend N={N}", sste, getattr(W, f"sst_ext_{N}"))
    assert np.abs(fe - getattr(W, f"f_ext_{N}")).max() < 1e-15

x = S.tone_1024
for tag, kw in [("morse (3,60)", dict(wavelet="morse")),
                ("morse (3,30)", dict(wavelet="morse", wavelet_parameters=(3.0, 30.0))),
                ("morse (2,40)", dict(wavelet="morse", wavelet_parameters=(2.0, 40.0))),
                ("bump", dict(wavelet="bump")),
                ("fs=100", dict(fs=100.0))]:
    sst, f = wsstpy.wsst(x, **kw)
    key = {"morse (3,60)": "morse", "morse (3,30)": "morse30",
           "morse (2,40)": "morse_g2", "bump": "bump", "fs=100": "fs"}[tag]
    cmp_matrix(f"wsst {tag}", sst, getattr(W, f"sst_{key}"))
    cmp_vector(f"  f {tag}", f, getattr(W, f"f_{key}"))

print("running cwt / icwt comparisons ...")
cfs, fc = wsstpy.cwt(x)
cmp_matrix("cwt Morse default", cfs, C.cfs_morse)
cmp_vector("  f cwt Morse default", fc, C.f_morse)
cfs_lim, f_lim = wsstpy.cwt(x, frequency_limits=(0.02, 0.25))
cmp_matrix("cwt Morse limits [0.02 0.25]", cfs_lim, C.cfs_lim)
cmp_vector("  f cwt limits", f_lim, C.f_lim)

fb = wsstpy.CWTFilterBank(1024, 32, "amor", 1.0, None, "reflection")
cfs_amor = wsstpy.cwt(x, wavelet="amor", voices_per_octave=32)[0]
cmp_matrix("cwt amor 32 voices", cfs_amor, C.cfs_amor)

cmp_vector("icwt full", wsstpy.icwt(cfs), C.icwt_full)
cmp_vector("icwt band [0.03 0.08]", wsstpy.icwt(cfs, f=fc, freq_range=(0.03, 0.08),
                                                wavelet="morse"), C.icwt_band)
cmp_vector("icwt vpo=10", wsstpy.icwt(cfs, voices_per_octave=10), C.icwt_vpo)
cmp_vector("icwt amor", wsstpy.icwt(cfs_amor, wavelet="amor"), C.icwt_amor)
cmp_vector("icwt Morse limits [0.03 0.1]",
           wsstpy.icwt(cfs_lim, f=f_lim, freq_range=(0.03, 0.1), wavelet="morse"), C.icwt_lim)

print("running iwsst comparisons ...")
sst = W.sst_1024
cmp_vector("iwsst full", wsstpy.iwsst(sst), I.full)
cmp_vector("iwsst band [0.03 0.08]", wsstpy.iwsst(sst, W.f_1024, (0.03, 0.08)), I.band)
fr, ir = wsstpy.wsstridge(sst, penalty=5, num_ridges=2, f=W.f_1024)
cmp_matrix("iwsst by ridge", wsstpy.iwsst(sst, ridge_indices=ir), I.ridge)
cmp_vector("iwsst Morse", wsstpy.iwsst(W.sst_morse, wavelet="morse"), I.morse)

print("running ridge comparisons ...")
for pen in (0, 0.5, 1, 5, 10, 100):
    r = getattr(G, f"pen{int(pen * 10)}")
    fr_r, ir_r = wsstpy.wsstridge(sst, penalty=pen, f=W.f_1024)
    cmp_exact(f"wsstridge penalty={pen:g}", ir_r, r.iridge)
    cmp_vector(f"  fridge penalty={pen:g}", fr_r, r.fridge)
fr_r, ir_r = wsstpy.wsstridge(sst, penalty=5, num_ridges=2, f=W.f_1024)
cmp_exact("wsstridge 2 ridges", ir_r, G.r2.iridge)
fr_r, ir_r = wsstpy.wsstridge(sst, penalty=5, num_ridges=3, num_frequency_bins=6, f=W.f_1024)
cmp_exact("wsstridge 3 ridges, 6 bins", ir_r, G.r3.iridge)
xc = S.chirp_1024
sstc, fsc = wsstpy.wsst(xc)
fr_r, ir_r = wsstpy.wsstridge(sstc, penalty=2, f=fsc)
cmp_exact("wsstridge chirp penalty=2", ir_r, G.chirp.iridge)

print()
print(f"MATLAB {ref['R'].matlab_version} | Python {sys.version.split()[0]} | "
      f"numpy {np.__version__} | scipy {__import__('scipy').__version__}")
print()
print("| comparison | max abs diff | max rel diff | RMS diff |")
print("| --- | --- | --- | --- |")
for name, a, b, c in rows:
    if name.startswith("  "):
        continue
    print(f"| {name} | {a} | {b} | {c} |")

with open("/home/hilal/dev/wsstpy/probe/results5.json", "w") as fh:
    json.dump([{"name": n, "max_abs": a, "max_rel": b, "rms": c} for n, a, b, c in rows], fh, indent=1)
print()
print("rows:", len(rows))
