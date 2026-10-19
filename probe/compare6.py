"""Second validation pass: every waveform, plus ridges on chirp and noise."""
import json
import os
import sys

import numpy as np
import scipy.io as sio

sys.path.insert(0, "/home/hilal/dev/wsstpy/src")
import wsstpy

DATA = "/home/hilal/dev/wsstpy/probe/data6"
ref = sio.loadmat(os.path.join(DATA, "ref6.mat"), squeeze_me=True, struct_as_record=False)
S, W, G = ref["S"], ref["W"], ref["G"]

rows = []


def add(name, d, scale):
    rows.append((name, f"{d.max():.3e}", f"{d.max() / scale:.3e}" if scale else "0"))


Ns = [256, 512, 1000, 1024, 1025, 4096]
for N in Ns:
    for sig in ("chirp", "noise", "short", "impulse"):
        x = getattr(S, f"{sig}_{N}")
        sst, f = wsstpy.wsst(x)
        R = getattr(W, f"sst_{sig}_{N}")
        Rd = np.asarray(R)
        add(f"wsst {sig} N={N}", np.abs(sst - Rd), np.abs(Rd).max())
        F = np.asarray(getattr(W, f"f_{sig}_{N}")).ravel()
        assert np.abs(f - F).max() < 1e-15, (sig, N, np.abs(f - F).max())

print()
print("| case | max abs diff | max rel diff |")
print("| --- | --- | --- |")
for name, a, b in rows:
    print(f"| {name} | {a} | {b} |")

print()
print("ridge checks")
for tag, sst_key, pen in [("chirp2", "sst_chirp_1024", 2), ("chirp5", "sst_chirp_1024", 5),
                          ("chirp0", "sst_chirp_1024", 0), ("noise", "sst_noise_1024", 3)]:
    N = 1024
    x = getattr(S, sst_key.split("_")[1] + "_1024")
    sst, f = wsstpy.wsst(x)
    ref_sst = np.asarray(getattr(W, sst_key))
    fr, ir = wsstpy.wsstridge(sst, penalty=pen, f=f)
    ir_ref = np.asarray(getattr(G, tag).iridge).ravel().astype(int)
    same = int(np.sum(ir[:, 0] == ir_ref))
    print(f"  {tag}: sst matches? {np.abs(sst - ref_sst).max():.2e}; "
          f"ridge {same}/{ir.size} identical")

with open("/home/hilal/dev/wsstpy/probe/results6.json", "w") as fh:
    json.dump(rows, fh, indent=1)
