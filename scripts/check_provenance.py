"""Check that this repository contains no MATLAB source code.

Two passes.

1. Marker scan: strings that only appear in MATLAB toolbox source (codegen
   pragmas, package-private names, MathWorks copyright lines, mex scaffolding).
2. Shared-token scan: every window of this repository's text files is compared
   against every window of every .m file in a MATLAB installation. A hit means
   a run of consecutive tokens is identical. Windows with fewer than four
   identifier tokens are ignored, because punctuation runs and bare syntax
   ("end for k = 1 : numel") match everywhere and mean nothing.

This script and PROVENANCE.md are skipped, since both quote the marker strings
and the example snippets deliberately.

Usage:
    python scripts/check_provenance.py [--matlab /usr/local/MATLAB/R2025b]
"""
from __future__ import annotations

import argparse
import hashlib
import os
import re
import sys

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

TEXT_EXT = {".py", ".m", ".md", ".toml", ".txt", ".yml", ".yaml", ".ini", ".cfg"}

MARKERS = [
    "%#codegen",
    "coder.",
    "MathWorks",
    "bsxfun",
    "accumarray",
    "validateattributes",
    "inputParser",
    "message('Wavelet",
    'message("Wavelet',
    "quadgk",
    "fzero",
    "gammainc(",
    "underlyingType",
    "coder.target",
    "coder.internal",
    "@cwtfilterbank",
    "end %",
    "%%",
]

TOKEN_RE = re.compile(r"[A-Za-z_]\w*|\d+(?:\.\d*)?(?:[eE][-+]?\d+)?|[^\s\w]")
IDENT_RE = re.compile(r"^[a-z_][a-z0-9_]*$")


def ident_count(tokens):
    return sum(1 for t in tokens if IDENT_RE.match(t))


def tokenize(text: str):
    return TOKEN_RE.findall(text.lower())


def window_hashes(tokens, n):
    for i in range(len(tokens) - n + 1):
        key = " ".join(tokens[i : i + n])
        yield hashlib.blake2b(key.encode(), digest_size=8).digest(), i


# Files that quote the marker strings or the example snippets on purpose.
SKIP = {
    os.path.join("scripts", "check_provenance.py"),
    "PROVENANCE.md",
}


def skip(rel):
    return rel in SKIP


def repo_files():
    out = []
    for root, dirs, files in os.walk(REPO):
        dirs[:] = [d for d in dirs if d not in {".git", ".venv", "__pycache__", "figures"}]
        if os.sep + "probe" + os.sep + "data" in root:
            continue
        for name in files:
            if os.path.splitext(name)[1] in TEXT_EXT:
                full = os.path.join(root, name)
                if skip(os.path.relpath(full, REPO)):
                    continue
                out.append(full)
    return sorted(out)


def matlab_files(root):
    for base in ("toolbox/wavelet", "toolbox/shared/signalwavelet", "toolbox/signal/signal"):
        start = os.path.join(root, base)
        if not os.path.isdir(start):
            continue
        for b, _, files in os.walk(start):
            for name in files:
                if name.endswith(".m"):
                    path = os.path.join(b, name)
                    try:
                        if os.path.getsize(path) < 4_000_000:
                            yield path
                    except OSError:
                        pass


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--matlab", default="/usr/local/MATLAB/R2025b")
    ap.add_argument("--window", type=int, default=10)
    args = ap.parse_args()

    files = repo_files()
    print(f"repository text files: {len(files)}")

    # ---- pass 1: markers
    marker_hits = []
    for path in files:
        rel = os.path.relpath(path, REPO)
        with open(path, encoding="utf-8", errors="replace") as fh:
            for lineno, line in enumerate(fh, 1):
                for marker in MARKERS:
                    if marker in line:
                        marker_hits.append((rel, lineno, marker, line.strip()[:120]))
    print(f"\npass 1, marker scan: {len(marker_hits)} hit(s)")
    for rel, lineno, marker, line in marker_hits:
        print(f"  {rel}:{lineno}  [{marker}]  {line}")

    # ---- pass 2: shared token windows
    print(f"\npass 2, building {args.window}-token window index of MATLAB sources ...")
    index = {}
    n_mfiles = 0
    for path in matlab_files(args.matlab):
        n_mfiles += 1
        try:
            text = open(path, encoding="utf-8", errors="replace").read()
        except OSError:
            continue
        tokens = tokenize(text)
        for digest, pos in window_hashes(tokens, args.window):
            index.setdefault(digest, (path, pos))
    print(f"  {n_mfiles} .m files, {len(index)} distinct windows")

    print("\n  scanning repository files ...")
    hits = []
    for path in files:
        rel = os.path.relpath(path, REPO)
        text = open(path, encoding="utf-8", errors="replace").read()
        tokens = tokenize(text)
        for digest, pos in window_hashes(tokens, args.window):
            if ident_count(tokens[pos : pos + args.window]) < max(4, args.window // 2):
                continue
            if digest in index:
                mpath, mpos = index[digest]
                hits.append((rel, pos, mpath, tokens[pos : pos + args.window]))

    print(f"\npass 2 result: {len(hits)} window(s) of {args.window} tokens shared with MATLAB")
    for rel, pos, mpath, toks in hits[:40]:
        print(f"  {rel} @ token {pos}: {' '.join(toks)}")
        print(f"      in {os.path.relpath(mpath, args.matlab)}")

    # ---- pass 3: longest shared run per repo file, against the files that
    # were actually consulted while writing the port
    print("\npass 3, longest shared token run against the consulted files ...")
    shortlist = [
        "toolbox/wavelet/core/wavelet/wsst.m",
        "toolbox/wavelet/core/wavelet/iwsst.m",
        "toolbox/wavelet/core/wavelet/wsstridge.m",
        "toolbox/wavelet/core/wavelet/cwt.m",
        "toolbox/wavelet/core/wavelet/icwt.m",
        "toolbox/wavelet/core/wavelet/@cwtfilterbank/cwtfilterbank.m",
        "toolbox/wavelet/core/wavelet/@cwtfilterbank/wt.m",
        "toolbox/wavelet/core/wavelet/+wavelet/+internal/+cwt/wavbpfilters.m",
        "toolbox/wavelet/core/wavelet/+wavelet/+internal/+cwt/morsebpfilters.m",
        "toolbox/wavelet/core/wavelet/+wavelet/+internal/+cwt/getCWTScales.m",
        "toolbox/wavelet/core/wavelet/+wavelet/+internal/+cwt/cwtfreqlimits.m",
        "toolbox/wavelet/core/wavelet/+wavelet/+internal/+cwt/wavCFandSD.m",
        "toolbox/wavelet/core/wavelet/+wavelet/+internal/+cwt/findFreqIndices.m",
        "toolbox/wavelet/core/wavelet/+wavelet/+internal/+cwt/morsenormconstant.m",
        "toolbox/wavelet/core/wavelet/+wavelet/+internal/+cwt/morsepeakfreq.m",
        "toolbox/wavelet/core/wavelet/+wavelet/+internal/+cwt/admConstant.m",
        "toolbox/wavelet/core/wavelet/+wavelet/+internal/+cwt/morseproperties.m",
        "toolbox/shared/signalwavelet/signalwavelet/+signalwavelet/+internal/+tfridge/extractRidges.m",
    ]

    corpus = []
    for rel_m in shortlist:
        path_m = os.path.join(args.matlab, rel_m)
        if os.path.isfile(path_m):
            corpus.append((rel_m, tokenize(open(path_m, encoding="utf-8", errors="replace").read())))
    print(f"  {len(corpus)} MATLAB files in the shortlist")

    SEED = 6

    def longest_run(toks, mtoks):
        """Longest common token run, seeded on SEED-token windows."""
        index = {}
        for j in range(len(mtoks) - SEED + 1):
            index.setdefault(tuple(mtoks[j : j + SEED]), []).append(j)
        best = (0, 0, 0)
        for i in range(len(toks) - SEED + 1):
            starts = index.get(tuple(toks[i : i + SEED]))
            if not starts:
                continue
            for j in starts:
                a, b = i, j
                while a > 0 and b > 0 and toks[a - 1] == mtoks[b - 1]:
                    a -= 1
                    b -= 1
                length = 0
                while a + length < len(toks) and b + length < len(mtoks) \
                        and toks[a + length] == mtoks[b + length]:
                    length += 1
                if length > best[0] and ident_count(toks[a : a + length]) >= 5:
                    best = (length, a, b)
        return best

    worst = []
    for path in files:
        rel = os.path.relpath(path, REPO)
        toks = tokenize(open(path, encoding="utf-8", errors="replace").read())
        best = (0, "", "")
        for mrel, mtoks in corpus:
            size, a, b = longest_run(toks, mtoks)
            if size > best[0]:
                best = (size, " ".join(toks[a : a + size]), mrel)
        if best[0] >= 5:
            worst.append((rel, best[0], best[1], best[2]))

    worst.sort(key=lambda r: -r[1])
    print(f"  {len(worst)} file(s) with a run of 5 or more shared tokens")
    for rel, length, snippet, mpath in worst[:25]:
        print(f"  {rel}: {length} tokens  '{snippet}'")
        print(f"      {mpath}")

    print("\nsummary")
    print(f"  markers: {len(marker_hits)}")
    print(f"  shared {args.window}-token windows: {len(hits)}")
    print(f"  longest shared run: {max([w[1] for w in worst], default=0)} tokens")
    return 0


if __name__ == "__main__":
    sys.exit(main())