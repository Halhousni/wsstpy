"""Tests for the inverse transform, the CWT pair and ridge extraction."""
from __future__ import annotations

import os

import numpy as np
import pytest

import wsstpy

DATA = os.path.join(os.path.dirname(__file__), "data", "reference.npz")


@pytest.fixture(scope="module")
def ref():
    return np.load(DATA)


# --- cwt / icwt ----------------------------------------------------------
def test_cwt_morse(ref):
    cfs, f = wsstpy.cwt(ref["sig_tone"], 1.0)
    assert np.abs(cfs - ref["cwt_morse"]).max() < 1e-12
    assert np.abs(f - ref["f_cwt_morse"]).max() < 1e-15


def test_cwt_amor_filterbank(ref):
    cfs, _ = wsstpy.cwt(ref["sig_tone"], wavelet="amor", voices_per_octave=32)
    assert np.abs(cfs - ref["cwt_amor"]).max() < 1e-12


def test_cwt_frequency_limits(ref):
    cfs, f = wsstpy.cwt(ref["sig_tone"], frequency_limits=(0.02, 0.25))
    assert np.abs(cfs - ref["cwt_limits"]).max() < 1e-12
    assert np.abs(f - ref["f_cwt_limits"]).max() < 1e-15


def test_icwt_matches_matlab(ref):
    cfs = ref["cwt_morse"]
    f = ref["f_cwt_morse"]
    assert np.abs(wsstpy.icwt(cfs) - ref["icwt_full"]).max() < 1e-12
    assert np.abs(wsstpy.icwt(cfs, voices_per_octave=10) - ref["icwt_vpo10"]).max() < 1e-12
    out = wsstpy.icwt(cfs, f=f, freq_range=(0.03, 0.08), wavelet="morse")
    assert np.abs(out - ref["icwt_band"]).max() < 1e-12


def test_icwt_band_selects_rows(ref):
    cfs = ref["cwt_morse"]
    f = ref["f_cwt_morse"]
    full = wsstpy.icwt(cfs, wavelet="morse")
    band = wsstpy.icwt(cfs, f=f, freq_range=(0.03, 0.08), wavelet="morse")
    assert np.abs(band).max() < np.abs(full).max()
    assert np.corrcoef(band, full)[0, 1] > 0.5


def test_icwt_amor(ref):
    out = wsstpy.icwt(ref["cwt_amor"], wavelet="amor")
    # MATLAB evaluates the admissibility integral with integral()'s default
    # relative tolerance of 1e-6, so its constant carries a 7e-8 error.
    assert np.abs(out - ref["icwt_amor"]).max() < 1e-6


def test_icwt_roundtrip():
    n = 1024
    x = np.cos(2 * np.pi * 0.05 * np.arange(n)) + np.cos(2 * np.pi * 0.2 * np.arange(n))
    cfs, f = wsstpy.cwt(x, wavelet="morse")
    rec = wsstpy.icwt(cfs, f=f, wavelet="morse")
    assert np.abs(rec - x).max() < 0.35
    assert np.corrcoef(rec, x)[0, 1] > 0.99



# --- iwsst ---------------------------------------------------------------
def test_iwsst_full(ref):
    sst = ref["wsst_tone"]
    assert np.abs(wsstpy.iwsst(sst) - ref["iwsst_full"]).max() < 1e-6


def test_iwsst_band(ref):
    sst = ref["wsst_tone"]
    out = wsstpy.iwsst(sst, ref["f_tone"], (0.03, 0.08))
    assert np.abs(out - ref["iwsst_band"]).max() < 1e-6


def test_iwsst_morse(ref):
    out = wsstpy.iwsst(ref["wsst_morse"], wavelet="morse")
    assert np.abs(out - ref["iwsst_morse"]).max() < 1e-12


def test_iwsst_ridge_reconstruction(ref):
    sst = ref["wsst_tone"]
    ir = ref["iridge_two"]
    out = wsstpy.iwsst(sst, ridge_indices=ir)
    assert out.shape == ref["iwsst_ridge"].shape
    assert np.abs(out - ref["iwsst_ridge"]).max() < 1e-6


# --- ridges --------------------------------------------------------------
@pytest.mark.parametrize("pen", [0, 1, 5, 100])
def test_ridge_indices_match_matlab(ref, pen):
    sst = ref["wsst_tone"]
    fridge, iridge = wsstpy.wsstridge(sst, penalty=pen, f=ref["f_tone"])
    assert np.array_equal(iridge[:, 0], ref[f"iridge_pen{pen}"])
    assert np.abs(fridge[:, 0] - ref[f"fridge_pen{pen}"]).max() == 0


def test_two_ridges(ref):
    _, ir = wsstpy.wsstridge(ref["wsst_tone"], penalty=5, num_ridges=2, f=ref["f_tone"])
    assert np.array_equal(ir, ref["iridge_two"])


def test_three_ridges_with_bins(ref):
    _, ir = wsstpy.wsstridge(ref["wsst_tone"], penalty=5, num_ridges=3,
                             num_frequency_bins=6, f=ref["f_tone"])
    assert np.array_equal(ir, ref["iridge_three_bins6"])


def test_ridge_on_chirp(ref):
    _, ir = wsstpy.wsstridge(ref["wsst_chirp"], penalty=2, f=ref["f_chirp"])
    assert np.array_equal(ir[:, 0], ref["iridge_chirp"])


def test_ridge_on_noise(ref):
    _, ir = wsstpy.wsstridge(ref["wsst_noise"], penalty=3, f=ref["f_noise"])
    assert np.array_equal(ir[:, 0], ref["iridge_noise"])


def test_ridge_tracks_a_tone():
    n = 512
    sst, f = wsstpy.wsst(np.cos(2 * np.pi * 0.12 * np.arange(n)))
    fridge, iridge = wsstpy.wsstridge(sst, penalty=2.0, f=f)
    assert abs(np.median(fridge) - 0.12) < 0.01
    assert iridge.min() >= 1 and iridge.max() <= sst.shape[0]
