"""Tests that compare wsstpy against MATLAB's output.

Every reference array in ``data/reference.npz`` was produced by MATLAB R2025b
with the Wavelet Toolbox and saved from the probe scripts in ``probe/``.
"""
from __future__ import annotations

import os

import numpy as np
import pytest

import wsstpy

DATA = os.path.join(os.path.dirname(__file__), "data", "reference.npz")


@pytest.fixture(scope="module")
def ref():
    return np.load(DATA)


SIGNALS = ["tone", "chirp", "noise", "short", "impulse"]


@pytest.mark.parametrize("name", SIGNALS + ["odd"])
def test_wsst_matches_matlab(ref, name):
    x = ref[f"sig_{name}"]
    sst, f = wsstpy.wsst(x)
    expected = ref[f"wsst_{name}"]
    assert sst.shape == expected.shape
    assert np.abs(sst - expected).max() < 1e-12
    assert np.abs(f - ref[f"f_{name}"]).max() < 1e-15


def test_wsst_voices_per_octave(ref):
    sst, f = wsstpy.wsst(ref["sig_tone"], voices_per_octave=10)
    assert np.abs(sst - ref["wsst_vpo10"]).max() < 1e-12
    assert np.abs(f - ref["f_vpo10"]).max() < 1e-15


def test_wsst_extend_signal(ref):
    sst, f = wsstpy.wsst(ref["sig_tone"], extend_signal=True)
    assert np.abs(sst - ref["wsst_extend"]).max() < 1e-12
    assert np.abs(f - ref["f_extend"]).max() < 1e-15


def test_wsst_morse(ref):
    sst, f = wsstpy.wsst(ref["sig_tone"], wavelet="morse")
    assert np.abs(sst - ref["wsst_morse"]).max() < 1e-12
    assert np.abs(f - ref["f_morse"]).max() < 1e-15


def test_wsst_morse_parameters(ref):
    sst, _ = wsstpy.wsst(ref["sig_tone"], wavelet="morse", wavelet_parameters=(3, 30))
    assert np.abs(sst - ref["wsst_morse30"]).max() < 1e-12


def test_wsst_bump(ref):
    sst, f = wsstpy.wsst(ref["sig_tone"], wavelet="bump")
    assert np.abs(sst - ref["wsst_bump"]).max() < 1e-12
    assert np.abs(f - ref["f_bump"]).max() < 1e-15


def test_wsst_sampling_frequency(ref):
    sst, f = wsstpy.wsst(ref["sig_tone"], 100)
    assert np.abs(sst - ref["wsst_fs"]).max() < 1e-12
    assert np.abs(f - ref["f_fs"]).max() < 1e-12


def test_wsst_rows_are_ordered_by_frequency():
    x = np.cos(2 * np.pi * 0.1 * np.arange(512))
    sst, f = wsstpy.wsst(x)
    assert np.all(np.diff(f) > 0)
    # the tone must land on the row whose frequency is closest to 0.1
    row = np.argmax(np.abs(sst[:, 100]))
    assert abs(f[row] - 0.1) < 0.01


def test_frequencies_scale_with_sampling_frequency():
    x = np.cos(2 * np.pi * 0.1 * np.arange(256))
    _, f_norm = wsstpy.wsst(x)
    _, f_hz = wsstpy.wsst(x, 100)
    assert np.allclose(f_hz, 100 * f_norm)
