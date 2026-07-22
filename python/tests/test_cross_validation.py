# Copyright (c) 2026 Pedro Lourenco. All rights reserved.
# This code is released under the MIT License. See LICENSE file in the
# project root for full license information.
#
# This module is part of the Open Source System Identification Toolbox (SID).
# https://github.com/pdlourenco/sid

"""Cross-language validation tests using MATLAB-generated reference data.

These tests load JSON reference vectors produced by
``testdata/generate_reference.m`` (run in MATLAB CI) and verify that the
Python implementation produces numerically identical results.

The JSON files contain the actual input data (not RNG seeds), so results
are deterministic and RNG-independent.  Tests are skipped when the JSON
files are not available (e.g. before the first MATLAB CI run).
"""

from __future__ import annotations

import json
import pathlib

import numpy as np
import pytest

import sid
from sid._internal.cov import sid_cov
from sid._internal.hann_win import hann_win
from sid._internal.windowed_dft import windowed_dft

TESTDATA = pathlib.Path(__file__).resolve().parent.parent.parent / "testdata"


def _load(name: str) -> dict:
    """Load a JSON reference file, or skip the test if missing."""
    path = TESTDATA / name
    if not path.exists():
        pytest.skip(f"Reference file {name} not found (run MATLAB CI first)")
    with open(path) as f:
        return json.load(f)


def _to_array(data, key: str) -> np.ndarray:
    """Convert a JSON field to a numpy array."""
    return np.array(data[key], dtype=np.float64)


def _to_complex(data, base_key: str) -> np.ndarray:
    """Reconstruct complex array from _real and _imag JSON fields."""
    re = np.array(data[base_key + "_real"], dtype=np.float64)
    im = np.array(data[base_key + "_imag"], dtype=np.float64)
    return re + 1j * im


def _tol(ref: dict, field: str) -> dict:
    """Look up ``(rtol, atol)`` for an output field as an ``assert_allclose``
    kwargs dict, reading the JSON ``tolerance`` block as the single
    authoritative source.

    Mirrors ``testdata/validate_reference.m`` exactly: the default is
    ``rtol=1e-6``, ``atol=0`` when the field carries no ``<field>_rel`` /
    ``<field>_atol`` entry, so both language consumers apply identical
    thresholds. Tests must not hardcode their own tolerances — tighten or
    loosen a comparison by editing the vector's tolerance block in
    ``generate_reference.m`` instead, where it binds every language.
    """
    tolerance = ref.get("tolerance", {})
    return {
        "rtol": tolerance.get(field + "_rel", 1e-6),
        "atol": tolerance.get(field + "_atol", 0.0),
    }


class TestCrossValidationSISO:
    """SISO Blackman-Tukey: reference_siso_bt.json."""

    def test_siso_bt_response(self):
        ref = _load("reference_siso_bt.json")
        y = _to_array(ref["input"], "y")
        u = _to_array(ref["input"], "u")
        ws = ref["params"]["WindowSize"]
        ts = ref["params"]["SampleTime"]

        result = sid.freq_bt(y, u, window_size=ws, sample_time=ts)

        expected_resp = _to_complex(ref["output"], "Response")
        np.testing.assert_allclose(
            result.response,
            expected_resp.ravel(),
            **_tol(ref, "Response"),
            err_msg="SISO BT response mismatch vs MATLAB reference",
        )

    def test_siso_bt_frequency(self):
        ref = _load("reference_siso_bt.json")
        y = _to_array(ref["input"], "y")
        u = _to_array(ref["input"], "u")
        ws = ref["params"]["WindowSize"]

        result = sid.freq_bt(y, u, window_size=ws)

        expected_freq = _to_array(ref["output"], "Frequency")
        np.testing.assert_allclose(
            result.frequency,
            expected_freq.ravel(),
            **_tol(ref, "Frequency"),
            err_msg="SISO BT frequency grid mismatch",
        )

    def test_siso_bt_noise_spectrum(self):
        ref = _load("reference_siso_bt.json")
        y = _to_array(ref["input"], "y")
        u = _to_array(ref["input"], "u")
        ws = ref["params"]["WindowSize"]

        result = sid.freq_bt(y, u, window_size=ws)

        expected_ns = _to_array(ref["output"], "NoiseSpectrum")
        np.testing.assert_allclose(
            result.noise_spectrum,
            expected_ns.ravel(),
            **_tol(ref, "NoiseSpectrum"),
            err_msg="SISO BT noise spectrum mismatch vs MATLAB reference",
        )

    def test_siso_bt_coherence(self):
        ref = _load("reference_siso_bt.json")
        y = _to_array(ref["input"], "y")
        u = _to_array(ref["input"], "u")
        ws = ref["params"]["WindowSize"]

        result = sid.freq_bt(y, u, window_size=ws)

        expected_coh = _to_array(ref["output"], "Coherence")
        np.testing.assert_allclose(
            result.coherence,
            expected_coh.ravel(),
            **_tol(ref, "Coherence"),
            err_msg="SISO BT coherence mismatch vs MATLAB reference",
        )


class TestCrossValidationSISOLargeM:
    """SISO Blackman-Tukey at M = 200: reference_siso_bt_large_M.json.

    Regression vector for SPEC.md S2.5.1.  Pins Python/MATLAB agreement
    in the previously broken FFT-fast-path window-size regime
    (``M >= nf = 128``), where the legacy hardcoded ``L = 2*nf = 256``
    caused silent positive/negative lag overlap in both languages.  Any
    future change that re-introduces the hardcode will flip every
    assertion in this class.
    """

    def test_siso_bt_large_M_response(self):
        ref = _load("reference_siso_bt_large_M.json")
        y = _to_array(ref["input"], "y")
        u = _to_array(ref["input"], "u")
        ws = ref["params"]["WindowSize"]
        ts = ref["params"]["SampleTime"]

        result = sid.freq_bt(y, u, window_size=ws, sample_time=ts)

        expected_resp = _to_complex(ref["output"], "Response")
        np.testing.assert_allclose(
            result.response,
            expected_resp.ravel(),
            **_tol(ref, "Response"),
            err_msg="SISO BT (large M) response mismatch vs MATLAB reference",
        )

    def test_siso_bt_large_M_frequency(self):
        ref = _load("reference_siso_bt_large_M.json")
        y = _to_array(ref["input"], "y")
        u = _to_array(ref["input"], "u")
        ws = ref["params"]["WindowSize"]

        result = sid.freq_bt(y, u, window_size=ws)

        expected_freq = _to_array(ref["output"], "Frequency")
        np.testing.assert_allclose(
            result.frequency,
            expected_freq.ravel(),
            **_tol(ref, "Frequency"),
            err_msg="SISO BT (large M) frequency grid mismatch",
        )

    def test_siso_bt_large_M_noise_spectrum(self):
        ref = _load("reference_siso_bt_large_M.json")
        y = _to_array(ref["input"], "y")
        u = _to_array(ref["input"], "u")
        ws = ref["params"]["WindowSize"]

        result = sid.freq_bt(y, u, window_size=ws)

        expected_ns = _to_array(ref["output"], "NoiseSpectrum")
        np.testing.assert_allclose(
            result.noise_spectrum,
            expected_ns.ravel(),
            **_tol(ref, "NoiseSpectrum"),
            err_msg="SISO BT (large M) noise spectrum mismatch vs MATLAB reference",
        )

    def test_siso_bt_large_M_coherence(self):
        ref = _load("reference_siso_bt_large_M.json")
        y = _to_array(ref["input"], "y")
        u = _to_array(ref["input"], "u")
        ws = ref["params"]["WindowSize"]

        result = sid.freq_bt(y, u, window_size=ws)

        expected_coh = _to_array(ref["output"], "Coherence")
        np.testing.assert_allclose(
            result.coherence,
            expected_coh.ravel(),
            **_tol(ref, "Coherence"),
            err_msg="SISO BT (large M) coherence mismatch vs MATLAB reference",
        )


class TestCrossValidationMIMO:
    """MIMO Blackman-Tukey: reference_mimo_bt.json."""

    def test_mimo_bt_response(self):
        ref = _load("reference_mimo_bt.json")
        y = _to_array(ref["input"], "y")
        u = _to_array(ref["input"], "u")
        ws = ref["params"]["WindowSize"]

        result = sid.freq_bt(y, u, window_size=ws)

        expected_resp = _to_complex(ref["output"], "Response")
        # MATLAB stores as (nf x ny x nu); JSON flattens row-major
        # Reshape to match
        nf = result.response.shape[0]
        ny = result.response.shape[1]
        nu = result.response.shape[2]
        expected_3d = expected_resp.reshape(nf, ny, nu)

        np.testing.assert_allclose(
            result.response,
            expected_3d,
            **_tol(ref, "Response"),
            err_msg="MIMO BT response mismatch vs MATLAB reference",
        )

    def test_mimo_bt_noise_spectrum(self):
        ref = _load("reference_mimo_bt.json")
        y = _to_array(ref["input"], "y")
        u = _to_array(ref["input"], "u")
        ws = ref["params"]["WindowSize"]

        result = sid.freq_bt(y, u, window_size=ws)

        expected_ns = _to_array(ref["output"], "NoiseSpectrum")
        nf = result.noise_spectrum.shape[0]
        ny = result.noise_spectrum.shape[1]
        expected_3d = expected_ns.reshape(nf, ny, ny)

        np.testing.assert_allclose(
            result.noise_spectrum,
            expected_3d,
            **_tol(ref, "NoiseSpectrum"),
            err_msg="MIMO BT noise spectrum mismatch vs MATLAB reference",
        )


class TestCrossValidationTimeSeries:
    """Time-series Blackman-Tukey: reference_timeseries_bt.json."""

    def test_timeseries_bt_noise_spectrum(self):
        ref = _load("reference_timeseries_bt.json")
        y = _to_array(ref["input"], "y")
        ws = ref["params"]["WindowSize"]

        result = sid.freq_bt(y, None, window_size=ws)

        expected_ns = _to_array(ref["output"], "NoiseSpectrum")
        np.testing.assert_allclose(
            result.noise_spectrum,
            expected_ns.ravel(),
            **_tol(ref, "NoiseSpectrum"),
            err_msg="Time-series BT noise spectrum mismatch vs MATLAB reference",
        )

    def test_timeseries_bt_noise_spectrum_std(self):
        ref = _load("reference_timeseries_bt.json")
        y = _to_array(ref["input"], "y")
        ws = ref["params"]["WindowSize"]

        result = sid.freq_bt(y, None, window_size=ws)

        expected_std = _to_array(ref["output"], "NoiseSpectrumStd")
        np.testing.assert_allclose(
            result.noise_spectrum_std,
            expected_std.ravel(),
            **_tol(ref, "NoiseSpectrumStd"),
            err_msg="Time-series BT noise spectrum std mismatch vs MATLAB reference",
        )


class TestCrossValidationInternals:
    """Internal helper functions: reference_internals.json."""

    def test_covariance_auto(self):
        ref = _load("reference_internals.json")
        x = _to_array(ref["input"], "x").ravel()
        M = int(ref["input"]["M"])

        x2d = x[:, np.newaxis]
        R = sid_cov(x2d, x2d, M)

        expected = _to_array(ref["output"], "R_xx").ravel()
        np.testing.assert_allclose(
            R,
            expected,
            **_tol(ref, "R_xx"),
            err_msg="Auto-covariance mismatch vs MATLAB reference",
        )

    def test_covariance_cross(self):
        ref = _load("reference_internals.json")
        x = _to_array(ref["input"], "x").ravel()
        z = _to_array(ref["input"], "z").ravel()
        M = int(ref["input"]["M"])

        R = sid_cov(x[:, np.newaxis], z[:, np.newaxis], M)

        expected = _to_array(ref["output"], "R_xz").ravel()
        np.testing.assert_allclose(
            R,
            expected,
            **_tol(ref, "R_xz"),
            err_msg="Cross-covariance mismatch vs MATLAB reference",
        )

    def test_hann_window(self):
        ref = _load("reference_internals.json")
        M = int(ref["input"]["M"])

        W = hann_win(M)

        expected = _to_array(ref["output"], "W").ravel()
        np.testing.assert_allclose(
            W,
            expected,
            **_tol(ref, "W"),
            err_msg="Hann window mismatch vs MATLAB reference",
        )

    def test_windowed_dft_auto(self):
        ref = _load("reference_internals.json")
        x = _to_array(ref["input"], "x").ravel()
        M = int(ref["input"]["M"])

        x2d = x[:, np.newaxis]
        R = sid_cov(x2d, x2d, M)
        W = hann_win(M)
        freqs = np.arange(1, 129) * np.pi / 128

        Phi = windowed_dft(R, W, freqs, True, R)

        expected = _to_complex(ref["output"], "Phi_xx")
        np.testing.assert_allclose(
            Phi,
            expected.ravel(),
            rtol=ref["tolerance"]["Phi_xx_real_rel"],
            err_msg="Auto-spectrum windowed DFT mismatch vs MATLAB reference",
        )

    def test_windowed_dft_cross(self):
        ref = _load("reference_internals.json")
        x = _to_array(ref["input"], "x").ravel()
        z = _to_array(ref["input"], "z").ravel()
        M = int(ref["input"]["M"])

        x2d = x[:, np.newaxis]
        z2d = z[:, np.newaxis]
        R_xz = sid_cov(x2d, z2d, M)
        R_zx = sid_cov(z2d, x2d, M)
        W = hann_win(M)
        freqs = np.arange(1, 129) * np.pi / 128

        Phi = windowed_dft(R_xz, W, freqs, True, R_zx)

        expected = _to_complex(ref["output"], "Phi_xz")
        np.testing.assert_allclose(
            Phi,
            expected.ravel(),
            rtol=ref["tolerance"]["Phi_xz_real_rel"],
            err_msg="Cross-spectrum windowed DFT mismatch vs MATLAB reference",
        )


class TestCrossValidationETFE:
    """SISO ETFE: reference_siso_etfe.json."""

    def test_etfe_response(self):
        ref = _load("reference_siso_etfe.json")
        y = _to_array(ref["input"], "y")
        u = _to_array(ref["input"], "u")
        smoothing = ref["params"]["Smoothing"]
        ts = ref["params"]["SampleTime"]

        result = sid.freq_etfe(y, u, smoothing=smoothing, sample_time=ts)

        expected_resp = _to_complex(ref["output"], "Response")
        np.testing.assert_allclose(
            result.response,
            expected_resp.ravel(),
            **_tol(ref, "Response"),
            err_msg="ETFE response mismatch vs MATLAB reference",
        )

    def test_etfe_noise_spectrum(self):
        ref = _load("reference_siso_etfe.json")
        y = _to_array(ref["input"], "y")
        u = _to_array(ref["input"], "u")
        smoothing = ref["params"]["Smoothing"]

        result = sid.freq_etfe(y, u, smoothing=smoothing)

        expected_ns = _to_array(ref["output"], "NoiseSpectrum")
        np.testing.assert_allclose(
            result.noise_spectrum,
            expected_ns.ravel(),
            **_tol(ref, "NoiseSpectrum"),
            err_msg="ETFE noise spectrum mismatch vs MATLAB reference",
        )


class TestCrossValidationBTFDR:
    """SISO BTFDR: reference_siso_btfdr.json."""

    def test_btfdr_response(self):
        ref = _load("reference_siso_btfdr.json")
        y = _to_array(ref["input"], "y")
        u = _to_array(ref["input"], "u")
        ts = ref["params"]["SampleTime"]

        result = sid.freq_btfdr(y, u, sample_time=ts)

        expected_resp = _to_complex(ref["output"], "Response")
        np.testing.assert_allclose(
            result.response,
            expected_resp.ravel(),
            **_tol(ref, "Response"),
            err_msg="BTFDR response mismatch vs MATLAB reference",
        )

    def test_btfdr_noise_spectrum(self):
        ref = _load("reference_siso_btfdr.json")
        y = _to_array(ref["input"], "y")
        u = _to_array(ref["input"], "u")

        result = sid.freq_btfdr(y, u)

        expected_ns = _to_array(ref["output"], "NoiseSpectrum")
        np.testing.assert_allclose(
            result.noise_spectrum,
            expected_ns.ravel(),
            **_tol(ref, "NoiseSpectrum"),
            err_msg="BTFDR noise spectrum mismatch vs MATLAB reference",
        )


class TestCrossValidationSpectrogram:
    """Spectrogram: reference_spectrogram.json."""

    def test_spectrogram_time(self):
        ref = _load("reference_spectrogram.json")
        x = _to_array(ref["input"], "x")
        wl = ref["params"]["WindowLength"]
        ov = ref["params"]["Overlap"]
        ts = ref["params"]["SampleTime"]

        from sid.spectrogram import spectrogram

        result = spectrogram(x, window_length=wl, overlap=ov, sample_time=ts)

        expected_time = _to_array(ref["output"], "Time")
        np.testing.assert_allclose(
            result.time,
            expected_time.ravel(),
            **_tol(ref, "Time"),
            err_msg="Spectrogram time vector mismatch vs MATLAB reference",
        )

    def test_spectrogram_power(self):
        ref = _load("reference_spectrogram.json")
        x = _to_array(ref["input"], "x")
        wl = ref["params"]["WindowLength"]
        ov = ref["params"]["Overlap"]
        ts = ref["params"]["SampleTime"]

        from sid.spectrogram import spectrogram

        result = spectrogram(x, window_length=wl, overlap=ov, sample_time=ts)

        expected_power = _to_array(ref["output"], "Power")
        np.testing.assert_allclose(
            result.power.ravel(),
            expected_power.ravel(),
            **_tol(ref, "Power"),
            err_msg="Spectrogram power mismatch vs MATLAB reference",
        )

    def test_spectrogram_frequency(self):
        ref = _load("reference_spectrogram.json")
        x = _to_array(ref["input"], "x")
        wl = ref["params"]["WindowLength"]
        ov = ref["params"]["Overlap"]
        ts = ref["params"]["SampleTime"]

        from sid.spectrogram import spectrogram

        result = spectrogram(x, window_length=wl, overlap=ov, sample_time=ts)

        expected_freq = _to_array(ref["output"], "Frequency")
        np.testing.assert_allclose(
            result.frequency,
            expected_freq.ravel(),
            **_tol(ref, "Frequency"),
            err_msg="Spectrogram frequency mismatch vs MATLAB reference",
        )


class TestCrossValidationFreqMap:
    """FreqMap BT: reference_freqmap_bt.json."""

    def test_freqmap_response(self):
        ref = _load("reference_freqmap_bt.json")
        y = _to_array(ref["input"], "y")
        u = _to_array(ref["input"], "u")
        seg_len = ref["params"]["SegmentLength"]
        ov = ref["params"]["Overlap"]
        ws = ref["params"]["WindowSize"]

        from sid.freq_map import freq_map

        result = freq_map(y, u, segment_length=seg_len, overlap=ov, window_size=ws, algorithm="bt")

        expected_resp = _to_complex(ref["output"], "Response")
        np.testing.assert_allclose(
            result.response.ravel(),
            expected_resp.ravel(),
            **_tol(ref, "Response"),
            err_msg="FreqMap BT response mismatch vs MATLAB reference",
        )

    def test_freqmap_noise_spectrum(self):
        ref = _load("reference_freqmap_bt.json")
        y = _to_array(ref["input"], "y")
        u = _to_array(ref["input"], "u")
        seg_len = ref["params"]["SegmentLength"]
        ov = ref["params"]["Overlap"]
        ws = ref["params"]["WindowSize"]

        from sid.freq_map import freq_map

        result = freq_map(y, u, segment_length=seg_len, overlap=ov, window_size=ws, algorithm="bt")

        expected_ns = _to_array(ref["output"], "NoiseSpectrum")
        np.testing.assert_allclose(
            result.noise_spectrum.ravel(),
            expected_ns.ravel(),
            **_tol(ref, "NoiseSpectrum"),
            err_msg="FreqMap BT noise spectrum mismatch vs MATLAB reference",
        )

    def test_freqmap_coherence(self):
        ref = _load("reference_freqmap_bt.json")
        y = _to_array(ref["input"], "y")
        u = _to_array(ref["input"], "u")
        seg_len = ref["params"]["SegmentLength"]
        ov = ref["params"]["Overlap"]
        ws = ref["params"]["WindowSize"]

        from sid.freq_map import freq_map

        result = freq_map(y, u, segment_length=seg_len, overlap=ov, window_size=ws, algorithm="bt")

        expected_coh = _to_array(ref["output"], "Coherence")
        np.testing.assert_allclose(
            result.coherence.ravel(),
            expected_coh.ravel(),
            **_tol(ref, "Coherence"),
            err_msg="FreqMap BT coherence mismatch vs MATLAB reference",
        )

    def test_freqmap_time(self):
        ref = _load("reference_freqmap_bt.json")
        y = _to_array(ref["input"], "y")
        u = _to_array(ref["input"], "u")
        seg_len = ref["params"]["SegmentLength"]
        ov = ref["params"]["Overlap"]
        ws = ref["params"]["WindowSize"]

        from sid.freq_map import freq_map

        result = freq_map(y, u, segment_length=seg_len, overlap=ov, window_size=ws, algorithm="bt")

        expected_time = _to_array(ref["output"], "Time")
        np.testing.assert_allclose(
            result.time,
            expected_time.ravel(),
            **_tol(ref, "Time"),
            err_msg="FreqMap BT time vector mismatch vs MATLAB reference",
        )


class TestCrossValidationLTVCosmic:
    """LTV COSMIC: reference_ltv_cosmic.json."""

    def test_cosmic_a(self):
        ref = _load("reference_ltv_cosmic.json")
        X = _to_array(ref["input"], "X")
        U = _to_array(ref["input"], "U")
        lam = ref["params"]["Lambda"]

        from sid.ltv_disc import ltv_disc

        result = ltv_disc(X, U, lambda_=lam)

        expected_A = _to_array(ref["output"], "A")
        np.testing.assert_allclose(
            result.a.ravel(),
            expected_A.ravel(),
            **_tol(ref, "A"),
            err_msg="COSMIC A matrices mismatch vs MATLAB reference",
        )

    def test_cosmic_b(self):
        ref = _load("reference_ltv_cosmic.json")
        X = _to_array(ref["input"], "X")
        U = _to_array(ref["input"], "U")
        lam = ref["params"]["Lambda"]

        from sid.ltv_disc import ltv_disc

        result = ltv_disc(X, U, lambda_=lam)

        expected_B = _to_array(ref["output"], "B")
        np.testing.assert_allclose(
            result.b.ravel(),
            expected_B.ravel(),
            **_tol(ref, "B"),
            err_msg="COSMIC B matrices mismatch vs MATLAB reference",
        )

    def test_cosmic_cost(self):
        ref = _load("reference_ltv_cosmic.json")
        X = _to_array(ref["input"], "X")
        U = _to_array(ref["input"], "U")
        lam = ref["params"]["Lambda"]

        from sid.ltv_disc import ltv_disc

        result = ltv_disc(X, U, lambda_=lam)

        expected_cost = _to_array(ref["output"], "Cost")
        np.testing.assert_allclose(
            result.cost,
            expected_cost.ravel(),
            **_tol(ref, "Cost"),
            err_msg="COSMIC cost mismatch vs MATLAB reference",
        )


class TestCrossValidationTestMSD:
    """MSD plant: reference_test_msd.json.

    Verifies that the Python ``util_msd`` example helper produces the
    same discrete state-space matrices as the MATLAB port for the
    reference input parameters.
    """

    def test_msd_matrices(self):
        ref = _load("reference_test_msd.json")
        # util_msd lives in python/examples/util_msd.py; conftest.py
        # prepends that directory to sys.path so this import works
        # without sys.path manipulation here.
        from util_msd import util_msd

        m = _to_array(ref["input"], "m").ravel()
        k_spring = _to_array(ref["input"], "k_spring").ravel()
        c_damp = _to_array(ref["input"], "c_damp").ravel()
        F = _to_array(ref["input"], "F")
        Ts = float(ref["input"]["Ts"])

        Ad, Bd = util_msd(m, k_spring, c_damp, F, Ts)

        expected_Ad = _to_array(ref["output"], "Ad")
        expected_Bd = _to_array(ref["output"], "Bd")
        np.testing.assert_allclose(
            Ad, expected_Ad, err_msg="MSD Ad mismatch vs MATLAB", **_tol(ref, "Ad")
        )
        np.testing.assert_allclose(
            Bd, expected_Bd, err_msg="MSD Bd mismatch vs MATLAB", **_tol(ref, "Bd")
        )


class TestCrossValidationCosmicInternals:
    """COSMIC internals: reference_cosmic_internals.json."""

    def test_cosmic_internals(self):
        ref = _load("reference_cosmic_internals.json")
        X = _to_array(ref["input"], "X")
        U = _to_array(ref["input"], "U")
        if U.ndim == 1:
            U = U[:, np.newaxis]
        lam_val = float(ref["input"]["lambda"])

        from sid._internal.ltv_build_data_matrices import build_data_matrices
        from sid._internal.ltv_build_block_terms import build_block_terms
        from sid._internal.ltv_cosmic_solve import cosmic_solve
        from sid._internal.ltv_evaluate_cost import evaluate_cost
        from sid._internal.ltv_uncertainty_backward_pass import (
            uncertainty_backward_pass,
        )

        N = X.shape[0] - 1
        p = X.shape[1]
        q = U.shape[1]
        d = p + q
        L = 1
        X3 = X[:, :, np.newaxis]
        U3 = U[:, :, np.newaxis]
        lam = lam_val * np.ones(N - 1)

        D, Xl = build_data_matrices(X3, U3, N, p, q, L)
        S, T = build_block_terms(D, Xl, lam, N, p, q)
        C, _ = cosmic_solve(S, T, lam, N, p, q)

        # Extract A, B
        A_est = C[:p, :, :].transpose(1, 0, 2)
        B_est = C[p:, :, :].transpose(1, 0, 2)
        cost, fid, reg, _ = evaluate_cost(A_est, B_est, D, Xl, lam, N, p, q)

        # Bayesian posterior blocks. Pass S directly, matching production
        # (ltv_disc.py:314): S already carries the 1/sqrt(N) normalization.
        # The stored reference is generated the same way (issue #145a); the
        # 3-D arrays (d,d,N) flatten in matching C-order on both sides.
        P = uncertainty_backward_pass(S, lam, N, d)

        # Intermediate matrices: data (D, Xl), block terms (S, T), and the
        # COSMIC solution (C) — each an exact-arithmetic cross-check of the
        # whole pipeline, not just the scalar cost triple.
        for name, actual in [
            ("D", D),
            ("Xl", Xl),
            ("S", S),
            ("T", T),
            ("C", C),
            ("P", P),
        ]:
            expected = _to_array(ref["output"], name)
            np.testing.assert_allclose(
                actual.ravel(),
                expected.ravel(),
                err_msg=f"COSMIC {name} mismatch vs MATLAB reference",
                **_tol(ref, name),
            )

        expected_cost = float(ref["output"]["cost"])
        expected_fid = float(ref["output"]["fidelity"])
        expected_reg = float(ref["output"]["regularization"])

        np.testing.assert_allclose(
            cost, expected_cost, err_msg="COSMIC cost mismatch", **_tol(ref, "cost")
        )
        np.testing.assert_allclose(
            fid, expected_fid, err_msg="COSMIC fidelity mismatch", **_tol(ref, "fidelity")
        )
        np.testing.assert_allclose(
            reg,
            expected_reg,
            err_msg="COSMIC regularization mismatch",
            **_tol(ref, "regularization"),
        )


class TestCrossValidationLTVFrozen:
    """LTV frozen TF: reference_ltv_frozen.json."""

    def test_frozen_response(self):
        ref = _load("reference_ltv_frozen.json")
        X = _to_array(ref["input"], "X")
        U = _to_array(ref["input"], "U")
        if U.ndim == 1:
            U = U[:, np.newaxis]
        lam = ref["params"]["Lambda"]
        # MATLAB 1-based time steps -> Python 0-based
        time_steps_matlab = np.array(ref["params"]["frozen_TimeSteps"], dtype=int)
        time_steps = time_steps_matlab - 1

        from sid.ltv_disc import ltv_disc
        from sid.ltv_disc_frozen import ltv_disc_frozen

        ltv = ltv_disc(X, U, lambda_=lam)
        frz = ltv_disc_frozen(ltv, time_steps=time_steps)

        expected_resp = _to_complex(ref["output"], "Response")
        np.testing.assert_allclose(
            frz.response.ravel(),
            expected_resp.ravel(),
            **_tol(ref, "Response"),
            err_msg="Frozen TF response mismatch vs MATLAB",
        )


class TestCrossValidationModelOrder:
    """Model order estimation: reference_model_order.json."""

    def test_model_order(self):
        ref = _load("reference_model_order.json")
        y = _to_array(ref["input"], "y")
        u = _to_array(ref["input"], "u")
        if u.ndim == 1:
            u = u[:, np.newaxis]
        horizon = int(ref["params"]["Horizon"])
        bt_ws = int(ref["params"]["bt_WindowSize"])

        from sid.freq_bt import freq_bt
        from sid.model_order import model_order

        r = freq_bt(y, u, window_size=bt_ws)
        n_est, sv = model_order(r, horizon=horizon)

        expected_sv = _to_array(ref["output"], "SingularValues").ravel()
        np.testing.assert_allclose(
            sv["singular_values"][: len(expected_sv)],
            expected_sv,
            **_tol(ref, "SingularValues"),
            err_msg="Model order SVs mismatch vs MATLAB",
        )


class TestCrossValidationLTVStateEst:
    """LTV state estimation: reference_ltv_state_est.json."""

    def test_state_est(self):
        ref = _load("reference_ltv_state_est.json")
        Y = _to_array(ref["input"], "Y")
        U = _to_array(ref["input"], "U")
        if U.ndim == 1:
            U = U[:, np.newaxis]
        A = _to_array(ref["input"], "A")
        B = _to_array(ref["input"], "B")
        H = _to_array(ref["input"], "H")

        from sid.ltv_state_est import ltv_state_est

        X_hat = ltv_state_est(Y, U, A, B, H)
        if isinstance(X_hat, np.ndarray) and X_hat.ndim == 3 and X_hat.shape[2] == 1:
            X_hat = X_hat[:, :, 0]

        expected = _to_array(ref["output"], "X_hat")
        np.testing.assert_allclose(
            X_hat,
            expected,
            **_tol(ref, "X_hat"),
            err_msg="LTV state estimation mismatch vs MATLAB reference",
        )


class TestCrossValidationLTIFreqIO:
    """LTI freq IO: reference_lti_freq_io.json."""

    def test_lti_freq_io(self):
        ref = _load("reference_lti_freq_io.json")
        Y = _to_array(ref["input"], "Y")
        U = _to_array(ref["input"], "U")
        if U.ndim == 1:
            U = U[:, np.newaxis]
        H = _to_array(ref["input"], "H")

        from sid.lti_freq_io import lti_freq_io

        A0, B0 = lti_freq_io(Y, U, H)

        expected_A0 = _to_array(ref["output"], "A0")
        expected_B0 = _to_array(ref["output"], "B0")
        if expected_B0.ndim == 1:
            expected_B0 = expected_B0[:, np.newaxis]

        np.testing.assert_allclose(
            A0,
            expected_A0,
            **_tol(ref, "A0"),
            err_msg="LTI A0 mismatch vs MATLAB reference",
        )
        np.testing.assert_allclose(
            B0,
            expected_B0,
            **_tol(ref, "B0"),
            err_msg="LTI B0 mismatch vs MATLAB reference",
        )


class TestCrossValidationResidual:
    """Residual analysis: reference_residual.json."""

    def test_residual(self):
        ref = _load("reference_residual.json")
        y = _to_array(ref["input"], "y")
        u = _to_array(ref["input"], "u")
        if u.ndim == 1:
            u = u[:, np.newaxis]
        bt_ws = int(ref["params"]["bt_WindowSize"])
        max_lag = int(ref["params"]["MaxLag"])

        from sid.freq_bt import freq_bt
        from sid.residual import residual

        model = freq_bt(y, u, window_size=bt_ws)
        result = residual(model, y, u, max_lag=max_lag)

        expected_resid = _to_array(ref["output"], "Residual")
        expected_ac = _to_array(ref["output"], "AutoCorr")
        expected_cc = _to_array(ref["output"], "CrossCorr")

        np.testing.assert_allclose(
            result.residual.ravel(),
            expected_resid.ravel(),
            **_tol(ref, "Residual"),
            err_msg="Residual mismatch vs MATLAB",
        )
        np.testing.assert_allclose(
            result.auto_corr.ravel(),
            expected_ac.ravel(),
            **_tol(ref, "AutoCorr"),
            err_msg="AutoCorr mismatch vs MATLAB",
        )
        np.testing.assert_allclose(
            result.cross_corr.ravel(),
            expected_cc.ravel(),
            **_tol(ref, "CrossCorr"),
            err_msg="CrossCorr mismatch vs MATLAB",
        )


class TestCrossValidationCompare:
    """Model comparison: reference_compare.json."""

    def test_compare(self):
        ref = _load("reference_compare.json")
        X = _to_array(ref["input"], "X")
        U = _to_array(ref["input"], "U")
        if U.ndim == 1:
            U = U[:, np.newaxis]
        lam = ref["params"]["Lambda"]

        from sid.ltv_disc import ltv_disc
        from sid.compare import compare

        model = ltv_disc(X, U, lambda_=lam)
        result = compare(model, X, U)

        expected_pred = _to_array(ref["output"], "Predicted")
        expected_fit = _to_array(ref["output"], "Fit")

        np.testing.assert_allclose(
            result.predicted.ravel(),
            expected_pred.ravel(),
            **_tol(ref, "Predicted"),
            err_msg="Compare predicted mismatch vs MATLAB",
        )
        np.testing.assert_allclose(
            result.fit.ravel(),
            expected_fit.ravel(),
            **_tol(ref, "Fit"),
            err_msg="Compare fit mismatch vs MATLAB",
        )


# ---------------------------------------------------------------------------
# Previously-orphaned reference vectors (issue #145b).
#
# These four JSONs are generated by generate_reference.m and validated by
# validate_reference.m under Octave, but until now no Python test consumed
# them, leaving the Python ports of detrend, freq_domain_sim, ltv_disc_io,
# and the frequency-domain uncertainty helper with no cross-language check.
# Each test mirrors the corresponding branch of validate_reference.m's
# callSidFunction() so both consumers exercise the same call.
# ---------------------------------------------------------------------------


class TestCrossValidationDetrend:
    """Polynomial detrending: reference_detrend.json."""

    def test_detrend(self):
        ref = _load("reference_detrend.json")
        x = _to_array(ref["input"], "x")
        order = int(ref["params"]["Order"])

        from sid.detrend import detrend

        x_det, trend = detrend(x, order=order)

        np.testing.assert_allclose(
            x_det.ravel(),
            _to_array(ref["output"], "x_detrended").ravel(),
            err_msg="Detrend x_detrended mismatch vs MATLAB",
            **_tol(ref, "x_detrended"),
        )
        np.testing.assert_allclose(
            trend.ravel(),
            _to_array(ref["output"], "trend").ravel(),
            err_msg="Detrend trend mismatch vs MATLAB",
            **_tol(ref, "trend"),
        )


class TestCrossValidationFreqDomainSim:
    """Frequency-domain simulation: reference_freq_domain_sim.json.

    Mirrors validate_reference.m: estimate a BT model of the noiseless
    plant, then filter the input through it via the IFFT path.
    """

    def test_freq_domain_sim(self):
        ref = _load("reference_freq_domain_sim.json")
        y_noiseless = _to_array(ref["input"], "y_noiseless")
        u = _to_array(ref["input"], "u")
        if u.ndim == 1:
            u = u[:, np.newaxis]
        bt_ws = int(ref["params"]["bt_WindowSize"])

        from sid.freq_bt import freq_bt
        from sid._internal.freq_domain_sim import freq_domain_sim

        r_bt = freq_bt(y_noiseless, u, window_size=bt_ws)
        y_pred = freq_domain_sim(r_bt.response, r_bt.frequency, u, u.shape[0])

        np.testing.assert_allclose(
            y_pred.ravel(),
            _to_array(ref["output"], "Y_pred").ravel(),
            err_msg="FreqDomainSim Y_pred mismatch vs MATLAB",
            **_tol(ref, "Y_pred"),
        )


class TestCrossValidationUncertainty:
    """Frequency-domain asymptotic uncertainty: reference_uncertainty.json.

    Exercises sid_uncertainty directly (via freq_bt's outputs), mirroring
    validate_reference.m's sidUncertainty branch.
    """

    def test_uncertainty_std(self):
        ref = _load("reference_uncertainty.json")
        y = _to_array(ref["input"], "y")
        u = _to_array(ref["input"], "u")
        if u.ndim == 1:
            u = u[:, np.newaxis]
        bt_ws = int(ref["params"]["bt_WindowSize"])

        from sid.freq_bt import freq_bt
        from sid._internal.hann_win import hann_win
        from sid._internal.uncertainty import sid_uncertainty

        r_bt = freq_bt(y, u, window_size=bt_ws)
        W = hann_win(bt_ws)
        g_std, phi_v_std = sid_uncertainty(
            r_bt.response, r_bt.noise_spectrum, r_bt.coherence, y.shape[0], W, 1
        )

        np.testing.assert_allclose(
            g_std.ravel(),
            _to_array(ref["output"], "GStd").ravel(),
            err_msg="Uncertainty GStd mismatch vs MATLAB",
            **_tol(ref, "GStd"),
        )
        np.testing.assert_allclose(
            phi_v_std.ravel(),
            _to_array(ref["output"], "PhiVStd").ravel(),
            err_msg="Uncertainty PhiVStd mismatch vs MATLAB",
            **_tol(ref, "PhiVStd"),
        )


class TestCrossValidationLTVdiscIO:
    """Output-COSMIC LTV-IO identification: reference_ltv_io.json."""

    def test_ltv_io(self):
        ref = _load("reference_ltv_io.json")
        Y = _to_array(ref["input"], "Y")
        U = _to_array(ref["input"], "U")
        if U.ndim == 1:
            U = U[:, np.newaxis]
        H = _to_array(ref["input"], "H")
        lam = float(ref["params"]["Lambda"])

        from sid.ltv_disc_io import ltv_disc_io

        result = ltv_disc_io(Y, U, H, lambda_=lam)

        np.testing.assert_allclose(
            result.a.ravel(),
            _to_array(ref["output"], "A").ravel(),
            err_msg="LTV-IO A mismatch vs MATLAB",
            **_tol(ref, "A"),
        )
        np.testing.assert_allclose(
            result.b.ravel(),
            _to_array(ref["output"], "B").ravel(),
            err_msg="LTV-IO B mismatch vs MATLAB",
            **_tol(ref, "B"),
        )
        np.testing.assert_allclose(
            np.asarray(result.cost).ravel(),
            _to_array(ref["output"], "Cost").ravel(),
            err_msg="LTV-IO cost history mismatch vs MATLAB",
            **_tol(ref, "Cost"),
        )
