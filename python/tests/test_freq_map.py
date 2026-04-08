# Copyright (c) 2026 Pedro Lourenco. All rights reserved.
# This code is released under the MIT License. See LICENSE file in the
# project root for full license information.
#
# This module is part of the Open Source System Identification Toolbox (SID).
# https://github.com/pdlourenco/sid-matlab

"""Tests for freq_map from sid.

Port of test_sidFreqMap.m (25 tests).
"""

from __future__ import annotations

import numpy as np
import pytest
from scipy.signal import lfilter

from sid._exceptions import SidError
from sid.freq_map import freq_map
from sid.spectrogram import spectrogram


class TestFreqMapBT:
    """Unit tests for freq_map with Blackman-Tukey (BT) algorithm."""

    # ------------------------------------------------------------------
    # Shared helpers
    # ------------------------------------------------------------------
    @pytest.fixture()
    def siso_data(self) -> tuple[np.ndarray, np.ndarray, int]:
        """Generate SISO dataset with known LTI system, N=2000."""
        rng = np.random.default_rng(42)
        N = 2000
        u = rng.standard_normal(N)
        y = lfilter([1], [1, -0.9], u) + 0.1 * rng.standard_normal(N)
        return y, u, N

    # ------------------------------------------------------------------
    # Test 1: BT result fields
    # ------------------------------------------------------------------
    def test_bt_result_fields(self, siso_data: tuple) -> None:
        """All FreqMapResult fields present, method=='freq_map', algorithm=='bt'."""
        y, u, _N = siso_data
        result = freq_map(y, u)

        required = [
            "time",
            "frequency",
            "frequency_hz",
            "response",
            "response_std",
            "noise_spectrum",
            "noise_spectrum_std",
            "coherence",
            "sample_time",
            "segment_length",
            "overlap",
            "window_size",
            "algorithm",
            "method",
        ]
        for field in required:
            assert hasattr(result, field), f"Missing field: {field}"
        assert result.method == "freq_map"
        assert result.algorithm == "bt"

    # ------------------------------------------------------------------
    # Test 2: BT default parameters
    # ------------------------------------------------------------------
    def test_bt_default_params(self, siso_data: tuple) -> None:
        """segment_length, overlap, window_size match defaults."""
        y, u, N = siso_data
        result = freq_map(y, u)

        L_exp = min(N // 4, 256)
        P_exp = L_exp // 2
        M_exp = min(L_exp // 10, 30)
        assert result.segment_length == L_exp
        assert result.overlap == P_exp
        assert result.window_size == M_exp

    # ------------------------------------------------------------------
    # Test 3: BT dimensions
    # ------------------------------------------------------------------
    def test_bt_dimensions(self, siso_data: tuple) -> None:
        """K segments, nf freqs, response shape (nf, K) SISO."""
        y, u, N = siso_data
        result = freq_map(y, u)

        L = result.segment_length
        P = result.overlap
        step = L - P
        K = (N - L) // step + 1
        nf = len(result.frequency)

        assert len(result.time) == K
        assert result.response.shape == (nf, K)
        assert result.coherence.shape == (nf, K)

    # ------------------------------------------------------------------
    # Test 4: BT time vector
    # ------------------------------------------------------------------
    def test_bt_time_vector(self, siso_data: tuple) -> None:
        """Time vector matches formula."""
        y, u, N = siso_data
        result = freq_map(y, u)

        L = result.segment_length
        P = result.overlap
        step = L - P
        K = (N - L) // step + 1
        expected_time = (np.arange(K) * step + L / 2) * result.sample_time
        np.testing.assert_allclose(result.time, expected_time, atol=1e-12)

    # ------------------------------------------------------------------
    # Test 5: BT coherence bounds
    # ------------------------------------------------------------------
    def test_bt_coherence_bounds(self, siso_data: tuple) -> None:
        """Coherence in [0, 1] (within 1e-10)."""
        y, u, _N = siso_data
        result = freq_map(y, u)

        assert np.all(result.coherence >= -1e-10)
        assert np.all(result.coherence <= 1.0 + 1e-10)

    # ------------------------------------------------------------------
    # Test 6: BT noise spectrum non-negative
    # ------------------------------------------------------------------
    def test_bt_noise_nonneg(self, siso_data: tuple) -> None:
        """noise_spectrum >= -1e-10."""
        y, u, _N = siso_data
        result = freq_map(y, u)
        assert np.all(result.noise_spectrum >= -1e-10)

    # ------------------------------------------------------------------
    # Test 7: BT LTI constancy
    # ------------------------------------------------------------------
    def test_bt_lti_constancy(self, siso_data: tuple) -> None:
        """LTI system -> map roughly constant (median CV < 0.5)."""
        y, u, _N = siso_data
        result = freq_map(y, u)

        mag_map = np.abs(result.response)
        mean_mag = np.mean(mag_map, axis=1)
        std_mag = np.std(mag_map, axis=1)
        cv = std_mag / np.maximum(mean_mag, np.finfo(float).eps)
        assert np.median(cv) < 0.5, (
            f"LTI map should be roughly constant, median CV = {np.median(cv):.2f}"
        )

    # ------------------------------------------------------------------
    # Test 8: BT time series mode
    # ------------------------------------------------------------------
    def test_bt_time_series(self) -> None:
        """u=None -> response is None, coherence is None, noise_spectrum exists."""
        rng = np.random.default_rng(42)
        y = rng.standard_normal(1000)
        result = freq_map(y, None)

        assert result.response is None
        assert result.coherence is None
        assert result.noise_spectrum is not None

    # ------------------------------------------------------------------
    # Test 9: BT custom parameters
    # ------------------------------------------------------------------
    def test_bt_custom_params(self, siso_data: tuple) -> None:
        """segment_length=128, overlap=64, window_size=10 reflected."""
        y, u, _N = siso_data
        result = freq_map(y, u, segment_length=128, overlap=64, window_size=10)

        assert result.segment_length == 128
        assert result.overlap == 64
        assert result.window_size == 10

    # ------------------------------------------------------------------
    # Test 10: BT time alignment with spectrogram
    # ------------------------------------------------------------------
    def test_bt_time_alignment_with_spectrogram(self) -> None:
        """Same L, P, Ts -> time vectors match to 1e-12."""
        rng = np.random.default_rng(7)
        N = 2000
        L = 200
        P = 100
        Ts = 0.01
        u = rng.standard_normal(N)
        y = lfilter([1], [1, -0.8], u) + 0.1 * rng.standard_normal(N)

        map_result = freq_map(y, u, segment_length=L, overlap=P, sample_time=Ts)
        spec_result = spectrogram(y, window_length=L, overlap=P, sample_time=Ts)

        np.testing.assert_allclose(map_result.time, spec_result.time, atol=1e-12)

    # ------------------------------------------------------------------
    # Test 11: BT custom frequencies
    # ------------------------------------------------------------------
    def test_bt_custom_frequencies(self, siso_data: tuple) -> None:
        """64 custom freqs preserved."""
        y, u, _N = siso_data
        freqs = np.linspace(0.1, np.pi, 64)
        result = freq_map(y, u, segment_length=128, frequencies=freqs)

        assert len(result.frequency) == 64
        np.testing.assert_allclose(result.frequency, freqs, atol=1e-12)

    # ------------------------------------------------------------------
    # Test 12: BT error segment too long
    # ------------------------------------------------------------------
    def test_bt_error_segment_too_long(self) -> None:
        """L > N -> SidError code='segment_too_long'."""
        rng = np.random.default_rng(42)
        x = rng.standard_normal(50)
        u = rng.standard_normal(50)
        with pytest.raises(SidError) as exc_info:
            freq_map(x, u, segment_length=100)
        assert exc_info.value.code == "segment_too_long"

    # ------------------------------------------------------------------
    # Test 13: BT error segment too short
    # ------------------------------------------------------------------
    def test_bt_error_segment_too_short(self) -> None:
        """L <= 2*M -> SidError code='segment_too_short'."""
        rng = np.random.default_rng(42)
        N = 2000
        u = rng.standard_normal(N)
        y = lfilter([1], [1, -0.9], u) + 0.1 * rng.standard_normal(N)
        with pytest.raises(SidError) as exc_info:
            freq_map(y, u, segment_length=20, window_size=15)
        assert exc_info.value.code == "segment_too_short"


class TestFreqMapWelch:
    """Unit tests for freq_map with Welch algorithm."""

    # ------------------------------------------------------------------
    # Shared helpers
    # ------------------------------------------------------------------
    @pytest.fixture()
    def welch_data(self) -> tuple[np.ndarray, np.ndarray, int]:
        """Generate SISO dataset for Welch tests, N=4000."""
        rng = np.random.default_rng(50)
        N = 4000
        u = rng.standard_normal(N)
        y = lfilter([1], [1, -0.9], u) + 0.1 * rng.standard_normal(N)
        return y, u, N

    # ------------------------------------------------------------------
    # Test 14: Welch runs
    # ------------------------------------------------------------------
    def test_welch_runs(self, welch_data: tuple) -> None:
        """algorithm='welch' -> correct method, algorithm, window_size, response, coherence."""
        y, u, _N = welch_data
        result = freq_map(y, u, algorithm="welch", segment_length=512)

        assert result.method == "freq_map"
        assert result.algorithm == "welch"
        assert result.window_size is None
        assert result.response is not None
        assert result.coherence is not None

    # ------------------------------------------------------------------
    # Test 15: Welch frequency grid
    # ------------------------------------------------------------------
    def test_welch_freq_grid(self, welch_data: tuple) -> None:
        """Correct number of bins (skip DC), freq <= pi."""
        y, u, _N = welch_data
        result = freq_map(y, u, algorithm="welch", segment_length=512)

        Lsub_default = int(np.floor(512 / 4.5))
        nfft_default = max(256, int(2 ** np.ceil(np.log2(Lsub_default))))
        expected_nf = nfft_default // 2

        assert len(result.frequency) == expected_nf
        assert result.frequency[0] > 0, "Welch frequencies should skip DC"
        assert result.frequency[-1] <= np.pi + 1e-10

    # ------------------------------------------------------------------
    # Test 16: Welch coherence bounds
    # ------------------------------------------------------------------
    def test_welch_coherence_bounds(self, welch_data: tuple) -> None:
        """Coherence in [0, 1]."""
        y, u, _N = welch_data
        result = freq_map(y, u, algorithm="welch", segment_length=512)

        assert np.all(result.coherence >= -1e-10)
        assert np.all(result.coherence <= 1.0 + 1e-10)

    # ------------------------------------------------------------------
    # Test 17: Welch noise spectrum non-negative
    # ------------------------------------------------------------------
    def test_welch_noise_nonneg(self, welch_data: tuple) -> None:
        """noise_spectrum >= -1e-10."""
        y, u, _N = welch_data
        result = freq_map(y, u, algorithm="welch", segment_length=512)
        assert np.all(result.noise_spectrum >= -1e-10)

    # ------------------------------------------------------------------
    # Test 18: Welch time series mode
    # ------------------------------------------------------------------
    def test_welch_time_series(self) -> None:
        """response is None, noise_spectrum exists."""
        rng = np.random.default_rng(42)
        y = rng.standard_normal(2000)
        result = freq_map(y, None, algorithm="welch", segment_length=256)

        assert result.response is None
        assert result.noise_spectrum is not None

    # ------------------------------------------------------------------
    # Test 19: Welch custom sub-segment parameters
    # ------------------------------------------------------------------
    def test_welch_custom_sub_params(self, welch_data: tuple) -> None:
        """sub_segment_length=128, nfft=512 -> 256 freq bins."""
        y, u, _N = welch_data
        result = freq_map(
            y,
            u,
            algorithm="welch",
            segment_length=512,
            sub_segment_length=128,
            sub_overlap=64,
            nfft=512,
            window="hamming",
        )
        assert len(result.frequency) == 256

    # ------------------------------------------------------------------
    # Test 20: Welch LTI constancy
    # ------------------------------------------------------------------
    def test_welch_lti_constancy(self, welch_data: tuple) -> None:
        """Median CV < 0.5."""
        y, u, _N = welch_data
        result = freq_map(y, u, algorithm="welch", segment_length=512)

        mag_map = np.abs(result.response)
        mean_mag = np.mean(mag_map, axis=1)
        std_mag = np.std(mag_map, axis=1)
        cv = std_mag / np.maximum(mean_mag, np.finfo(float).eps)
        assert np.median(cv) < 0.5, (
            f"Welch LTI map should be roughly constant, median CV = {np.median(cv):.2f}"
        )

    # ------------------------------------------------------------------
    # Test 21: Welch time alignment with spectrogram
    # ------------------------------------------------------------------
    def test_welch_time_alignment(self, welch_data: tuple) -> None:
        """Same L, P, Ts -> matches spectrogram time vector."""
        y, u, _N = welch_data
        L = 256
        P = 128
        Ts = 0.01
        map_result = freq_map(y, u, algorithm="welch", segment_length=L, overlap=P, sample_time=Ts)
        spec_result = spectrogram(y, window_length=L, overlap=P, sample_time=Ts)

        np.testing.assert_allclose(map_result.time, spec_result.time, atol=1e-12)


class TestFreqMapGeneral:
    """General tests for freq_map (cross-algorithm, error handling)."""

    # ------------------------------------------------------------------
    # Test 22: Invalid algorithm
    # ------------------------------------------------------------------
    def test_error_invalid_algorithm(self) -> None:
        """'foobar' -> SidError code='invalid_algorithm'."""
        rng = np.random.default_rng(42)
        N = 2000
        u = rng.standard_normal(N)
        y = lfilter([1], [1, -0.9], u) + 0.1 * rng.standard_normal(N)
        with pytest.raises(SidError) as exc_info:
            freq_map(y, u, algorithm="foobar")
        assert exc_info.value.code == "invalid_algorithm"

    # ------------------------------------------------------------------
    # Test 23: BT and Welch correlated on same LTI data
    # ------------------------------------------------------------------
    def test_bt_welch_correlated(self) -> None:
        """BT and Welch on same LTI data -> correlation > 0.8."""
        rng = np.random.default_rng(99)
        N = 4000
        u = rng.standard_normal(N)
        y = lfilter([1, 0.5], [1, -0.8], u) + 0.1 * rng.standard_normal(N)

        res_bt = freq_map(y, u, algorithm="bt", segment_length=512)
        res_wl = freq_map(y, u, algorithm="welch", segment_length=512)

        # Average across time segments (LTI so roughly constant)
        mag_bt = np.mean(np.abs(res_bt.response), axis=1)
        mag_wl = np.mean(np.abs(res_wl.response), axis=1)

        # Interpolate Welch onto BT frequency grid for comparison
        mag_wl_interp = np.interp(res_bt.frequency, res_wl.frequency, mag_wl)

        corr_mat = np.corrcoef(mag_bt, mag_wl_interp)
        assert corr_mat[0, 1] > 0.8, (
            f"BT and Welch magnitude shapes should correlate (r={corr_mat[0, 1]:.2f})"
        )

    # ------------------------------------------------------------------
    # Test 24: BT multi-trajectory
    # ------------------------------------------------------------------
    def test_bt_multi_trajectory(self) -> None:
        """L=6 trajectories, num_trajectories correct, multi-traj std < single-traj."""
        rng = np.random.default_rng(24)
        N = 4000
        L_traj = 6
        u = rng.standard_normal((N, 1, L_traj))
        y = np.zeros((N, 1, L_traj))
        for ll in range(L_traj):
            y[:, 0, ll] = lfilter([1], [1, -0.8], u[:, 0, ll]) + (0.2 * rng.standard_normal(N))

        res_mt = freq_map(y, u, segment_length=512, algorithm="bt")
        assert res_mt.num_trajectories == L_traj

        res_st = freq_map(y[:, :, 0], u[:, :, 0], segment_length=512, algorithm="bt")
        assert res_st.num_trajectories == 1

        # Multi-traj ResponseStd should be smaller
        std_mt = np.median(res_mt.response_std.ravel())
        std_st = np.median(res_st.response_std.ravel())
        assert std_mt < std_st, f"Multi-traj std {std_mt:.4f} should be < single-traj {std_st:.4f}"

    # ------------------------------------------------------------------
    # Test 25: Welch multi-trajectory
    # ------------------------------------------------------------------
    def test_welch_multi_trajectory(self) -> None:
        """L=6 trajectories, num_trajectories correct."""
        rng = np.random.default_rng(24)
        N = 4000
        L_traj = 6
        u = rng.standard_normal((N, 1, L_traj))
        y = np.zeros((N, 1, L_traj))
        for ll in range(L_traj):
            y[:, 0, ll] = lfilter([1], [1, -0.8], u[:, 0, ll]) + (0.2 * rng.standard_normal(N))

        result = freq_map(y, u, segment_length=512, algorithm="welch")
        assert result.num_trajectories == L_traj
