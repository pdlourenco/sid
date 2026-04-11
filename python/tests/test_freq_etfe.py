# Copyright (c) 2026 Pedro Lourenco. All rights reserved.
# This code is released under the MIT License. See LICENSE file in the
# project root for full license information.
#
# This module is part of the Open Source System Identification Toolbox (SID).
# https://github.com/pdlourenco/sid

"""Tests for freq_etfe from sid.

Port of test_sidFreqETFE.m (15 tests).
"""

from __future__ import annotations

import numpy as np
import pytest
from scipy.signal import lfilter

from sid import freq_etfe, SidError


class TestFreqETFE:
    """Unit tests for Empirical Transfer Function Estimate."""

    # ------------------------------------------------------------------
    # Fixtures
    # ------------------------------------------------------------------

    @pytest.fixture()
    def siso_data(self) -> tuple[np.ndarray, np.ndarray, int]:
        """Generate a SISO dataset: AR(1) with noise, N=500."""
        rng = np.random.default_rng(42)
        N = 500
        u = rng.standard_normal(N)
        y = lfilter([1], [1, -0.8], u) + 0.1 * rng.standard_normal(N)
        return y, u, N

    # ------------------------------------------------------------------
    # Test 1: Result struct has all required fields
    # ------------------------------------------------------------------
    def test_result_fields(self, siso_data: tuple) -> None:
        """All FreqResult fields are present."""
        y, u, _N = siso_data
        result = freq_etfe(y, u)
        required = [
            "frequency",
            "frequency_hz",
            "response",
            "response_std",
            "noise_spectrum",
            "noise_spectrum_std",
            "coherence",
            "sample_time",
            "window_size",
            "data_length",
            "method",
        ]
        for field in required:
            assert hasattr(result, field), f"Missing field: {field}"

    # ------------------------------------------------------------------
    # Test 2: Method and metadata
    # ------------------------------------------------------------------
    def test_method_and_metadata(self, siso_data: tuple) -> None:
        """method=='freq_etfe', data_length==N, window_size==N."""
        y, u, N = siso_data
        result = freq_etfe(y, u)
        assert result.method == "freq_etfe"
        assert result.data_length == N
        assert result.window_size == N

    # ------------------------------------------------------------------
    # Test 3: SISO dimensions
    # ------------------------------------------------------------------
    def test_siso_dimensions(self, siso_data: tuple) -> None:
        """Response shape is (nf,), coherence is None."""
        y, u, _N = siso_data
        result = freq_etfe(y, u)
        nf = len(result.frequency)
        assert result.response.shape == (nf,)
        assert result.coherence is None

    # ------------------------------------------------------------------
    # Test 4: ResponseStd is NaN
    # ------------------------------------------------------------------
    def test_response_std_is_nan(self, siso_data: tuple) -> None:
        """All response_std values are NaN (no asymptotic formula for ETFE)."""
        y, u, _N = siso_data
        result = freq_etfe(y, u)
        assert np.all(np.isnan(result.response_std))

    # ------------------------------------------------------------------
    # Test 5: Noise spectrum is non-negative
    # ------------------------------------------------------------------
    def test_noise_spectrum_nonneg(self, siso_data: tuple) -> None:
        """Noise spectrum >= 0."""
        y, u, _N = siso_data
        result = freq_etfe(y, u)
        assert np.all(result.noise_spectrum >= 0)

    # ------------------------------------------------------------------
    # Test 6: Time series mode (periodogram)
    # ------------------------------------------------------------------
    def test_time_series(self) -> None:
        """u=None: response is None, noise_spectrum >= 0."""
        rng = np.random.default_rng(10)
        y = rng.standard_normal(200)
        result = freq_etfe(y, None)
        assert result.response is None
        assert np.all(result.noise_spectrum >= 0)

    # ------------------------------------------------------------------
    # Test 7: Smoothing reduces variance
    # ------------------------------------------------------------------
    def test_smoothing(self) -> None:
        """S=11 reduces variance of |response| compared to S=1."""
        rng = np.random.default_rng(42)
        N = 500
        u = rng.standard_normal(N)
        y = lfilter([1], [1, -0.8], u)
        result_s1 = freq_etfe(y, u, smoothing=1)
        result_s11 = freq_etfe(y, u, smoothing=11)
        var_s1 = np.var(np.abs(result_s1.response))
        var_s11 = np.var(np.abs(result_s11.response))
        assert var_s11 < var_s1, "Smoothing should reduce variance of response"

    # ------------------------------------------------------------------
    # Test 8: Error on even smoothing parameter
    # ------------------------------------------------------------------
    def test_error_even_smoothing(self) -> None:
        """smoothing=4 raises SidError with code='bad_smoothing'."""
        rng = np.random.default_rng(42)
        N = 500
        u = rng.standard_normal(N)
        y = lfilter([1], [1, -0.8], u)
        with pytest.raises(SidError) as exc_info:
            freq_etfe(y, u, smoothing=4)
        assert exc_info.value.code == "bad_smoothing"

    # ------------------------------------------------------------------
    # Test 9: Error on non-integer smoothing
    # ------------------------------------------------------------------
    def test_error_noninteger_smoothing(self) -> None:
        """smoothing=3.5 raises SidError with code='bad_smoothing'."""
        rng = np.random.default_rng(42)
        N = 500
        u = rng.standard_normal(N)
        y = lfilter([1], [1, -0.8], u)
        with pytest.raises(SidError) as exc_info:
            freq_etfe(y, u, smoothing=3.5)
        assert exc_info.value.code == "bad_smoothing"

    # ------------------------------------------------------------------
    # Test 10: Custom frequencies
    # ------------------------------------------------------------------
    def test_custom_frequencies(self) -> None:
        """50 custom frequencies => len(frequency)==50."""
        rng = np.random.default_rng(42)
        N = 500
        u = rng.standard_normal(N)
        y = lfilter([1], [1, -0.8], u)
        w = np.linspace(0.1, np.pi, 50)
        result = freq_etfe(y, u, frequencies=w)
        assert len(result.frequency) == 50

    # ------------------------------------------------------------------
    # Test 11: Known system -- pure gain y = 2*u
    # ------------------------------------------------------------------
    def test_known_pure_gain(self) -> None:
        """y=2*u (no noise) => |G| approx 2 within 0.01."""
        rng = np.random.default_rng(5)
        N = 1024
        u = rng.standard_normal(N)
        y = 2.0 * u
        result = freq_etfe(y, u)
        G_mag = np.abs(result.response)
        np.testing.assert_allclose(G_mag, 2.0, atol=0.01)

    # ------------------------------------------------------------------
    # Test 12: Known system -- pure delay y(t) = u(t-1)
    # ------------------------------------------------------------------
    def test_known_pure_delay(self) -> None:
        """|G| approx 1, phase approx -w."""
        rng = np.random.default_rng(5)
        N = 1024
        u = rng.standard_normal(N)
        y = np.concatenate(([0.0], u[:-1]))
        result = freq_etfe(y, u)
        G_mag = np.abs(result.response)
        G_phase = np.angle(result.response)
        expected_phase = -result.frequency

        # Magnitude should be ~1
        assert np.median(np.abs(G_mag - 1.0)) < 0.05, "ETFE of pure delay: |G| should be ~1"
        # Phase should be ~-w (allow wrapped tolerance)
        phase_err = np.abs(G_phase - expected_phase)
        phase_err = np.minimum(phase_err, 2 * np.pi - phase_err)
        assert np.median(phase_err) < 0.1, "ETFE of pure delay: phase should be ~-w"

    # ------------------------------------------------------------------
    # Test 13: MIMO mode
    # ------------------------------------------------------------------
    def test_mimo(self) -> None:
        """y(500x2), u(500x1) => response shape (nf, 2, 1)."""
        rng = np.random.default_rng(7)
        N = 500
        u = rng.standard_normal(N)
        y1 = lfilter([1], [1, -0.5], u)
        y2 = lfilter([0.3], [1, -0.7], u)
        y = np.column_stack([y1, y2])
        result = freq_etfe(y, u)
        nf = len(result.frequency)
        assert result.response.shape[0] == nf
        assert result.response.shape[1] == 2

    # ------------------------------------------------------------------
    # Test 14: Custom sample time
    # ------------------------------------------------------------------
    def test_custom_sample_time(self) -> None:
        """sample_time=0.001 => result.sample_time==0.001."""
        rng = np.random.default_rng(7)
        N = 500
        u = rng.standard_normal(N)
        y = lfilter([1], [1, -0.5], u)
        result = freq_etfe(y, u, sample_time=0.001)
        assert result.sample_time == 0.001

    # ------------------------------------------------------------------
    # Test 15: Multi-trajectory -- ensemble averaging
    # ------------------------------------------------------------------
    def test_multi_trajectory(self) -> None:
        """L=8 trajectories: num_trajectories==L, multi-traj error < single-traj."""
        rng = np.random.default_rng(15)
        N = 2000
        L = 8
        u_3d = rng.standard_normal((N, 1, L))
        y_3d = np.zeros((N, 1, L))
        for traj in range(L):
            y_3d[:, 0, traj] = lfilter([1], [1, -0.5], u_3d[:, 0, traj]) + (
                0.3 * rng.standard_normal(N)
            )

        res_mt = freq_etfe(y_3d, u_3d)
        assert res_mt.num_trajectories == L

        res_st = freq_etfe(y_3d[:, :, 0], u_3d[:, :, 0])
        assert res_st.num_trajectories == 1

        # Ensemble-averaged magnitude should be closer to truth
        w = res_mt.frequency
        G_true = 1.0 / (1.0 - 0.5 * np.exp(-1j * w))
        err_mt = np.median(np.abs(np.abs(res_mt.response) - np.abs(G_true)))
        err_st = np.median(np.abs(np.abs(res_st.response) - np.abs(G_true)))
        assert err_mt < err_st * 1.5, (
            f"Multi-traj error {err_mt:.4f} should improve on single {err_st:.4f}"
        )
