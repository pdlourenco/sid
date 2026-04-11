# Copyright (c) 2026 Pedro Lourenco. All rights reserved.
# This code is released under the MIT License. See LICENSE file in the
# project root for full license information.
#
# This module is part of the Open Source System Identification Toolbox (SID).
# https://github.com/pdlourenco/sid-matlab

"""Tests for freq_bt from sid."""

from __future__ import annotations

import numpy as np
import pytest
from scipy.signal import lfilter

from sid import freq_bt, SidError


class TestFreqBT:
    """Unit tests for Blackman-Tukey spectral analysis."""

    @pytest.fixture()
    def siso_data(self) -> tuple[np.ndarray, np.ndarray]:
        """Generate a SISO dataset with a known system."""
        rng = np.random.default_rng(42)
        N = 500
        u = rng.standard_normal(N)
        y = lfilter([1, 0.5], [1, -0.8], u) + 0.1 * rng.standard_normal(N)
        return y, u

    def test_result_fields(self, siso_data: tuple) -> None:
        """Result has all required FreqResult fields."""
        y, u = siso_data
        result = freq_bt(y, u)
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

    def test_metadata(self, siso_data: tuple) -> None:
        """Correct metadata: data_length, method, sample_time, 128 freqs."""
        y, u = siso_data
        result = freq_bt(y, u)
        assert result.data_length == 500
        assert result.method == "freq_bt"
        assert result.sample_time == 1.0
        assert len(result.frequency) == 128

    def test_siso_dimensions(self, siso_data: tuple) -> None:
        """SISO output shapes are all (128,)."""
        y, u = siso_data
        result = freq_bt(y, u)
        nf = len(result.frequency)
        assert result.response.shape == (nf,)
        assert result.response_std.shape == (nf,)
        assert result.noise_spectrum.shape == (nf,)
        assert result.coherence.shape == (nf,)

    def test_coherence_bounds(self, siso_data: tuple) -> None:
        """Coherence is in [0, 1]."""
        y, u = siso_data
        result = freq_bt(y, u)
        assert np.all(result.coherence >= 0)
        assert np.all(result.coherence <= 1)

    def test_noise_spectrum_nonneg(self, siso_data: tuple) -> None:
        """Noise spectrum is non-negative."""
        y, u = siso_data
        result = freq_bt(y, u)
        assert np.all(result.noise_spectrum >= 0)

    def test_time_series(self) -> None:
        """Time series mode: response=None, coherence=None, spectrum >= 0."""
        rng = np.random.default_rng(42)
        y = rng.standard_normal(300)
        result = freq_bt(y, None)
        assert result.response is None
        assert result.coherence is None
        assert np.all(result.noise_spectrum >= -1e-10)

    def test_custom_window_size(self, siso_data: tuple) -> None:
        """Custom window_size=20 is reflected in result."""
        y, u = siso_data
        result = freq_bt(y, u, window_size=20)
        assert result.window_size == 20

    def test_custom_frequencies(self, siso_data: tuple) -> None:
        """50 custom frequencies are returned correctly."""
        y, u = siso_data
        w = np.linspace(0.1, np.pi, 50)
        result = freq_bt(y, u, window_size=20, frequencies=w)
        assert len(result.frequency) == 50
        np.testing.assert_allclose(result.frequency, w, atol=1e-12)

    def test_custom_sample_time(self, siso_data: tuple) -> None:
        """frequency_hz = frequency / (2*pi*Ts)."""
        y, u = siso_data
        Ts = 0.01
        result = freq_bt(y, u, sample_time=Ts)
        expected_hz = result.frequency / (2 * np.pi * Ts)
        np.testing.assert_allclose(result.frequency_hz, expected_hz, atol=1e-12)

    def test_mimo_2x1(self) -> None:
        """MIMO: y(500x2), u(500x1) gives response shape (128, 2, 1)."""
        rng = np.random.default_rng(7)
        N = 500
        u = rng.standard_normal(N)
        y1 = lfilter([1], [1, -0.8], u)
        y2 = lfilter([0.5], [1, -0.5], u)
        y = np.column_stack([y1, y2]) + 0.1 * rng.standard_normal((N, 2))
        result = freq_bt(y, u)
        nf = len(result.frequency)
        assert result.response.shape[0] == nf
        assert result.response.shape[1] == 2
        assert result.coherence is None

    def test_siso_uncertainty(self, siso_data: tuple) -> None:
        """SISO response_std and noise_spectrum_std are finite and >= 0."""
        y, u = siso_data
        result = freq_bt(y, u)
        assert np.all(np.isfinite(result.response_std))
        assert np.all(result.response_std >= 0)
        assert np.all(np.isfinite(result.noise_spectrum_std))

    def test_known_first_order(self) -> None:
        """Noiseless first-order system: magnitude matches at 4 freqs.

        G(z) = 1 / (1 - 0.9*z^{-1}).
        """
        rng = np.random.default_rng(1)
        N = 5000
        u = rng.standard_normal(N)
        y = lfilter([1], [1, -0.9], u)
        result = freq_bt(y, u, window_size=50)
        w = result.frequency
        G_true = 1.0 / (1.0 - 0.9 * np.exp(-1j * w))
        idx = [9, 29, 59, 99]  # 0-based versions of MATLAB [10,30,60,100]
        for k in idx:
            rel_err = abs(abs(result.response[k]) - abs(G_true[k])) / abs(G_true[k])
            assert rel_err < 0.15, f"Magnitude mismatch at freq index {k}: relErr={rel_err:.3f}"

    def test_multi_trajectory(self) -> None:
        """Multi-trajectory: num_trajectories==L, error < single-traj."""
        rng = np.random.default_rng(14)
        N = 2000
        L = 8
        u_3d = rng.standard_normal((N, 1, L))
        y_3d = np.zeros((N, 1, L))
        for traj in range(L):
            y_3d[:, 0, traj] = lfilter([1], [1, -0.9], u_3d[:, 0, traj]) + (
                0.3 * rng.standard_normal(N)
            )

        res_mt = freq_bt(y_3d, u_3d, window_size=50)
        assert res_mt.num_trajectories == L

        res_st = freq_bt(y_3d[:, :, 0], u_3d[:, :, 0], window_size=50)
        assert res_st.num_trajectories == 1

        w = res_mt.frequency
        G_true = 1.0 / (1.0 - 0.9 * np.exp(-1j * w))
        err_mt = np.median(np.abs(np.abs(res_mt.response) - np.abs(G_true)))
        err_st = np.median(np.abs(np.abs(res_st.response) - np.abs(G_true)))
        assert err_mt < err_st, f"Multi-traj error {err_mt:.4f} should be < single {err_st:.4f}"

    def test_error_m_less_than_2(self) -> None:
        """window_size=1 raises SidError."""
        rng = np.random.default_rng(15)
        N = 100
        u = rng.standard_normal(N)
        y = lfilter([1], [1, -0.5], u)
        with pytest.raises(SidError):
            freq_bt(y, u, window_size=1)

    def test_near_zero_input(self) -> None:
        """Near-zero input does not crash."""
        rng = np.random.default_rng(16)
        N = 500
        u = 1e-15 * rng.standard_normal(N)
        y = rng.standard_normal(N)
        result = freq_bt(y, u, window_size=20)
        assert hasattr(result, "response")

    @pytest.mark.parametrize("M", [128, 200, 256, 500])
    def test_large_window_size(self, M: int) -> None:
        """freq_bt produces finite, default-grid output for M in the old crash region.

        This is a regression test for the FFT-fast-path bug where L was hardcoded
        to 2*nf = 256, causing silent positive/negative lag overlap (M in
        [128, 255]) and an IndexError crash (M >= 256). See SPEC.md S2.5.1.
        """
        rng = np.random.default_rng(123)
        N = 2 * M + 100  # ensures freq_bt does not clamp M to N//2
        u = rng.standard_normal(N)
        y = lfilter([1, 0.5], [1, -0.8], u) + 0.05 * rng.standard_normal(N)

        result = freq_bt(y, u, window_size=M)

        # Spec-required output shape on the default 128-point grid
        assert result.window_size == M
        assert result.response.shape == (128,)
        assert np.all(np.isfinite(result.response))
        assert np.all(np.isfinite(result.noise_spectrum))
        assert np.all(result.noise_spectrum >= -1e-10)

        # Cross-check against the direct DFT path. is_default_freqs returns
        # True for any vector that matches the default grid to within 1e-12,
        # so we offset by -1e-11 to defeat that check and force the direct
        # path. The shift must be negative -- a positive shift would push
        # the last element above pi and trip freq_bt's frequency validation.
        # Do NOT "clean this up": the offset is what makes the cross-check
        # exercise the two code paths.
        w_direct = np.arange(1, 129) * np.pi / 128 - 1e-11
        result_direct = freq_bt(y, u, window_size=M, frequencies=w_direct)
        rel_err = np.max(np.abs(result.response - result_direct.response)) / np.max(
            np.abs(result_direct.response)
        )
        assert rel_err < 1e-8, f"M={M}: FFT path vs direct path relErr={rel_err:.2e}"
