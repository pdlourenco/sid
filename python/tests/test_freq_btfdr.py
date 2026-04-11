# Copyright (c) 2026 Pedro Lourenco. All rights reserved.
# This code is released under the MIT License. See LICENSE file in the
# project root for full license information.
#
# This module is part of the Open Source System Identification Toolbox (SID).
# https://github.com/pdlourenco/sid

"""Tests for freq_btfdr from sid.

Port of test_sidFreqBTFDR.m (16 tests).
"""

from __future__ import annotations

import numpy as np
import pytest
from scipy.signal import lfilter

from sid import freq_btfdr, SidError


class TestFreqBTFDR:
    """Unit tests for Blackman-Tukey with frequency-dependent resolution."""

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
        result = freq_btfdr(y, u)
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
    # Test 2: Method identifier
    # ------------------------------------------------------------------
    def test_method(self, siso_data: tuple) -> None:
        """method=='freq_btfdr'."""
        y, u, _N = siso_data
        result = freq_btfdr(y, u)
        assert result.method == "freq_btfdr"

    # ------------------------------------------------------------------
    # Test 3: WindowSize is a vector (per-frequency)
    # ------------------------------------------------------------------
    def test_window_size_is_vector(self, siso_data: tuple) -> None:
        """window_size has shape (nf,)."""
        y, u, _N = siso_data
        result = freq_btfdr(y, u)
        nf = len(result.frequency)
        assert isinstance(result.window_size, np.ndarray)
        assert result.window_size.shape == (nf,)

    # ------------------------------------------------------------------
    # Test 4: Default resolution gives uniform window size
    # ------------------------------------------------------------------
    def test_default_uniform_window(self, siso_data: tuple) -> None:
        """Default R => all Mk are the same."""
        y, u, _N = siso_data
        result = freq_btfdr(y, u)
        Mk = result.window_size
        assert np.all(Mk == Mk[0]), "Default resolution should give uniform window size"

    # ------------------------------------------------------------------
    # Test 5: Custom scalar resolution
    # ------------------------------------------------------------------
    def test_custom_scalar_resolution(self, siso_data: tuple) -> None:
        """R=0.1 => larger M than R=1.0."""
        y, u, _N = siso_data
        result_fine = freq_btfdr(y, u, resolution=0.1)
        result_coarse = freq_btfdr(y, u, resolution=1.0)
        assert result_fine.window_size[0] > result_coarse.window_size[0], (
            "Finer resolution should use larger window"
        )

    # ------------------------------------------------------------------
    # Test 6: Per-frequency resolution vector
    # ------------------------------------------------------------------
    def test_per_frequency_resolution(self, siso_data: tuple) -> None:
        """R_vec from 0.2 to 2.0 => Mk[0] >= Mk[-1]."""
        y, u, _N = siso_data
        nf = 128
        R_vec = np.linspace(0.2, 2.0, nf)
        result = freq_btfdr(y, u, resolution=R_vec)
        Mk = result.window_size
        assert Mk[0] >= Mk[-1], "Larger resolution should give smaller window"

    # ------------------------------------------------------------------
    # Test 7: Error on negative resolution
    # ------------------------------------------------------------------
    def test_error_negative_resolution(self, siso_data: tuple) -> None:
        """R=-0.5 raises SidError with code='bad_resolution'."""
        y, u, _N = siso_data
        with pytest.raises(SidError) as exc_info:
            freq_btfdr(y, u, resolution=-0.5)
        assert exc_info.value.code == "bad_resolution"

    # ------------------------------------------------------------------
    # Test 8: Error on mismatched resolution vector length
    # ------------------------------------------------------------------
    def test_error_mismatched_resolution_length(self, siso_data: tuple) -> None:
        """R=[0.1, 0.2, 0.3] (wrong length) raises SidError."""
        y, u, _N = siso_data
        with pytest.raises(SidError) as exc_info:
            freq_btfdr(y, u, resolution=np.array([0.1, 0.2, 0.3]))
        assert exc_info.value.code == "bad_resolution"

    # ------------------------------------------------------------------
    # Test 9: Time series mode
    # ------------------------------------------------------------------
    def test_time_series(self) -> None:
        """u=None => response is None, coherence is None, noise_spectrum len 128."""
        rng = np.random.default_rng(42)
        y = rng.standard_normal(300)
        result = freq_btfdr(y, None)
        assert result.response is None
        assert result.coherence is None
        assert len(result.noise_spectrum) == 128

    # ------------------------------------------------------------------
    # Test 10: Coherence bounds [0, 1]
    # ------------------------------------------------------------------
    def test_coherence_bounds(self, siso_data: tuple) -> None:
        """All coherence values in [0, 1]."""
        y, u, _N = siso_data
        result = freq_btfdr(y, u)
        assert np.all(result.coherence >= 0)
        assert np.all(result.coherence <= 1)

    # ------------------------------------------------------------------
    # Test 11: Noise spectrum is non-negative
    # ------------------------------------------------------------------
    def test_noise_spectrum_nonneg(self, siso_data: tuple) -> None:
        """Noise spectrum >= 0."""
        y, u, _N = siso_data
        result = freq_btfdr(y, u)
        assert np.all(result.noise_spectrum >= 0)

    # ------------------------------------------------------------------
    # Test 12: MIMO mode
    # ------------------------------------------------------------------
    def test_mimo(self) -> None:
        """y(500x2), u(500x2) => response shape (nf, 2, 2), coherence is None."""
        rng = np.random.default_rng(7)
        N = 500
        u = rng.standard_normal((N, 2))
        y = np.column_stack(
            [
                u[:, 0] + 0.5 * u[:, 1],
                0.3 * u[:, 0] + u[:, 1],
            ]
        ) + 0.1 * rng.standard_normal((N, 2))
        result = freq_btfdr(y, u)
        nf = len(result.frequency)
        assert result.response.shape == (nf, 2, 2)
        assert result.coherence is None

    # ------------------------------------------------------------------
    # Test 13: Custom frequencies
    # ------------------------------------------------------------------
    def test_custom_frequencies(self) -> None:
        """50 custom frequencies => len==50."""
        rng = np.random.default_rng(7)
        N = 500
        u = rng.standard_normal(N)
        y = lfilter([1], [1, -0.8], u)
        w = np.linspace(0.1, np.pi, 50)
        result = freq_btfdr(y, u, frequencies=w)
        assert len(result.frequency) == 50

    # ------------------------------------------------------------------
    # Test 14: Window sizes clamped to [2, N//2]
    # ------------------------------------------------------------------
    def test_window_size_clamped(self) -> None:
        """N=20, very fine R => Mk in [2, N//2]."""
        rng = np.random.default_rng(42)
        N = 20
        u = rng.standard_normal(N)
        y = rng.standard_normal(N)
        result = freq_btfdr(y, u, resolution=0.01)
        assert np.all(result.window_size <= N // 2)
        assert np.all(result.window_size >= 2)

    # ------------------------------------------------------------------
    # Test 15: Known first-order system
    # ------------------------------------------------------------------
    def test_known_first_order(self) -> None:
        """N=5000, noiseless AR(1) a=0.9, relErr < 0.25 at selected freqs."""
        rng = np.random.default_rng(1)
        N = 5000
        u = rng.standard_normal(N)
        y = lfilter([1], [1, -0.9], u)
        result = freq_btfdr(y, u, resolution=0.2)
        w = result.frequency
        G_true = 1.0 / (1.0 - 0.9 * np.exp(-1j * w))
        # 0-based indices corresponding to MATLAB [10, 30, 60, 100]
        idx = [9, 29, 59, 99]
        for k in idx:
            rel_err = abs(abs(result.response[k]) - abs(G_true[k])) / abs(G_true[k])
            assert rel_err < 0.25, f"Magnitude at freq index {k}: relErr={rel_err:.3f}"

    # ------------------------------------------------------------------
    # Test 16: Multi-trajectory -- variance reduction
    # ------------------------------------------------------------------
    def test_multi_trajectory(self) -> None:
        """L=8: num_trajectories==L, variance reduction ratio reasonable."""
        rng = np.random.default_rng(16)
        N = 2000
        L = 8
        u_3d = rng.standard_normal((N, 1, L))
        y_3d = np.zeros((N, 1, L))
        for traj in range(L):
            y_3d[:, 0, traj] = lfilter([1], [1, -0.9], u_3d[:, 0, traj]) + (
                0.3 * rng.standard_normal(N)
            )

        res_mt = freq_btfdr(y_3d, u_3d, resolution=0.2)
        assert res_mt.num_trajectories == L

        res_st = freq_btfdr(y_3d[:, :, 0], u_3d[:, :, 0], resolution=0.2)
        assert res_st.num_trajectories == 1

        ratio = np.median(res_mt.response_std) / np.median(res_st.response_std)
        expected = 1.0 / np.sqrt(L)
        assert ratio < expected * 2.5 and ratio > expected * 0.2, (
            f"Variance reduction ratio {ratio:.3f} should be near {expected:.3f}"
        )
