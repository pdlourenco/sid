# Copyright (c) 2026 Pedro Lourenco. All rights reserved.
# This code is released under the MIT License. See LICENSE file in the
# project root for full license information.
#
# This module is part of the Open Source System Identification Toolbox (SID).
# https://github.com/pdlourenco/sid-matlab

"""Tests for spectrogram from sid.

Port of test_sidSpectrogram.m (18 tests).
"""

from __future__ import annotations

import numpy as np
import pytest

from sid._exceptions import SidError
from sid.spectrogram import spectrogram


class TestSpectrogram:
    """Unit tests for short-time FFT spectrogram."""

    # ------------------------------------------------------------------
    # Test 1: Result struct has all required fields
    # ------------------------------------------------------------------
    def test_result_fields(self) -> None:
        """All SpectrogramResult fields present, method=='spectrogram'."""
        rng = np.random.default_rng(42)
        N = 1000
        x = rng.standard_normal(N)
        result = spectrogram(x)

        required = [
            "time",
            "frequency",
            "frequency_rad",
            "power",
            "power_db",
            "complex_stft",
            "sample_time",
            "window_length",
            "overlap",
            "nfft",
            "method",
        ]
        for field in required:
            assert hasattr(result, field), f"Missing field: {field}"
        assert result.method == "spectrogram"

    # ------------------------------------------------------------------
    # Test 2: Default parameters
    # ------------------------------------------------------------------
    def test_default_parameters(self) -> None:
        """Default window_length=256, overlap=128, nfft=256, sample_time=1.0."""
        rng = np.random.default_rng(42)
        x = rng.standard_normal(1000)
        result = spectrogram(x)

        assert result.window_length == 256
        assert result.overlap == 128
        assert result.nfft == 256
        assert result.sample_time == 1.0

    # ------------------------------------------------------------------
    # Test 3: Output dimensions
    # ------------------------------------------------------------------
    def test_output_dimensions(self) -> None:
        """K segments, n_bins freq bins, shapes match."""
        rng = np.random.default_rng(42)
        N = 1000
        x = rng.standard_normal(N)
        result = spectrogram(x)

        L = 256
        P = 128
        nfft = 256
        step = L - P
        K = (N - L) // step + 1
        n_bins = nfft // 2 + 1

        assert len(result.time) == K
        assert len(result.frequency) == n_bins
        assert result.power.shape[:2] == (n_bins, K)
        assert result.power_db.shape[:2] == (n_bins, K)
        assert result.complex_stft.shape[:2] == (n_bins, K)

    # ------------------------------------------------------------------
    # Test 4: Frequency starts at DC
    # ------------------------------------------------------------------
    def test_frequency_starts_at_dc(self) -> None:
        """frequency[0] == 0."""
        rng = np.random.default_rng(42)
        x = rng.standard_normal(1000)
        result = spectrogram(x)
        assert result.frequency[0] == 0

    # ------------------------------------------------------------------
    # Test 5: Power non-negative
    # ------------------------------------------------------------------
    def test_power_nonneg(self) -> None:
        """All power >= 0."""
        rng = np.random.default_rng(42)
        x = rng.standard_normal(1000)
        result = spectrogram(x)
        assert np.all(result.power >= 0)

    # ------------------------------------------------------------------
    # Test 6: PowerDB consistent with Power
    # ------------------------------------------------------------------
    def test_power_db_consistent(self) -> None:
        """power_db == 10*log10(max(power, eps))."""
        rng = np.random.default_rng(42)
        x = rng.standard_normal(1000)
        result = spectrogram(x)

        expected_db = 10.0 * np.log10(np.maximum(result.power, np.finfo(float).eps))
        np.testing.assert_allclose(result.power_db, expected_db, atol=1e-10)

    # ------------------------------------------------------------------
    # Test 7: Known sinusoid peak at correct frequency
    # ------------------------------------------------------------------
    def test_known_sinusoid(self) -> None:
        """Fs=1000, f0=100 Hz, peak at correct freq for each segment."""
        Fs = 1000.0
        Ts = 1.0 / Fs
        N = 4000
        f0 = 100.0
        t = np.arange(N) * Ts
        x = np.sin(2.0 * np.pi * f0 * t)

        L = 256
        result = spectrogram(x, window_length=L, sample_time=Ts)

        for k in range(len(result.time)):
            peak_idx = np.argmax(result.power[:, k])
            peak_freq = result.frequency[peak_idx]
            assert abs(peak_freq - f0) < Fs / L, (
                f"Segment {k}: peak at {peak_freq} Hz, expected near {f0} Hz"
            )

    # ------------------------------------------------------------------
    # Test 8: Custom parameters
    # ------------------------------------------------------------------
    def test_custom_params(self) -> None:
        """window_length=64, overlap=32, nfft=128 reflected in result."""
        rng = np.random.default_rng(42)
        x = rng.standard_normal(500)
        result = spectrogram(x, window_length=64, overlap=32, nfft=128)

        assert result.window_length == 64
        assert result.overlap == 32
        assert result.nfft == 128
        n_bins = 128 // 2 + 1
        assert result.power.shape[0] == n_bins

    # ------------------------------------------------------------------
    # Test 9: Multi-channel data
    # ------------------------------------------------------------------
    def test_multi_channel(self) -> None:
        """x(1000x3) -> power shape has 3 channels."""
        rng = np.random.default_rng(42)
        x = rng.standard_normal((1000, 3))
        result = spectrogram(x, window_length=128)

        assert result.power.shape[2] == 3
        assert result.complex_stft.shape[2] == 3

    # ------------------------------------------------------------------
    # Test 10: Hamming window
    # ------------------------------------------------------------------
    def test_hamming_window(self) -> None:
        """Power non-negative with 'hamming' window."""
        rng = np.random.default_rng(42)
        x = rng.standard_normal(500)
        result = spectrogram(x, window_length=64, window="hamming")
        assert np.all(result.power >= 0)

    # ------------------------------------------------------------------
    # Test 11: Rectangular window
    # ------------------------------------------------------------------
    def test_rect_window(self) -> None:
        """Power non-negative with 'rect' window."""
        rng = np.random.default_rng(42)
        x = rng.standard_normal(500)
        result = spectrogram(x, window_length=64, window="rect")
        assert np.all(result.power >= 0)

    # ------------------------------------------------------------------
    # Test 12: Custom window vector
    # ------------------------------------------------------------------
    def test_custom_window_vector(self) -> None:
        """w=0.5*ones(64) -> power non-negative."""
        rng = np.random.default_rng(42)
        x = rng.standard_normal(500)
        w = 0.5 * np.ones(64)
        result = spectrogram(x, window_length=64, window=w)
        assert np.all(result.power >= 0)

    # ------------------------------------------------------------------
    # Test 13: Time vector correctness
    # ------------------------------------------------------------------
    def test_time_vector(self) -> None:
        """L=100, P=50, Ts=0.01 -> time matches formula."""
        rng = np.random.default_rng(42)
        L = 100
        P = 50
        Ts = 0.01
        N = 500
        x = rng.standard_normal(N)
        result = spectrogram(x, window_length=L, overlap=P, sample_time=Ts)

        step = L - P
        K = (N - L) // step + 1
        expected_time = (np.arange(K) * step + L / 2) * Ts
        np.testing.assert_allclose(result.time, expected_time, atol=1e-12)

    # ------------------------------------------------------------------
    # Test 14: Error on short data
    # ------------------------------------------------------------------
    def test_error_short_data(self) -> None:
        """N=10, window_length=256 -> SidError code='too_short'."""
        rng = np.random.default_rng(42)
        x = rng.standard_normal(10)
        with pytest.raises(SidError) as exc_info:
            spectrogram(x, window_length=256)
        assert exc_info.value.code == "too_short"

    # ------------------------------------------------------------------
    # Test 15: Error on invalid overlap
    # ------------------------------------------------------------------
    def test_error_invalid_overlap(self) -> None:
        """P >= L -> SidError code='invalid_overlap'."""
        rng = np.random.default_rng(42)
        x = rng.standard_normal(500)
        with pytest.raises(SidError) as exc_info:
            spectrogram(x, window_length=64, overlap=64)
        assert exc_info.value.code == "invalid_overlap"

    # ------------------------------------------------------------------
    # Test 16: FrequencyRad consistent with Frequency
    # ------------------------------------------------------------------
    def test_frequency_rad_consistent(self) -> None:
        """frequency_rad == 2*pi*frequency."""
        rng = np.random.default_rng(42)
        x = rng.standard_normal(1000)
        result = spectrogram(x)
        np.testing.assert_allclose(result.frequency_rad, 2.0 * np.pi * result.frequency, atol=1e-12)

    # ------------------------------------------------------------------
    # Test 17: Multi-trajectory reduces variance
    # ------------------------------------------------------------------
    def test_multi_trajectory(self) -> None:
        """L=10 trajectories, ensemble PSD variance at signal freq < single-traj."""
        rng = np.random.default_rng(17)
        N = 4000
        L_traj = 10
        t = np.arange(N, dtype=np.float64)
        sig = np.sin(2.0 * np.pi * 0.1 * t)

        y = np.zeros((N, 1, L_traj))
        for ll in range(L_traj):
            y[:, 0, ll] = sig + rng.standard_normal(N)

        wlen = 256
        res_mt = spectrogram(y, window_length=wlen, overlap=128)
        res_st = spectrogram(y[:, :, 0], window_length=wlen, overlap=128)

        # At the signal frequency, variance across time segments should be
        # lower for ensemble average (noise cancels, signal is deterministic)
        fbin = np.argmin(np.abs(res_mt.frequency - 0.1))
        var_mt = np.var(res_mt.power[fbin, :])
        var_st = np.var(res_st.power[fbin, :])
        assert var_mt < var_st, (
            f"Ensemble PSD variance {var_mt:.4f} should be < single-traj {var_st:.4f}"
        )

    # ------------------------------------------------------------------
    # Test 18: PSD normalization (integral approx sigma^2)
    # ------------------------------------------------------------------
    def test_psd_normalization(self) -> None:
        """N=1024, rect window, single segment -> integral of PSD ~ sigma^2 (15%)."""
        rng = np.random.default_rng(1800)
        N = 1024
        sigma = 2.0
        y = sigma * rng.standard_normal(N)

        result = spectrogram(y, window_length=N, window="rect", overlap=0)

        # Only one time segment
        assert result.power.shape[1] == 1

        # Frequency resolution: Fs / N, default Fs = 1
        df = 1.0 / N
        psd_integral = np.sum(result.power[:, 0]) * df

        # Should approximate sigma^2 = 4.0
        assert abs(psd_integral - sigma**2) / sigma**2 < 0.15, (
            f"PSD integral {psd_integral:.4f} should approximate sigma^2={sigma**2:.1f}"
        )

        # DC bin should NOT be doubled - verify it's smaller than the mean
        # of interior bins (which are doubled)
        interior_mean = np.mean(result.power[1:-1, 0])
        dc_val = result.power[0, 0]
        assert dc_val < interior_mean, (
            f"DC ({dc_val:.4f}) should be < interior mean ({interior_mean:.4f}) "
            "due to one-sided doubling"
        )
