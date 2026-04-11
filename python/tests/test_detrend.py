# Copyright (c) 2026 Pedro Lourenco. All rights reserved.
# This code is released under the MIT License. See LICENSE file in the
# project root for full license information.
#
# This module is part of the Open Source System Identification Toolbox (SID).
# https://github.com/pdlourenco/sid

"""Tests for detrend from sid.

Port of test_sidDetrend.m (11 tests).
"""

from __future__ import annotations

import numpy as np

from sid.detrend import detrend


class TestDetrend:
    """Unit tests for polynomial detrending."""

    # ------------------------------------------------------------------
    # Test 1: Linear trend removal
    # ------------------------------------------------------------------
    def test_linear_trend_removal(self) -> None:
        """slope=2.5, intercept=10: detrended mean < 0.5, reconstruction exact."""
        rng = np.random.default_rng(42)
        N = 500
        t = np.arange(N, dtype=np.float64)
        slope = 2.5
        intercept = 10.0
        trend_true = slope * t + intercept
        noise = 0.1 * rng.standard_normal(N)
        x = trend_true + noise

        x_dt, trend_est = detrend(x)

        assert abs(np.mean(x_dt)) < 0.5, "Detrended mean should be near zero"
        # Reconstruction: x = x_dt + trend exactly
        np.testing.assert_allclose(x, x_dt + trend_est, atol=1e-12)

    # ------------------------------------------------------------------
    # Test 2: Mean removal (order=0)
    # ------------------------------------------------------------------
    def test_mean_removal(self) -> None:
        """order=0, x=42+randn => mean(x_dt) < 1e-10."""
        rng = np.random.default_rng(5001)
        N = 1000
        x = 42.0 + rng.standard_normal(N)

        x_dm, _ = detrend(x, order=0)
        assert abs(np.mean(x_dm)) < 1e-10, "Mean should be removed to machine precision"

    # ------------------------------------------------------------------
    # Test 3: Quadratic trend
    # ------------------------------------------------------------------
    def test_quadratic_trend(self) -> None:
        """order=2, known quadratic => trend matches within 0.5."""
        rng = np.random.default_rng(42)
        N = 500
        t = np.arange(N, dtype=np.float64)
        trend_true = 0.001 * t**2 - 0.5 * t + 100.0
        x = trend_true + 0.01 * rng.standard_normal(N)

        x_dt, trend_est = detrend(x, order=2)
        assert np.max(np.abs(trend_est - trend_true)) < 0.5, (
            "Quadratic trend should be well approximated"
        )

    # ------------------------------------------------------------------
    # Test 4: Multi-channel
    # ------------------------------------------------------------------
    def test_multi_channel(self) -> None:
        """x(300x3) => output (300,3), reconstruction exact."""
        rng = np.random.default_rng(5002)
        N = 300
        t = np.arange(N, dtype=np.float64)
        x = np.column_stack(
            [
                3.0 * t + 10.0 + rng.standard_normal(N),
                -2.0 * t + 50.0 + rng.standard_normal(N),
                0.5 * t**2 + rng.standard_normal(N),
            ]
        )

        x_dt, trend = detrend(x)
        assert x_dt.shape == (N, 3), "Output should be N x 3"
        np.testing.assert_allclose(x, x_dt + trend, atol=1e-12)

    # ------------------------------------------------------------------
    # Test 5: Segment-wise detrending
    # ------------------------------------------------------------------
    def test_segment_wise(self) -> None:
        """segment_length=300, piecewise linear => each segment mean < 2."""
        rng = np.random.default_rng(5003)
        N = 600
        # Piecewise linear trend: slope changes at t=300
        trend_true = np.zeros(N)
        trend_true[:300] = 2.0 * np.arange(300, dtype=np.float64)
        trend_true[300:] = trend_true[299] - 1.5 * np.arange(300, dtype=np.float64)
        x = trend_true + 0.5 * rng.standard_normal(N)

        x_dt, _ = detrend(x, segment_length=300)
        assert abs(np.mean(x_dt[:300])) < 2, "First segment mean should be small"
        assert abs(np.mean(x_dt[300:])) < 2, "Second segment mean should be small"

    # ------------------------------------------------------------------
    # Test 6: Multi-trajectory (3D)
    # ------------------------------------------------------------------
    def test_multi_trajectory_3d(self) -> None:
        """x(200x1x4) 3D => output preserves shape, reconstruction exact."""
        rng = np.random.default_rng(5004)
        N = 200
        L = 4
        t = np.arange(N, dtype=np.float64)
        x3 = np.zeros((N, 1, L))
        for ll in range(L):
            x3[:, 0, ll] = ((ll + 1) * 0.5) * t + 10.0 * (ll + 1) + rng.standard_normal(N)

        x3_dt, trend3 = detrend(x3)
        assert x3_dt.shape == (N, 1, L), "Output should preserve 3D shape"
        np.testing.assert_allclose(x3, x3_dt + trend3, atol=1e-12)

    # ------------------------------------------------------------------
    # Test 7: Already zero-mean data
    # ------------------------------------------------------------------
    def test_already_zero_mean(self) -> None:
        """order=0, zero-mean data => unchanged within 1e-10."""
        rng = np.random.default_rng(5005)
        N = 200
        x = rng.standard_normal(N)
        x = x - np.mean(x)  # exact zero mean

        x_dt, _ = detrend(x, order=0)
        np.testing.assert_allclose(x, x_dt, atol=1e-10)

    # ------------------------------------------------------------------
    # Test 8: Trend reconstruction
    # ------------------------------------------------------------------
    def test_trend_reconstruction(self) -> None:
        """x = x_dt + trend exactly (1e-12)."""
        rng = np.random.default_rng(5006)
        N = 100
        t = np.arange(N, dtype=np.float64)
        x = 5.0 * t + rng.standard_normal(N)

        x_dt, trend = detrend(x)
        np.testing.assert_allclose(x, x_dt + trend, atol=1e-12)

    # ------------------------------------------------------------------
    # Test 9: High polynomial degree (order=5)
    # ------------------------------------------------------------------
    def test_high_polynomial(self) -> None:
        """order=5, 5th-degree trend => recovered within 0.1."""
        rng = np.random.default_rng(5009)
        N = 200
        t = np.linspace(0.0, 1.0, N)
        trend_true = 3.0 * t**5 - 2.0 * t**4 + t**3 - 0.5 * t**2 + 0.1 * t + 7.0
        x = trend_true + 0.001 * rng.standard_normal(N)

        x_dt, trend_est = detrend(x, order=5)
        assert np.max(np.abs(trend_est - trend_true)) < 0.1, (
            "Order-5 trend should be well recovered"
        )
        np.testing.assert_allclose(x, x_dt + trend_est, atol=1e-12)

    # ------------------------------------------------------------------
    # Test 10: Order >= N clamped gracefully
    # ------------------------------------------------------------------
    def test_order_ge_n_clamped(self) -> None:
        """N=5, order=10 => fits perfectly, detrended near zero."""
        rng = np.random.default_rng(5010)
        N = 5
        x = rng.standard_normal(N)

        x_dt, trend = detrend(x, order=10)
        assert x_dt.shape == (N,), "Output size should be N"
        # With order clamped to N-1, polynomial fits data perfectly
        np.testing.assert_allclose(x_dt, 0.0, atol=1e-8)

    # ------------------------------------------------------------------
    # Test 11: Segment length not dividing N evenly
    # ------------------------------------------------------------------
    def test_segment_not_dividing_n(self) -> None:
        """N=100, segment=30 => last short segment detrended."""
        rng = np.random.default_rng(5011)
        N = 100
        t = np.arange(N, dtype=np.float64)
        x = 2.0 * t + rng.standard_normal(N)

        x_dt, trend = detrend(x, segment_length=30)
        assert x_dt.shape == (N,), "Output size correct"
        np.testing.assert_allclose(x, x_dt + trend, atol=1e-12)
        # Last segment (10 samples, indices 90-99) should still be detrended
        assert abs(np.mean(x_dt[90:])) < 5, "Last short segment should be detrended"
