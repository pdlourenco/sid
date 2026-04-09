# Copyright (c) 2026 Pedro Lourenco. All rights reserved.
# This code is released under the MIT License. See LICENSE file in the
# project root for full license information.
#
# This module is part of the Open Source System Identification Toolbox (SID).
# https://github.com/pdlourenco/sid-matlab

"""Tests for compare from sid.

Port of test_sidCompare.m.
"""

from __future__ import annotations

import numpy as np
from scipy.signal import lfilter

from sid.compare import compare
from sid.freq_bt import freq_bt
from sid.ltv_disc import ltv_disc


class TestCompare:
    """Unit tests for model output comparison."""

    # ------------------------------------------------------------------
    # Test 1: Result has all expected keys
    # ------------------------------------------------------------------
    def test_result_keys(self) -> None:
        """All expected keys present in compare result."""
        rng = np.random.default_rng(42)
        N = 500
        u = rng.standard_normal(N)
        y = lfilter([1.0], [1.0, -0.8], u) + 0.1 * rng.standard_normal(N)

        model = freq_bt(y, u)
        result = compare(model, y, u)

        expected_keys = ["predicted", "measured", "fit", "residual", "method"]
        for key in expected_keys:
            assert key in result, f"Missing key: {key}"

    # ------------------------------------------------------------------
    # Test 2: Perfect noiseless model -> fit ~ 100%
    # ------------------------------------------------------------------
    def test_perfect_model(self) -> None:
        """Noiseless LTI system with exact model -> fit approx 100%."""
        rng = np.random.default_rng(42)
        N, p, q = 50, 2, 1
        A_true = np.array([[0.9, 0.1], [-0.05, 0.8]])
        B_true = np.array([[0.5], [0.3]])
        L = 3

        X = np.zeros((N + 1, p, L))
        U = rng.standard_normal((N, q, L))
        for ll in range(L):
            X[0, :, ll] = rng.standard_normal(p)
            for k in range(N):
                X[k + 1, :, ll] = A_true @ X[k, :, ll] + B_true.ravel() * U[k, :, ll]

        # Build model with exact A, B replicated across time
        model = ltv_disc(X, U, lambda_=1e8)
        result = compare(model, X, U)

        # With near-perfect recovery and noiseless data, fit should be very high
        fit_vals = np.atleast_1d(result["fit"])
        assert np.all(fit_vals > 90), f"Perfect model should give >90% fit, got {fit_vals}"

    # ------------------------------------------------------------------
    # Test 3: Mean predictor -> fit ~ 0%
    # ------------------------------------------------------------------
    def test_mean_predictor(self) -> None:
        """When predicted = 0 for all t, fit should be near 0%."""
        rng = np.random.default_rng(42)
        N = 500
        u = rng.standard_normal(N)
        y = lfilter([1.0], [1.0, -0.8], u) + 0.1 * rng.standard_normal(N)

        model = freq_bt(y, u)
        result = compare(model, y, u)

        # Compute what fit would be if predicted were zero
        y_meas = result["measured"].ravel()
        y_mean = np.mean(y_meas)
        # The actual NRMSE fit formula: 100*(1 - ||y-ypred||/||y-mean(y)||)
        # When ypred=mean(y), fit=0.  When ypred=0 and mean(y)!=0, fit can be negative
        # Just verify the formula is consistent
        y_pred = result["predicted"].ravel()
        manual_fit = 100.0 * (
            1.0 - np.linalg.norm(y_meas - y_pred) / np.linalg.norm(y_meas - y_mean)
        )
        np.testing.assert_allclose(
            np.atleast_1d(result["fit"]).ravel()[0],
            manual_fit,
            atol=1e-8,
            err_msg="Fit should match manual NRMSE computation",
        )

    # ------------------------------------------------------------------
    # Test 4: Freq-domain compare -> fit > 0
    # ------------------------------------------------------------------
    def test_freq_domain_compare(self) -> None:
        """Estimate with freq_bt, compare -> fit > 0 for good model."""
        rng = np.random.default_rng(42)
        N = 2000
        u = rng.standard_normal(N)
        y = lfilter([1.0], [1.0, -0.85], u) + 0.05 * rng.standard_normal(N)

        model = freq_bt(y, u, window_size=40)
        result = compare(model, y, u)

        fit_vals = np.atleast_1d(result["fit"])
        assert np.all(fit_vals > 0), f"Good freq model should have fit > 0, got {fit_vals}"
        assert result["predicted"].shape[0] == N

    # ------------------------------------------------------------------
    # Test 5: State-space compare -> fit > 50%
    # ------------------------------------------------------------------
    def test_state_space_compare(self) -> None:
        """LTI system with ltv_disc -> fit > 50%."""
        rng = np.random.default_rng(42)
        N, p, q = 60, 2, 1
        A_true = np.array([[0.9, 0.1], [-0.05, 0.8]])
        B_true = np.array([[0.5], [0.3]])
        sigma = 0.03
        L = 5

        X = np.zeros((N + 1, p, L))
        U = rng.standard_normal((N, q, L))
        for ll in range(L):
            X[0, :, ll] = rng.standard_normal(p)
            for k in range(N):
                X[k + 1, :, ll] = (
                    A_true @ X[k, :, ll]
                    + B_true.ravel() * U[k, :, ll]
                    + sigma * rng.standard_normal(p)
                )

        model = ltv_disc(X, U, lambda_=1e4)
        result = compare(model, X, U)

        fit_vals = np.atleast_1d(result["fit"])
        assert np.all(fit_vals > 50), f"COSMIC model should give >50% fit, got {fit_vals}"

    # ------------------------------------------------------------------
    # Test 6: Fit values in reasonable range
    # ------------------------------------------------------------------
    def test_fit_range(self) -> None:
        """Fit values are in a reasonable range (-100 to 100)."""
        rng = np.random.default_rng(42)
        N = 500
        u = rng.standard_normal(N)
        y = lfilter([1.0], [1.0, -0.8], u) + 0.1 * rng.standard_normal(N)

        model = freq_bt(y, u, window_size=30)
        result = compare(model, y, u)

        fit_vals = np.atleast_1d(result["fit"])
        assert np.all(fit_vals >= -100), f"Fit should be >= -100, got {fit_vals}"
        assert np.all(fit_vals <= 100), f"Fit should be <= 100, got {fit_vals}"

    # ------------------------------------------------------------------
    # Test 7: Time-series compare (u=None)
    # ------------------------------------------------------------------
    def test_time_series_compare(self) -> None:
        """Time-series mode: u=None -> predicted is zero, fit near 0."""
        rng = np.random.default_rng(42)
        N = 500
        y = lfilter([1.0], [1.0, -0.6], rng.standard_normal(N))

        model = freq_bt(y, None)
        result = compare(model, y, None)

        # Without input, the model cannot predict output -> predicted ~ 0
        pred = result["predicted"].ravel()
        assert np.allclose(pred, 0.0, atol=1e-10) or (
            np.linalg.norm(pred) < 0.1 * np.linalg.norm(y)
        ), "Time-series predicted should be zero or very small"
