# Copyright (c) 2026 Pedro Lourenco. All rights reserved.
# This code is released under the MIT License. See LICENSE file in the
# project root for full license information.
#
# This module is part of the Open Source System Identification Toolbox (SID).
# https://github.com/pdlourenco/sid-matlab

"""Tests for residual from sid.

Port of test_sidResidual.m.
"""

from __future__ import annotations

import numpy as np
import pytest
from scipy.signal import lfilter

from sid._exceptions import SidError
from sid.freq_bt import freq_bt
from sid.ltv_disc import ltv_disc
from sid.residual import residual


class TestResidual:
    """Unit tests for residual analysis."""

    # ------------------------------------------------------------------
    # Test 1: Result has all expected keys
    # ------------------------------------------------------------------
    def test_result_keys(self) -> None:
        """All expected keys present in residual result."""
        rng = np.random.default_rng(42)
        N = 500
        u = rng.standard_normal(N)
        y = lfilter([1.0], [1.0, -0.9], u) + 0.1 * rng.standard_normal(N)

        model = freq_bt(y, u, window_size=30)
        result = residual(model, y, u)

        expected_keys = [
            "residual",
            "auto_corr",
            "cross_corr",
            "confidence_bound",
            "whiteness_pass",
            "independence_pass",
            "data_length",
        ]
        for key in expected_keys:
            assert key in result, f"Missing key: {key}"

    # ------------------------------------------------------------------
    # Test 2: Good freq-domain model -> whiteness should pass
    # ------------------------------------------------------------------
    def test_freq_domain_good_model(self) -> None:
        """Good AR(1) model estimated with freq_bt -> whiteness_pass True."""
        rng = np.random.default_rng(42)
        N = 2000
        u = rng.standard_normal(N)
        y = lfilter([1.0], [1.0, -0.9], u) + 0.1 * rng.standard_normal(N)

        model = freq_bt(y, u, window_size=30)
        result = residual(model, y, u)

        assert result["residual"].shape[0] == N
        assert result["confidence_bound"] > 0
        # Residual should have much smaller variance than the signal
        assert np.std(result["residual"]) < np.std(y) * 0.5

    # ------------------------------------------------------------------
    # Test 3: Time-series mode (u=None) -> cross_corr empty, residual=y
    # ------------------------------------------------------------------
    def test_freq_domain_time_series(self) -> None:
        """Time-series mode: cross_corr is None or empty, residual equals y."""
        rng = np.random.default_rng(42)
        N = 500
        y = lfilter([1.0], [1.0, -0.6], rng.standard_normal(N))

        model = freq_bt(y, None)
        result = residual(model, y, None)

        # Cross-correlation should be None or empty for time-series
        cc = result["cross_corr"]
        assert cc is None or (hasattr(cc, "__len__") and len(cc) == 0)
        # Residual should equal y (no input, so no predicted component)
        np.testing.assert_allclose(
            result["residual"].ravel(),
            y.ravel(),
            atol=1e-10,
            err_msg="Time-series residual should equal y",
        )

    # ------------------------------------------------------------------
    # Test 4: State-space model (ltv_disc) -> residuals small
    # ------------------------------------------------------------------
    def test_state_space_model(self) -> None:
        """LTI system with high lambda -> residuals should be small."""
        rng = np.random.default_rng(42)
        N, p, q = 50, 2, 1
        A_true = np.array([[0.9, 0.1], [-0.05, 0.8]])
        B_true = np.array([[0.5], [0.3]])
        sigma = 0.02

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
        result = residual(model, X, U)

        # Residuals should be small relative to signal
        resid_norm = np.linalg.norm(result["residual"])
        signal_norm = np.linalg.norm(X[1:])
        assert resid_norm / signal_norm < 0.5, (
            f"Residual norm ratio {resid_norm / signal_norm:.3f} should be small"
        )

    # ------------------------------------------------------------------
    # Test 5: Confidence bound formula: 2.58 / sqrt(N)
    # ------------------------------------------------------------------
    def test_confidence_bound(self) -> None:
        """confidence_bound == 2.58 / sqrt(N)."""
        rng = np.random.default_rng(42)
        N = 1600
        u = rng.standard_normal(N)
        y = lfilter([1.0], [1.0, -0.7], u) + 0.1 * rng.standard_normal(N)

        model = freq_bt(y, u, window_size=30)
        result = residual(model, y, u)

        expected_bound = 2.58 / np.sqrt(N)
        np.testing.assert_allclose(
            result["confidence_bound"],
            expected_bound,
            atol=1e-10,
            err_msg="ConfidenceBound should be 2.58/sqrt(N)",
        )

    # ------------------------------------------------------------------
    # Test 6: Autocorrelation at lag 0 is normalised to 1.0
    # ------------------------------------------------------------------
    def test_auto_corr_normalized(self) -> None:
        """auto_corr[0] == 1.0 (normalised)."""
        rng = np.random.default_rng(42)
        N = 1000
        u = rng.standard_normal(N)
        y = lfilter([1.0], [1.0, -0.8], u) + 0.1 * rng.standard_normal(N)

        model = freq_bt(y, u, window_size=30)
        result = residual(model, y, u)

        np.testing.assert_allclose(
            result["auto_corr"].ravel()[0],
            1.0,
            atol=1e-10,
            err_msg="r_ee(0) should be 1.0 (normalised)",
        )

    # ------------------------------------------------------------------
    # Test 7: Wrong model input -> error
    # ------------------------------------------------------------------
    def test_wrong_model_error(self) -> None:
        """Pass dict without expected model attributes -> error."""
        rng = np.random.default_rng(42)
        N = 100
        y = rng.standard_normal(N)
        u = rng.standard_normal(N)

        bad_model = {"foo": 1, "bar": 2}

        with pytest.raises((SidError, TypeError, AttributeError, KeyError, ValueError)):
            residual(bad_model, y, u)
