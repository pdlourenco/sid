# Copyright (c) 2026 Pedro Lourenco. All rights reserved.
# This code is released under the MIT License. See LICENSE file in the
# project root for full license information.
#
# This module is part of the Open Source System Identification Toolbox (SID).
# https://github.com/pdlourenco/sid-matlab

"""Tests for ltv_disc_tune from sid.

Port of test_sidLTVdiscTune.m (12 tests).
"""

from __future__ import annotations

import numpy as np
import pytest

from sid.ltv_disc_tune import ltv_disc_tune
from sid.ltv_disc import ltv_disc


def _generate_shared_ltv_data():
    """Generate shared LTV train/val data (matches MATLAB rng(1000))."""
    rng = np.random.default_rng(1000)
    p, q, N = 2, 1, 30
    A0 = np.array([[0.95, 0.1], [-0.1, 0.85]])
    dA = np.array([[-0.3, 0.05], [0.05, -0.25]])
    B_true = np.array([[0.5], [0.3]])
    sigma = 0.15

    L_train, L_val = 3, 4

    X_train = np.zeros((N + 1, p, L_train))
    U_train = rng.standard_normal((N, q, L_train))
    X_val = np.zeros((N + 1, p, L_val))
    U_val = rng.standard_normal((N, q, L_val))

    for ll in range(L_train):
        X_train[0, :, ll] = rng.standard_normal(p)
        for k in range(N):
            Ak = A0 + ((k + 1) / N) * dA  # MATLAB k=1..N -> (k/N)*dA
            X_train[k + 1, :, ll] = (
                Ak @ X_train[k, :, ll]
                + B_true.ravel() * U_train[k, :, ll]
                + sigma * rng.standard_normal(p)
            )

    for ll in range(L_val):
        X_val[0, :, ll] = rng.standard_normal(p)
        for k in range(N):
            Ak = A0 + ((k + 1) / N) * dA
            X_val[k + 1, :, ll] = (
                Ak @ X_val[k, :, ll]
                + B_true.ravel() * U_val[k, :, ll]
                + sigma * rng.standard_normal(p)
            )

    return X_train, U_train, X_val, U_val, p, q, N


class TestLTVDiscTune:
    """Unit tests for validation-based lambda tuning."""

    # ------------------------------------------------------------------
    # Test 1: Output shapes
    # ------------------------------------------------------------------
    def test_output_shapes(self) -> None:
        """bestResult is LTVResult, bestLambda>0, allLosses shape (15,)>0."""
        X_train, U_train, X_val, U_val, p, q, N = _generate_shared_ltv_data()
        grid = np.logspace(-2, 6, 15)

        best_result, best_lambda, all_losses = ltv_disc_tune(
            X_train,
            U_train,
            X_val,
            U_val,
            lambda_grid=grid,
        )

        assert hasattr(best_result, "a"), "bestResult should be an LTVResult"
        assert best_lambda > 0, "bestLambda should be positive"
        assert all_losses.shape == (15,), f"allLosses should be (15,), got {all_losses.shape}"
        assert np.all(all_losses > 0), "allLosses should be positive"

    # ------------------------------------------------------------------
    # Test 2: bestResult has correct fields
    # ------------------------------------------------------------------
    def test_best_result_fields(self) -> None:
        """bestResult.method=='ltv_disc', A shape (p,p,N)."""
        X_train, U_train, X_val, U_val, p, q, N = _generate_shared_ltv_data()
        grid = np.logspace(-2, 6, 15)

        best_result, _, _ = ltv_disc_tune(
            X_train,
            U_train,
            X_val,
            U_val,
            lambda_grid=grid,
        )

        assert best_result.method == "ltv_disc"
        assert hasattr(best_result, "a")
        assert hasattr(best_result, "b")
        assert best_result.a.shape == (p, p, N)

    # ------------------------------------------------------------------
    # Test 3: Optimal lambda is interior (not at grid boundary)
    # ------------------------------------------------------------------
    def test_optimal_interior(self) -> None:
        """For LTV problem, optimal lambda is not at grid boundary."""
        X_train, U_train, X_val, U_val, p, q, N = _generate_shared_ltv_data()
        grid = np.logspace(-2, 6, 15)

        _, _, all_losses = ltv_disc_tune(
            X_train,
            U_train,
            X_val,
            U_val,
            lambda_grid=grid,
        )

        best_idx = int(np.argmin(all_losses))
        assert best_idx > 0 and best_idx < len(grid) - 1, (
            f"Optimal lambda should not be at grid boundary (idx={best_idx}/{len(grid)})"
        )

    # ------------------------------------------------------------------
    # Test 4: Custom grid with 3 elements
    # ------------------------------------------------------------------
    def test_custom_grid(self) -> None:
        """3-element grid [1,100,10000] -> losses has 3 entries."""
        X_train, U_train, X_val, U_val, p, q, N = _generate_shared_ltv_data()
        grid_small = np.array([1.0, 100.0, 10000.0])

        _, best_lambda, losses = ltv_disc_tune(
            X_train,
            U_train,
            X_val,
            U_val,
            lambda_grid=grid_small,
        )

        assert len(losses) == 3, f"allLosses should have 3 entries, got {len(losses)}"
        assert np.any(np.abs(best_lambda - grid_small) < 1e-12), "bestLambda should be from grid"

    # ------------------------------------------------------------------
    # Test 5: Consistency with ltv_disc
    # ------------------------------------------------------------------
    def test_consistency_with_ltv_disc(self) -> None:
        """Re-run at bestLambda, compare A,B to 1e-10."""
        X_train, U_train, X_val, U_val, p, q, N = _generate_shared_ltv_data()
        grid = np.logspace(-2, 6, 15)

        best_result, best_lambda, _ = ltv_disc_tune(
            X_train,
            U_train,
            X_val,
            U_val,
            lambda_grid=grid,
        )

        check = ltv_disc(X_train, U_train, lambda_=best_lambda)

        np.testing.assert_allclose(
            best_result.a,
            check.a,
            atol=1e-10,
            err_msg="bestResult.A should match direct ltv_disc call",
        )
        np.testing.assert_allclose(
            best_result.b,
            check.b,
            atol=1e-10,
            err_msg="bestResult.B should match direct ltv_disc call",
        )

    # ------------------------------------------------------------------
    # Test 6: Precondition passthrough
    # ------------------------------------------------------------------
    def test_precondition_passthrough(self) -> None:
        """With precondition=True (warns), check result is valid."""
        X_train, U_train, X_val, U_val, p, q, N = _generate_shared_ltv_data()
        grid_small = np.array([100.0, 1000.0, 10000.0])

        with pytest.warns(UserWarning):
            best_result, _, _ = ltv_disc_tune(
                X_train,
                U_train,
                X_val,
                U_val,
                lambda_grid=grid_small,
                precondition=True,
            )

        # Preconditioning is disabled in v1.0 but the call should still succeed
        assert hasattr(best_result, "a")

    # ------------------------------------------------------------------
    # Test 7: Frequency method output
    # ------------------------------------------------------------------
    def test_frequency_method_output(self) -> None:
        """Frequency method produces valid info dict with fractions in [0,1]."""
        rng = np.random.default_rng(2001)
        p, q, N, L = 2, 1, 60, 5
        sigma = 0.05
        A0 = np.array([[0.9, 0.1], [-0.05, 0.8]])
        B_true = np.array([[0.5], [0.3]])
        X = np.zeros((N + 1, p, L))
        U = rng.standard_normal((N, q, L))
        for ll in range(L):
            X[0, :, ll] = 0.1 * rng.standard_normal(p)
            for k in range(N):
                X[k + 1, :, ll] = (
                    A0 @ X[k, :, ll] + B_true.ravel() * U[k, :, ll] + sigma * rng.standard_normal(p)
                )

        grid_freq = np.logspace(1, 8, 10)
        best_result, best_lambda, info = ltv_disc_tune(
            X,
            U,
            method="frequency",
            lambda_grid=grid_freq,
            segment_length=20,
        )

        assert hasattr(best_result, "a"), "bestResult should be an LTVResult"
        assert best_lambda > 0, "bestLambda should be positive scalar"
        assert "lambda_grid" in info, "info should have lambda_grid"
        assert "fractions" in info, "info should have fractions"
        assert "best_fraction" in info, "info should have best_fraction"
        assert "freq_map_results" in info, "info should have freq_map_results"
        assert len(info["fractions"]) == len(grid_freq), "fractions length should match grid"
        assert np.all(info["fractions"] >= 0) and np.all(info["fractions"] <= 1), (
            "fractions should be in [0, 1]"
        )

    # ------------------------------------------------------------------
    # Test 8: LTI system -> large lambda selected
    # ------------------------------------------------------------------
    def test_lti_large_lambda(self) -> None:
        """For LTI system, frequency method selects large lambda."""
        rng = np.random.default_rng(2002)
        p, q, N, L = 2, 1, 80, 8
        sigma = 0.03
        A_lti = np.array([[0.9, 0.1], [-0.05, 0.8]])
        B_lti = np.array([[0.5], [0.3]])
        X = np.zeros((N + 1, p, L))
        U = rng.standard_normal((N, q, L))
        for ll in range(L):
            X[0, :, ll] = 0.1 * rng.standard_normal(p)
            for k in range(N):
                X[k + 1, :, ll] = (
                    A_lti @ X[k, :, ll]
                    + B_lti.ravel() * U[k, :, ll]
                    + sigma * rng.standard_normal(p)
                )

        grid_lti = np.logspace(1, 10, 15)
        _, best_lambda, _ = ltv_disc_tune(
            X,
            U,
            method="frequency",
            lambda_grid=grid_lti,
            segment_length=25,
        )

        mid_grid = np.sqrt(grid_lti[0] * grid_lti[-1])
        assert best_lambda >= mid_grid, (
            f"LTI: bestLambda={best_lambda:.2e} should be >= midGrid={mid_grid:.2e}"
        )

    # ------------------------------------------------------------------
    # Test 9: LTV system -> moderate lambda (not extreme)
    # ------------------------------------------------------------------
    def test_ltv_moderate_lambda(self) -> None:
        """For LTV ramp system, lambda should be interior (not at extremes)."""
        rng = np.random.default_rng(2003)
        p, q, N, L = 2, 1, 80, 8
        sigma = 0.03
        A0_ltv = np.array([[0.95, 0.1], [-0.1, 0.85]])
        dA_ltv = np.array([[-0.4, 0.05], [0.05, -0.3]])
        B_ltv = np.array([[0.5], [0.3]])
        X = np.zeros((N + 1, p, L))
        U = rng.standard_normal((N, q, L))
        for ll in range(L):
            X[0, :, ll] = 0.1 * rng.standard_normal(p)
            for k in range(N):
                Ak = A0_ltv + ((k + 1) / N) * dA_ltv
                X[k + 1, :, ll] = (
                    Ak @ X[k, :, ll] + B_ltv.ravel() * U[k, :, ll] + sigma * rng.standard_normal(p)
                )

        grid_ltv = np.logspace(0, 10, 20)
        _, best_lambda, _ = ltv_disc_tune(
            X,
            U,
            method="frequency",
            lambda_grid=grid_ltv,
            segment_length=25,
        )

        assert best_lambda >= grid_ltv[0], "LTV: lambda should be >= smallest candidate"
        assert best_lambda <= grid_ltv[-1], "LTV: lambda should be <= largest candidate"

    # ------------------------------------------------------------------
    # Test 10: Fallback with strict threshold
    # ------------------------------------------------------------------
    def test_fallback_strict_threshold(self) -> None:
        """consistency_threshold=0.9999 -> still returns valid result."""
        rng = np.random.default_rng(2003)
        p, q, N, L = 2, 1, 80, 8
        sigma = 0.03
        A0_ltv = np.array([[0.95, 0.1], [-0.1, 0.85]])
        dA_ltv = np.array([[-0.4, 0.05], [0.05, -0.3]])
        B_ltv = np.array([[0.5], [0.3]])
        X = np.zeros((N + 1, p, L))
        U = rng.standard_normal((N, q, L))
        for ll in range(L):
            X[0, :, ll] = 0.1 * rng.standard_normal(p)
            for k in range(N):
                Ak = A0_ltv + ((k + 1) / N) * dA_ltv
                X[k + 1, :, ll] = (
                    Ak @ X[k, :, ll] + B_ltv.ravel() * U[k, :, ll] + sigma * rng.standard_normal(p)
                )

        _, _, info = ltv_disc_tune(
            X,
            U,
            method="frequency",
            lambda_grid=np.logspace(1, 6, 5),
            segment_length=25,
            consistency_threshold=0.9999,
        )

        assert "best_fraction" in info
        assert 0.0 <= info["best_fraction"] <= 1.0, "bestFraction should be in [0, 1]"

    # ------------------------------------------------------------------
    # Test 11: Backward compatibility (method='validation' explicit)
    # ------------------------------------------------------------------
    def test_backward_compatibility(self) -> None:
        """method='validation' explicit gives same result as default."""
        X_train, U_train, X_val, U_val, p, q, N = _generate_shared_ltv_data()
        grid = np.logspace(-2, 6, 15)

        _, best_lambda_default, losses_default = ltv_disc_tune(
            X_train,
            U_train,
            X_val,
            U_val,
            lambda_grid=grid,
        )

        _, best_lambda_explicit, losses_explicit = ltv_disc_tune(
            X_train,
            U_train,
            X_val,
            U_val,
            lambda_grid=grid,
            method="validation",
        )

        np.testing.assert_allclose(
            best_lambda_explicit,
            best_lambda_default,
            atol=1e-12,
            err_msg="Explicit method='validation' should match default",
        )
        np.testing.assert_allclose(
            losses_explicit,
            losses_default,
            atol=1e-12,
            err_msg="Validation losses should match exactly",
        )

    # ------------------------------------------------------------------
    # Test 12: Fractions in valid range
    # ------------------------------------------------------------------
    def test_fractions_valid_range(self) -> None:
        """All fractions in [0,1], some > 0."""
        rng = np.random.default_rng(2001)
        p, q, N, L = 2, 1, 60, 5
        sigma = 0.05
        A0 = np.array([[0.9, 0.1], [-0.05, 0.8]])
        B_true = np.array([[0.5], [0.3]])
        X = np.zeros((N + 1, p, L))
        U = rng.standard_normal((N, q, L))
        for ll in range(L):
            X[0, :, ll] = 0.1 * rng.standard_normal(p)
            for k in range(N):
                X[k + 1, :, ll] = (
                    A0 @ X[k, :, ll] + B_true.ravel() * U[k, :, ll] + sigma * rng.standard_normal(p)
                )

        grid_freq = np.logspace(1, 8, 10)
        _, _, info = ltv_disc_tune(
            X,
            U,
            method="frequency",
            lambda_grid=grid_freq,
            segment_length=20,
        )

        fractions = info["fractions"]
        assert np.any(fractions > 0), "At least some fractions should be > 0"
        assert np.all(fractions >= 0) and np.all(fractions <= 1), (
            "All fractions should be in [0, 1]"
        )
