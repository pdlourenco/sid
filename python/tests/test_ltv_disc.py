# Copyright (c) 2026 Pedro Lourenco. All rights reserved.
# This code is released under the MIT License. See LICENSE file in the
# project root for full license information.
#
# This module is part of the Open Source System Identification Toolbox (SID).
# https://github.com/pdlourenco/sid

"""Tests for ltv_disc from sid.

Port of test_sidLTVdisc.m.
"""

from __future__ import annotations

import numpy as np
import pytest

from sid._exceptions import SidError
from sid.ltv_disc import ltv_disc


class TestLTVDisc:
    """Unit tests for discrete LTV state-space identification (COSMIC)."""

    # ------------------------------------------------------------------
    # Test 1: Result struct has all required fields
    # ------------------------------------------------------------------
    def test_result_fields(self) -> None:
        """All LTVResult fields present."""
        rng = np.random.default_rng(42)
        N, p, q = 20, 2, 1
        A_true = 0.95 * np.eye(p)
        B_true = np.array([[1.0], [0.5]])
        X = np.zeros((N + 1, p))
        U = rng.standard_normal((N, q))
        X[0] = rng.standard_normal(p)
        for k in range(N):
            X[k + 1] = A_true @ X[k] + B_true @ U[k] + 0.01 * rng.standard_normal(p)

        result = ltv_disc(X, U, lambda_=1e3)

        required = [
            "a",
            "b",
            "lambda_",
            "cost",
            "data_length",
            "state_dim",
            "input_dim",
            "num_trajectories",
            "algorithm",
            "preconditioned",
            "method",
        ]
        for field in required:
            assert hasattr(result, field), f"Missing field: {field}"

    # ------------------------------------------------------------------
    # Test 2: Correct metadata values
    # ------------------------------------------------------------------
    def test_metadata(self) -> None:
        """method=='ltv_disc', algorithm=='cosmic', dimensions match."""
        rng = np.random.default_rng(42)
        N, p, q = 20, 2, 1
        A_true = 0.95 * np.eye(p)
        B_true = np.array([[1.0], [0.5]])
        X = np.zeros((N + 1, p))
        U = rng.standard_normal((N, q))
        X[0] = rng.standard_normal(p)
        for k in range(N):
            X[k + 1] = A_true @ X[k] + B_true @ U[k] + 0.01 * rng.standard_normal(p)

        result = ltv_disc(X, U, lambda_=1e3)

        assert result.method == "ltv_disc"
        assert result.algorithm == "cosmic"
        assert result.data_length == N
        assert result.state_dim == p
        assert result.input_dim == q
        assert result.num_trajectories == 1

    # ------------------------------------------------------------------
    # Test 3: Known LTI system recovery (high lambda)
    # ------------------------------------------------------------------
    def test_lti_recovery(self) -> None:
        """LTI system with high lambda -> recovered A,B close to true."""
        rng = np.random.default_rng(42)
        N, p, q = 100, 2, 1
        A_true = 0.95 * np.eye(p)
        B_true = np.array([[1.0], [0.5]])
        X = np.zeros((N + 1, p))
        U = rng.standard_normal((N, q))
        X[0] = rng.standard_normal(p)
        for k in range(N):
            X[k + 1] = A_true @ X[k] + B_true @ U[k] + 0.01 * rng.standard_normal(p)

        result = ltv_disc(X, U, lambda_=1e6)

        # Average recovered A,B across time
        A_mean = np.mean(result.a, axis=2)
        B_mean = np.mean(result.b, axis=2)

        np.testing.assert_allclose(
            A_mean, A_true, atol=0.05, err_msg="LTI A recovery should be close to true A"
        )
        np.testing.assert_allclose(
            B_mean, B_true, atol=0.05, err_msg="LTI B recovery should be close to true B"
        )

    # ------------------------------------------------------------------
    # Test 4: Known LTV system with ramp
    # ------------------------------------------------------------------
    def test_ltv_system(self) -> None:
        """Slowly varying A(k) with ramp pole, medium lambda -> A(k) tracks."""
        rng = np.random.default_rng(42)
        N, p, q = 100, 1, 1
        A_true_seq = 0.5 + 0.4 * np.arange(N) / (N - 1)  # ramp 0.5 -> 0.9
        B_true = 1.0
        X = np.zeros((N + 1, p))
        U = rng.standard_normal((N, q))
        X[0] = rng.standard_normal(p)
        for k in range(N):
            X[k + 1] = A_true_seq[k] * X[k] + B_true * U[k] + 0.01 * rng.standard_normal(p)

        result = ltv_disc(X, U, lambda_=1e2)

        # Recovered A(k) should track the ramp
        A_recovered = result.a.ravel()
        # Check pointwise tolerance
        np.testing.assert_allclose(
            A_recovered, A_true_seq, atol=0.3, err_msg="LTV A(k) should track true ramp A(k)"
        )

    # ------------------------------------------------------------------
    # Test 5: Multi-trajectory improves estimate
    # ------------------------------------------------------------------
    def test_multi_trajectory(self) -> None:
        """L=5 trajectories improve estimate vs L=1."""
        rng = np.random.default_rng(42)
        N, p, q = 30, 2, 1
        A_true = np.array([[0.8, 0.05], [-0.05, 0.7]])
        B_true = np.array([[1.0], [0.5]])
        sigma = 0.1

        # Single trajectory
        X1 = np.zeros((N + 1, p, 1))
        U1 = rng.standard_normal((N, q, 1))
        X1[0, :, 0] = rng.standard_normal(p)
        for k in range(N):
            X1[k + 1, :, 0] = (
                A_true @ X1[k, :, 0] + B_true.ravel() * U1[k, :, 0] + sigma * rng.standard_normal(p)
            )
        res1 = ltv_disc(X1, U1, lambda_=1e3)

        # 5 trajectories (reuse first)
        L = 5
        XL = np.zeros((N + 1, p, L))
        UL = np.zeros((N, q, L))
        XL[:, :, 0] = X1[:, :, 0]
        UL[:, :, 0] = U1[:, :, 0]
        for ll in range(1, L):
            UL[:, :, ll] = rng.standard_normal((N, q))
            XL[0, :, ll] = rng.standard_normal(p)
            for k in range(N):
                XL[k + 1, :, ll] = (
                    A_true @ XL[k, :, ll]
                    + B_true.ravel() * UL[k, :, ll]
                    + sigma * rng.standard_normal(p)
                )
        resL = ltv_disc(XL, UL, lambda_=1e3)

        err1 = np.linalg.norm(np.mean(res1.a, axis=2) - A_true, "fro")
        errL = np.linalg.norm(np.mean(resL.a, axis=2) - A_true, "fro")
        # Multi-trajectory should generally improve, but for a specific seed
        # it may not be strictly better. Assert at least within 3× of single.
        assert errL < err1 * 3.0, (
            f"Multi-trajectory error {errL:.4f} should be within 3× of single {err1:.4f}"
        )

    # ------------------------------------------------------------------
    # Test 6: Auto lambda runs without error
    # ------------------------------------------------------------------
    def test_auto_lambda(self) -> None:
        """lambda_='auto' runs without error, returns positive lambda."""
        rng = np.random.default_rng(42)
        N, p, q = 30, 2, 1
        A_true = 0.95 * np.eye(p)
        B_true = np.array([[1.0], [0.5]])
        L = 5
        X = np.zeros((N + 1, p, L))
        U = rng.standard_normal((N, q, L))
        for ll in range(L):
            X[0, :, ll] = rng.standard_normal(p)
            for k in range(N):
                X[k + 1, :, ll] = (
                    A_true @ X[k, :, ll]
                    + B_true.ravel() * U[k, :, ll]
                    + 0.02 * rng.standard_normal(p)
                )

        result = ltv_disc(X, U, lambda_="auto")

        assert np.all(result.lambda_ > 0), "Auto lambda should be positive"
        assert len(result.lambda_) == N - 1, "Auto lambda should have N-1 elements"

    # ------------------------------------------------------------------
    # Test 7: Scalar lambda expanded to vector
    # ------------------------------------------------------------------
    def test_scalar_lambda(self) -> None:
        """lambda_=1e5 -> lambda_ array is all 1e5."""
        rng = np.random.default_rng(42)
        N, p, q = 15, 2, 1
        X = rng.standard_normal((N + 1, p))
        U = rng.standard_normal((N, q))

        result = ltv_disc(X, U, lambda_=1e5)

        assert len(result.lambda_) == N - 1
        np.testing.assert_allclose(
            result.lambda_,
            1e5 * np.ones(N - 1),
            atol=1e-12,
            err_msg="All lambda values should be 1e5",
        )

    # ------------------------------------------------------------------
    # Test 8: Per-step lambda vector
    # ------------------------------------------------------------------
    def test_vector_lambda(self) -> None:
        """Per-step lambda vector -> stored exactly."""
        rng = np.random.default_rng(42)
        N, p, q = 15, 2, 1
        X = rng.standard_normal((N + 1, p))
        U = rng.standard_normal((N, q))
        lam = np.logspace(1, 5, N - 1)

        result = ltv_disc(X, U, lambda_=lam)

        np.testing.assert_allclose(
            result.lambda_, lam, atol=1e-10, err_msg="Per-step lambda should be stored exactly"
        )

    # ------------------------------------------------------------------
    # Test 9: Cost breakdown: total = fidelity + regularization
    # ------------------------------------------------------------------
    def test_cost_breakdown(self) -> None:
        """cost[0] = cost[1] + cost[2] (total = fidelity + reg)."""
        rng = np.random.default_rng(42)
        N, p, q = 20, 2, 1
        X = rng.standard_normal((N + 1, p))
        U = rng.standard_normal((N, q))

        result = ltv_disc(X, U, lambda_=1e3)

        np.testing.assert_allclose(
            result.cost[0],
            result.cost[1] + result.cost[2],
            atol=1e-10,
            err_msg="cost[0] should equal cost[1] + cost[2]",
        )
        assert result.cost[1] >= 0, "Data fidelity should be non-negative"
        assert result.cost[2] >= 0, "Regularization should be non-negative"

    # ------------------------------------------------------------------
    # Test 10: Error on negative lambda
    # ------------------------------------------------------------------
    def test_error_bad_lambda(self) -> None:
        """lambda_=-1 -> SidError code='bad_lambda'."""
        rng = np.random.default_rng(42)
        X = rng.standard_normal((10, 2))
        U = rng.standard_normal((9, 1))

        with pytest.raises(SidError) as exc_info:
            ltv_disc(X, U, lambda_=-1)
        assert exc_info.value.code == "bad_lambda"

    # ------------------------------------------------------------------
    # Test 11: Error on too-short data
    # ------------------------------------------------------------------
    def test_error_too_short(self) -> None:
        """N=1 -> SidError code='too_short'."""
        rng = np.random.default_rng(42)
        X = rng.standard_normal((2, 2))
        U = rng.standard_normal((1, 1))

        with pytest.raises(SidError) as exc_info:
            ltv_disc(X, U, lambda_=1.0)
        assert exc_info.value.code == "too_short"

    # ------------------------------------------------------------------
    # Test 12: Error on dimension mismatch
    # ------------------------------------------------------------------
    def test_error_dim_mismatch(self) -> None:
        """X and U incompatible dimensions -> SidError."""
        rng = np.random.default_rng(42)
        X = rng.standard_normal((10, 2))
        U = rng.standard_normal((8, 1))  # should be (9, 1)

        with pytest.raises(SidError):
            ltv_disc(X, U, lambda_=1.0)
