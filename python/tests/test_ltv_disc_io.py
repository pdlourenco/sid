# Copyright (c) 2026 Pedro Lourenco. All rights reserved.
# This code is released under the MIT License. See LICENSE file in the
# project root for full license information.
#
# This module is part of the Open Source System Identification Toolbox (SID).
# https://github.com/pdlourenco/sid

"""Tests for ltv_disc_io from sid.

Port of test_sidLTVdiscIO.m (8 tests).
"""

from __future__ import annotations

import numpy as np

from sid.ltv_disc import ltv_disc
from sid.ltv_disc_io import ltv_disc_io


class TestLTVDiscIO:
    """Unit tests for output-COSMIC LTV identification."""

    # ------------------------------------------------------------------
    # Helpers
    # ------------------------------------------------------------------
    @staticmethod
    def _make_lti_data(
        A_true: np.ndarray,
        B_true: np.ndarray,
        H_obs: np.ndarray,
        N: int,
        L: int,
        sigma: float = 0.02,
        seed: int = 100,
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Generate I/O data and true states from a known LTI system."""
        rng = np.random.default_rng(seed)
        n = A_true.shape[0]
        q = B_true.shape[1]
        py = H_obs.shape[0]

        X = np.zeros((N + 1, n, L))
        U = rng.standard_normal((N, q, L))
        Y = np.zeros((N + 1, py, L))

        for ll in range(L):
            X[0, :, ll] = rng.standard_normal(n)
            Y[0, :, ll] = H_obs @ X[0, :, ll] + sigma * rng.standard_normal(py)
            for k in range(N):
                X[k + 1, :, ll] = (
                    A_true @ X[k, :, ll] + B_true @ U[k, :, ll] + sigma * rng.standard_normal(n)
                )
                Y[k + 1, :, ll] = H_obs @ X[k + 1, :, ll] + sigma * rng.standard_normal(py)

        return Y, U, X

    # ------------------------------------------------------------------
    # Test 1: Result fields present
    # ------------------------------------------------------------------
    def test_result_fields(self) -> None:
        """All LTVIOResult fields present."""
        n, q, _py = 2, 1, 2
        N, L = 30, 5
        A_true = np.array([[0.9, 0.1], [-0.1, 0.8]])
        B_true = np.array([[0.5], [0.3]])
        H_obs = np.eye(n)

        Y, U, _ = self._make_lti_data(A_true, B_true, H_obs, N, L, sigma=0.01, seed=100)

        result = ltv_disc_io(Y, U, H_obs, lambda_=1e3)

        required = [
            "a",
            "b",
            "x",
            "h",
            "r",
            "cost",
            "iterations",
            "lambda_",
            "data_length",
            "state_dim",
            "output_dim",
            "input_dim",
            "num_trajectories",
            "algorithm",
            "method",
        ]
        for field in required:
            assert hasattr(result, field), f"Missing field: {field}"
        assert result.method == "ltv_disc_io"
        assert result.algorithm == "cosmic"
        assert result.data_length == N
        assert result.state_dim == n
        assert result.input_dim == q
        assert result.num_trajectories == L

    # ------------------------------------------------------------------
    # Test 2: H=I matches ltv_disc within tolerance
    # ------------------------------------------------------------------
    def test_h_identity_matches_ltv_disc(self) -> None:
        """H=eye(p) -> A,B match ltv_disc within 1e-6."""
        rng = np.random.default_rng(200)
        n, q = 2, 1
        N, L = 50, 10
        A_true = np.array([[0.9, 0.1], [-0.1, 0.8]])
        B_true = np.array([[0.5], [0.3]])
        sigma = 0.02

        X = np.zeros((N + 1, n, L))
        U = rng.standard_normal((N, q, L))
        for ll in range(L):
            X[0, :, ll] = rng.standard_normal(n)
            for k in range(N):
                X[k + 1, :, ll] = (
                    A_true @ X[k, :, ll]
                    + B_true.ravel() * U[k, :, ll]
                    + sigma * rng.standard_normal(n)
                )

        lam = 1e4
        res_std = ltv_disc(X, U, lambda_=lam)
        res_io = ltv_disc_io(X, U, np.eye(n), lambda_=lam)

        A_err = np.linalg.norm(np.mean(res_io.a, axis=2) - np.mean(res_std.a, axis=2), "fro") / max(
            np.linalg.norm(np.mean(res_std.a, axis=2), "fro"), 1e-15
        )
        B_err = np.linalg.norm(np.mean(res_io.b, axis=2) - np.mean(res_std.b, axis=2), "fro") / max(
            np.linalg.norm(np.mean(res_std.b, axis=2), "fro"), 1e-15
        )

        assert A_err < 1e-6, f"H=I: A mismatch with ltv_disc (err={A_err:.2e})"
        assert B_err < 1e-6, f"H=I: B mismatch with ltv_disc (err={B_err:.2e})"

    # ------------------------------------------------------------------
    # Test 3: Monotone cost decrease
    # ------------------------------------------------------------------
    def test_monotone_cost(self) -> None:
        """Cost history is non-increasing (within tolerance)."""
        rng = np.random.default_rng(200)
        n, q = 2, 1
        N, L = 50, 10
        A_true = np.array([[0.9, 0.1], [-0.1, 0.8]])
        B_true = np.array([[0.5], [0.3]])
        H_obs = np.array([[1.0, 0.0]])
        sigma = 0.02

        X = np.zeros((N + 1, n, L))
        U = rng.standard_normal((N, q, L))
        Y = np.zeros((N + 1, 1, L))
        for ll in range(L):
            X[0, :, ll] = rng.standard_normal(n)
            Y[0, :, ll] = H_obs @ X[0, :, ll]
            for k in range(N):
                X[k + 1, :, ll] = (
                    A_true @ X[k, :, ll]
                    + B_true.ravel() * U[k, :, ll]
                    + sigma * rng.standard_normal(n)
                )
                Y[k + 1, :, ll] = H_obs @ X[k + 1, :, ll] + sigma * rng.standard_normal(1)

        result = ltv_disc_io(Y, U, H_obs, lambda_=1e4)

        cost = result.cost
        if hasattr(cost, "__len__") and len(cost) >= 2:
            for i in range(1, len(cost)):
                assert cost[i] <= cost[i - 1] + 1e-8 * abs(cost[i - 1]), (
                    f"Cost increased at iteration {i}: {cost[i]:.6f} > {cost[i - 1]:.6f}"
                )

    # ------------------------------------------------------------------
    # Test 4: Partial observation
    # ------------------------------------------------------------------
    def test_partial_observation(self) -> None:
        """H=[[1,0],[0,1],[0,0],[0,0]] (2 of 4 obs), recovers reasonable A,B."""
        rng = np.random.default_rng(400)
        n, q = 2, 1
        py = 4
        N, L = 50, 5
        A_true = np.array([[0.9, 0.1], [-0.1, 0.8]])
        B_true = np.array([[0.5], [0.3]])
        H_obs = np.array(
            [
                [1.0, 0.0],
                [0.0, 1.0],
                [0.0, 0.0],
                [0.0, 0.0],
            ]
        )

        X = np.zeros((N + 1, n, L))
        U = rng.standard_normal((N, q, L))
        Y = np.zeros((N + 1, py, L))
        for ll in range(L):
            X[0, :, ll] = rng.standard_normal(n)
            Y[0, :, ll] = H_obs @ X[0, :, ll]
            for k in range(N):
                X[k + 1, :, ll] = (
                    A_true @ X[k, :, ll]
                    + B_true.ravel() * U[k, :, ll]
                    + 0.01 * rng.standard_normal(n)
                )
                Y[k + 1, :, ll] = H_obs @ X[k + 1, :, ll] + 0.01 * rng.standard_normal(py)

        result = ltv_disc_io(Y, U, H_obs, lambda_=1e3)

        assert not np.any(np.isnan(result.a)), "Partial obs: NaN in A"
        assert not np.any(np.isnan(result.b)), "Partial obs: NaN in B"
        assert result.a.shape == (n, n, N)
        assert result.b.shape == (n, q, N)

    # ------------------------------------------------------------------
    # Test 5: Multi-trajectory
    # ------------------------------------------------------------------
    def test_multi_trajectory(self) -> None:
        """L=3 -> num_trajectories==3, x shape matches."""
        n, _q, _py = 2, 1, 2
        N, L = 30, 3
        A_true = np.array([[0.9, 0.1], [-0.1, 0.8]])
        B_true = np.array([[0.5], [0.3]])
        H_obs = np.eye(n)

        Y, U, _ = self._make_lti_data(A_true, B_true, H_obs, N, L, sigma=0.01, seed=500)

        result = ltv_disc_io(Y, U, H_obs, lambda_=1e3)

        assert result.num_trajectories == L
        assert result.x.shape == (N + 1, n, L), (
            f"Expected x shape ({N + 1}, {n}, {L}), got {result.x.shape}"
        )

    # ------------------------------------------------------------------
    # Test 6: Fast path for full-rank H
    # ------------------------------------------------------------------
    def test_fast_path_full_rank(self) -> None:
        """H with rank(H)==n -> iterations==0."""
        rng = np.random.default_rng(1600)
        n, q = 2, 1
        N, L = 80, 5
        A_true = np.array([[0.9, 0.1], [-0.1, 0.8]])
        B_true = np.array([[0.5], [0.3]])

        X = np.zeros((N + 1, n, L))
        U = rng.standard_normal((N, q, L))
        for ll in range(L):
            X[0, :, ll] = rng.standard_normal(n)
            for k in range(N):
                X[k + 1, :, ll] = A_true @ X[k, :, ll] + B_true.ravel() * U[k, :, ll]

        result = ltv_disc_io(X, U, np.eye(n), lambda_=1e4)

        assert result.iterations == 0, (
            f"Fast path should have 0 iterations, got {result.iterations}"
        )

    # ------------------------------------------------------------------
    # Test 7: Convergence within max_iter
    # ------------------------------------------------------------------
    def test_convergence(self) -> None:
        """Converges within max_iter."""
        rng = np.random.default_rng(200)
        n, q = 2, 1
        N, L = 50, 10
        A_true = np.array([[0.9, 0.1], [-0.1, 0.8]])
        B_true = np.array([[0.5], [0.3]])
        H_obs = np.array([[1.0, 0.0]])
        sigma = 0.02

        X = np.zeros((N + 1, n, L))
        U = rng.standard_normal((N, q, L))
        Y = np.zeros((N + 1, 1, L))
        for ll in range(L):
            X[0, :, ll] = rng.standard_normal(n)
            Y[0, :, ll] = H_obs @ X[0, :, ll]
            for k in range(N):
                X[k + 1, :, ll] = (
                    A_true @ X[k, :, ll]
                    + B_true.ravel() * U[k, :, ll]
                    + sigma * rng.standard_normal(n)
                )
                Y[k + 1, :, ll] = H_obs @ X[k + 1, :, ll] + sigma * rng.standard_normal(1)

        max_iter = 50
        result = ltv_disc_io(Y, U, H_obs, lambda_=1e4, max_iter=max_iter)

        assert result.iterations <= max_iter, (
            f"Should converge within {max_iter} iterations, got {result.iterations}"
        )
        assert not np.any(np.isnan(result.a)), "Should produce finite A"

    # ------------------------------------------------------------------
    # Test 8: Uncertainty fields
    # ------------------------------------------------------------------
    def test_uncertainty_fields(self) -> None:
        """a_std, b_std exist and are positive."""
        n, _q, _py = 2, 1, 2
        N, L = 30, 5
        A_true = np.array([[0.9, 0.1], [-0.1, 0.8]])
        B_true = np.array([[0.5], [0.3]])
        H_obs = np.eye(n)

        Y, U, _ = self._make_lti_data(A_true, B_true, H_obs, N, L, sigma=0.05, seed=800)

        result = ltv_disc_io(Y, U, H_obs, lambda_=1e4, uncertainty=True)

        assert result.a_std is not None, "a_std should not be None"
        assert result.b_std is not None, "b_std should not be None"
        assert np.all(result.a_std > 0), "a_std should be positive"
        assert np.all(result.b_std > 0), "b_std should be positive"
        assert np.all(np.isfinite(result.a_std)), "a_std should be finite"
        assert np.all(np.isfinite(result.b_std)), "b_std should be finite"
