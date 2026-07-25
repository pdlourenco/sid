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

    # ------------------------------------------------------------------
    # Trust-region two-level schedule (SPEC S8.12.4, issue #138)
    # ------------------------------------------------------------------
    @staticmethod
    def _make_partial_obs_data() -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Deterministic partial-observation I/O data (no RNG, cross-platform).

        Partial observation (H = [1 0]) forces the alternating EM path rather
        than the full-rank fast path, so the trust-region schedule is exercised.
        """
        n, q, py, N, L = 2, 1, 1, 30, 8
        A = np.array([[0.9, 0.2], [-0.2, 0.85]])
        B = np.array([[1.0], [0.5]])
        H = np.array([[1.0, 0.0]])
        U = np.zeros((N, q, L))
        Y = np.zeros((N + 1, py, L))
        X = np.zeros((N + 1, n, L))
        for ll in range(L):
            f = (ll + 1) / (3 * N)
            for k in range(N):
                U[k, 0, ll] = np.sin(2 * np.pi * f * (k + 1)) + (1.0 if (k + 1) % 5 < 2 else 0.0)
        for ll in range(L):
            X[0, :, ll] = [0.4 * (ll + 1) / L, -0.3 * (ll + 1) / L]
            Y[0, :, ll] = H @ X[0, :, ll]
            for k in range(N):
                X[k + 1, :, ll] = A @ X[k, :, ll] + B.ravel() * U[k, 0, ll]
                Y[k + 1, :, ll] = H @ X[k + 1, :, ll]
        return Y, U, H

    # A budget large enough for the homotopy stages to make progress. The
    # trust-region benefit is threshold-dependent: too small a per-stage
    # MaxIter and the stages never converge, so TR can match or even trail off
    # (it is a fallback for hard cases, not a free lunch — hence off is the
    # default). At this budget TR beats off by ~2.5 decades on this data.
    _TR_MAX_ITER = 40
    _TR_MU_TOL = 1e-3

    def test_trust_region_helps_on_hard_case(self) -> None:
        """TrustRegion=1 markedly lowers the achieved cost on a hard partial-obs
        case where plain (mu=0) alternation stalls in a poor local minimum.

        Revert-check: the pre-#138 fused loop left TrustRegion=1 ~two orders of
        magnitude WORSE than off, so ``cost_tr < cost_off`` fails against it and
        passes against the corrected two-level schedule. (Note: ``cost_tr <=
        cost_off`` is *not* a general invariant — with too small a per-stage
        MaxIter the homotopy cannot converge and TR may trail off.)
        """
        Y, U, H = self._make_partial_obs_data()

        off = ltv_disc_io(Y, U, H, lambda_=1e2, trust_region="off", max_iter=self._TR_MAX_ITER)
        tr = ltv_disc_io(
            Y,
            U,
            H,
            lambda_=1e2,
            trust_region=1.0,
            max_iter=self._TR_MAX_ITER,
            trust_region_tol=self._TR_MU_TOL,
        )

        cost_off = float(np.min(off.cost))
        cost_tr = float(np.min(tr.cost))

        assert np.all(np.isfinite(tr.a)), "TrustRegion result must be finite"
        assert not np.any(np.isnan(tr.a)), "TrustRegion result must not be NaN"
        assert cost_tr < 0.5 * cost_off, (
            f"TrustRegion should lower the cost markedly: {cost_tr:.3e} vs off {cost_off:.3e}"
        )

    def test_trust_region_mu_advances_and_terminates(self) -> None:
        """The outer loop reduces mu across stages and terminates within budget.

        Iterations must exceed a single per-stage cap (proving mu advanced past
        the initial stage) and stay within the normative worst-case count
        ``MaxIter x (ceil(log2(1/eps_mu)) + 2)`` from SPEC S8.12.4.
        """
        Y, U, H = self._make_partial_obs_data()
        max_iter, mu_tol = self._TR_MAX_ITER, self._TR_MU_TOL

        tr = ltv_disc_io(
            Y,
            U,
            H,
            lambda_=1e2,
            trust_region=1.0,
            max_iter=max_iter,
            trust_region_tol=mu_tol,
        )

        n_stages_worst = int(np.ceil(np.log2(1.0 / mu_tol))) + 2
        bound = max_iter * n_stages_worst
        assert tr.iterations > max_iter, (
            f"Outer loop should run multiple mu stages, got {tr.iterations}"
        )
        assert tr.iterations <= bound, (
            f"Iterations {tr.iterations} exceed worst-case budget {bound}"
        )

    def test_max_iter_rejects_non_positive(self) -> None:
        """MaxIter must be a positive integer (guards the empty-loop crash)."""
        import pytest

        from sid._exceptions import SidError

        Y, U, H = self._make_partial_obs_data()
        with pytest.raises(SidError):
            ltv_disc_io(Y, U, H, lambda_=1e2, max_iter=0)

    def test_trust_region_out_of_range(self) -> None:
        """trust_region must lie in [0, 1] (guards mu>1 dynamics extrapolation)."""
        import pytest

        from sid._exceptions import SidError

        Y, U, H = self._make_partial_obs_data()
        with pytest.raises(SidError):
            ltv_disc_io(Y, U, H, lambda_=1e2, trust_region=2.0)
