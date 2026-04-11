# Copyright (c) 2026 Pedro Lourenco. All rights reserved.
# This code is released under the MIT License. See LICENSE file in the
# project root for full license information.
#
# This module is part of the Open Source System Identification Toolbox (SID).
# https://github.com/pdlourenco/sid

"""Tests for ltv_state_est from sid.

Port of test_sidLTVStateEst.m (6 tests).
"""

from __future__ import annotations

import numpy as np

# util_msd is the spring-mass-damper plant helper from python/examples/.
# It is made importable here by python/tests/conftest.py, which prepends
# the examples directory to sys.path.
from util_msd import util_msd as _test_msd  # noqa: E402
from sid.ltv_state_est import ltv_state_est


class TestLTVStateEst:
    """Unit tests for batch LTV state estimation (RTS smoother)."""

    # ------------------------------------------------------------------
    # Test 1: Double integrator, full observation
    # ------------------------------------------------------------------
    def test_double_integrator_full_obs(self) -> None:
        """H=I(2), known trajectory with noise, recovered states close."""
        rng = np.random.default_rng(42)
        n, q = 2, 1
        N = 50
        Ts = 1.0
        A_di = np.array([[1.0, Ts], [0.0, 1.0]])
        B_di = np.array([[0.5 * Ts**2], [Ts]])
        H_full = np.eye(2)

        # Known deterministic input trajectory
        U = np.zeros((N, q))
        for k in range(N):
            if k < 20:
                U[k, 0] = 1.0
            elif k < 40:
                U[k, 0] = 0.0
            else:
                U[k, 0] = -1.0

        # True states (noiseless)
        X_true = np.zeros((N + 1, n))
        X_true[0] = [0.0, 0.0]
        for k in range(N):
            X_true[k + 1] = A_di @ X_true[k] + B_di.ravel() * U[k, 0]

        # Noisy measurements
        sigma_meas = 0.05
        Y = X_true @ H_full.T + sigma_meas * rng.standard_normal((N + 1, n))

        # Replicate A, B for all time steps
        A_rep = np.tile(A_di[:, :, np.newaxis], (1, 1, N))
        B_rep = np.tile(B_di[:, :, np.newaxis], (1, 1, N))

        X_hat = ltv_state_est(Y, U, A_rep, B_rep, H_full)
        if isinstance(X_hat, np.ndarray) and X_hat.ndim == 3 and X_hat.shape[2] == 1:
            X_hat = X_hat[:, :, 0]

        # Check recovery
        err = np.linalg.norm(X_hat - X_true, "fro") / max(np.linalg.norm(X_true, "fro"), 1e-10)
        assert err < 0.1, f"Double integrator full obs: relative error {err:.4f} too large"

    # ------------------------------------------------------------------
    # Test 2: Double integrator, partial observation (position only)
    # ------------------------------------------------------------------
    def test_double_integrator_partial_obs(self) -> None:
        """H=[[1,0]] (position only), still recovers velocity."""
        n, q = 2, 1
        N = 50
        Ts = 1.0
        A_di = np.array([[1.0, Ts], [0.0, 1.0]])
        B_di = np.array([[0.5 * Ts**2], [Ts]])
        H_pos = np.array([[1.0, 0.0]])

        # Known deterministic input trajectory
        U = np.zeros((N, q))
        for k in range(N):
            if k < 20:
                U[k, 0] = 1.0
            elif k < 40:
                U[k, 0] = 0.0
            else:
                U[k, 0] = -1.0

        # True states (noiseless)
        X_true = np.zeros((N + 1, n))
        X_true[0] = [0.0, 0.0]
        for k in range(N):
            X_true[k + 1] = A_di @ X_true[k] + B_di.ravel() * U[k, 0]

        # Noiseless position measurements
        Y_pos = X_true @ H_pos.T

        A_rep = np.tile(A_di[:, :, np.newaxis], (1, 1, N))
        B_rep = np.tile(B_di[:, :, np.newaxis], (1, 1, N))

        X_hat = ltv_state_est(Y_pos, U, A_rep, B_rep, H_pos)
        if isinstance(X_hat, np.ndarray) and X_hat.ndim == 3 and X_hat.shape[2] == 1:
            X_hat = X_hat[:, :, 0]

        # Position should be closely recovered
        pos_err = np.linalg.norm(X_hat[:, 0] - X_true[:, 0]) / max(
            np.linalg.norm(X_true[:, 0]), 1e-10
        )
        assert pos_err < 0.01, f"Partial obs: position error {pos_err:.4f}"

        # Velocity should also be recovered via dynamics coupling
        vel_err = np.linalg.norm(X_hat[:, 1] - X_true[:, 1]) / max(
            np.linalg.norm(X_true[:, 1]), 1.0
        )
        assert vel_err < 0.05, f"Partial obs: velocity error {vel_err:.4f}"

    # ------------------------------------------------------------------
    # Test 3: 3-mass MSD system, full observation
    # ------------------------------------------------------------------
    def test_msd_full_obs(self) -> None:
        """3-mass MSD system, H=eye(6), recovered states close to true."""
        m = np.array([1.0, 1.5, 1.0])
        k_spring = np.array([100.0, 80.0, 60.0])
        c_damp = np.array([2.0, 1.5, 1.0])
        F = np.array([[1.0], [0.0], [0.0]])
        Ts = 0.01

        Ad, Bd = _test_msd(m, k_spring, c_damp, F, Ts)
        n = 6
        q = 1
        N = 200
        H_full = np.eye(n)

        # Generate trajectory with sinusoidal input
        U = np.zeros((N, q))
        for k in range(N):
            U[k, 0] = 5.0 * np.sin(2.0 * np.pi * (k + 1) / 50.0)

        X_true = np.zeros((N + 1, n))
        X_true[0] = [0.1, 0.0, -0.05, 0.0, 0.0, 0.0]
        for k in range(N):
            X_true[k + 1] = Ad @ X_true[k] + Bd.ravel() * U[k, 0]

        Y = X_true @ H_full.T

        A_rep = np.tile(Ad[:, :, np.newaxis], (1, 1, N))
        B_rep = np.tile(Bd[:, :, np.newaxis], (1, 1, N))

        X_hat = ltv_state_est(Y, U, A_rep, B_rep, H_full)
        if isinstance(X_hat, np.ndarray) and X_hat.ndim == 3 and X_hat.shape[2] == 1:
            X_hat = X_hat[:, :, 0]

        err = np.linalg.norm(X_hat - X_true, "fro") / max(np.linalg.norm(X_true, "fro"), 1e-10)
        assert err < 1e-4, f"MSD full obs: error {err:.2e}"

    # ------------------------------------------------------------------
    # Test 4: Noisy measurements
    # ------------------------------------------------------------------
    def test_noisy_measurements(self) -> None:
        """Measurement noise with known R; states closer to true than raw."""
        rng = np.random.default_rng(500)
        m = np.array([1.0, 1.5, 1.0])
        k_spring = np.array([100.0, 80.0, 60.0])
        c_damp = np.array([2.0, 1.5, 1.0])
        F = np.array([[1.0], [0.0], [0.0]])
        Ts = 0.01

        Ad, Bd = _test_msd(m, k_spring, c_damp, F, Ts)
        n = 6
        q = 1
        py = 3
        N = 200
        H_pos = np.hstack([np.eye(3), np.zeros((3, 3))])

        U = np.zeros((N, q))
        for k in range(N):
            U[k, 0] = 5.0 * np.sin(2.0 * np.pi * (k + 1) / 50.0)

        X_true = np.zeros((N + 1, n))
        X_true[0] = [0.1, 0.0, -0.05, 0.0, 0.0, 0.0]
        for k in range(N):
            X_true[k + 1] = Ad @ X_true[k] + Bd.ravel() * U[k, 0]

        sigma_meas = 0.01
        Y_noisy = X_true @ H_pos.T + sigma_meas * rng.standard_normal((N + 1, py))
        R_meas = sigma_meas**2 * np.eye(py)
        Q_small = 1e-6 * np.eye(n)

        A_rep = np.tile(Ad[:, :, np.newaxis], (1, 1, N))
        B_rep = np.tile(Bd[:, :, np.newaxis], (1, 1, N))

        X_hat = ltv_state_est(Y_noisy, U, A_rep, B_rep, H_pos, R=R_meas, Q=Q_small)
        if isinstance(X_hat, np.ndarray) and X_hat.ndim == 3 and X_hat.shape[2] == 1:
            X_hat = X_hat[:, :, 0]

        # Check positions are closer to true than raw measurements
        pos_err_hat = np.linalg.norm(X_hat[:, :3] - X_true[:, :3], "fro")
        pos_err_raw = np.linalg.norm(Y_noisy - X_true @ H_pos.T, "fro")
        assert pos_err_hat < pos_err_raw, (
            f"Estimated states ({pos_err_hat:.4f}) should be closer "
            f"to true than raw measurements ({pos_err_raw:.4f})"
        )

    # ------------------------------------------------------------------
    # Test 5: Multi-trajectory
    # ------------------------------------------------------------------
    def test_multi_trajectory(self) -> None:
        """L=3 trajectories -> output shape (N+1, n, 3)."""
        rng = np.random.default_rng(700)
        n, q = 2, 1
        N = 50
        L = 3
        A_true = np.array([[0.9, 0.1], [-0.1, 0.8]])
        B_true = np.array([[0.5], [0.3]])
        H = np.array([[1.0, 0.0]])

        X_true = np.zeros((N + 1, n, L))
        U = rng.standard_normal((N, q, L))
        Y = np.zeros((N + 1, 1, L))

        for ll in range(L):
            X_true[0, :, ll] = rng.standard_normal(n)
            Y[0, :, ll] = H @ X_true[0, :, ll]
            for k in range(N):
                X_true[k + 1, :, ll] = A_true @ X_true[k, :, ll] + B_true.ravel() * U[k, :, ll]
                Y[k + 1, :, ll] = H @ X_true[k + 1, :, ll]

        A_rep = np.tile(A_true[:, :, np.newaxis], (1, 1, N))
        B_rep = np.tile(B_true[:, :, np.newaxis], (1, 1, N))

        X_hat = ltv_state_est(Y, U, A_rep, B_rep, H)
        if isinstance(X_hat, np.ndarray) and X_hat.ndim == 3 and X_hat.shape[2] == 1:
            X_hat = X_hat[:, :, 0]

        assert X_hat.shape == (N + 1, n, L), (
            f"Expected shape ({N + 1}, {n}, {L}), got {X_hat.shape}"
        )

    # ------------------------------------------------------------------
    # Test 6: Default R=I, Q=I produces same result as explicit
    # ------------------------------------------------------------------
    def test_default_r_q(self) -> None:
        """Default R=I, Q=I produces same result as explicit."""
        rng = np.random.default_rng(1100)
        n, q = 2, 1
        N = 30
        A_true = np.array([[0.9, 0.1], [-0.1, 0.8]])
        B_true = np.array([[0.5], [0.3]])
        H = np.array([[1.0, 0.0]])

        U = rng.standard_normal((N, q))
        X_true = np.zeros((N + 1, n))
        X_true[0] = rng.standard_normal(n)
        for k in range(N):
            X_true[k + 1] = A_true @ X_true[k] + B_true.ravel() * U[k, 0]
        Y = X_true @ H.T

        A_rep = np.tile(A_true[:, :, np.newaxis], (1, 1, N))
        B_rep = np.tile(B_true[:, :, np.newaxis], (1, 1, N))

        X_default = ltv_state_est(Y, U, A_rep, B_rep, H)
        X_explicit = ltv_state_est(Y, U, A_rep, B_rep, H, R=np.eye(1), Q=np.eye(n))

        np.testing.assert_allclose(
            X_default,
            X_explicit,
            atol=1e-12,
            err_msg="Default R=I, Q=I should match explicit",
        )
