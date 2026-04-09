# Copyright (c) 2026 Pedro Lourenco. All rights reserved.
# This code is released under the MIT License. See LICENSE file in the
# project root for full license information.
#
# This module is part of the Open Source System Identification Toolbox (SID).
# https://github.com/pdlourenco/sid-matlab

"""Tests for lti_freq_io from sid.

Port of test_sidLTIfreqIO.m (5 tests).
"""

from __future__ import annotations

import numpy as np

from sid.lti_freq_io import lti_freq_io


class TestLTIFreqIO:
    """Unit tests for LTI frequency-domain I/O identification."""

    # ------------------------------------------------------------------
    # Helpers
    # ------------------------------------------------------------------
    @staticmethod
    def _make_data(
        A_true: np.ndarray,
        B_true: np.ndarray,
        H_obs: np.ndarray,
        N: int,
        L: int,
        sigma: float = 0.0,
        seed: int = 200,
    ) -> tuple[np.ndarray, np.ndarray]:
        """Generate I/O data from a known LTI system."""
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
                    A_true @ X[k, :, ll]
                    + B_true @ U[k, :, ll]
                    + sigma * 0.01 * rng.standard_normal(n)
                )
                Y[k + 1, :, ll] = H_obs @ X[k + 1, :, ll] + sigma * rng.standard_normal(py)

        return Y, U

    # ------------------------------------------------------------------
    # Test 1: Output shapes
    # ------------------------------------------------------------------
    def test_output_shapes(self) -> None:
        """A0 shape (n,n), B0 shape (n,q)."""
        n, q, _py = 3, 2, 2
        N, L = 200, 5
        A_true = np.array(
            [
                [0.8, 0.1, 0.0],
                [-0.1, 0.9, 0.05],
                [0.0, -0.05, 0.85],
            ]
        )
        B_true = np.array([[0.5, 0.0], [0.0, 0.3], [0.1, 0.2]])
        H_obs = np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])

        Y, U = self._make_data(A_true, B_true, H_obs, N, L, seed=100)

        A0, B0 = lti_freq_io(Y, U, H_obs)

        assert A0.shape == (n, n), f"A0 should be ({n},{n}), got {A0.shape}"
        assert B0.shape == (n, q), f"B0 should be ({n},{q}), got {B0.shape}"
        assert np.all(np.isfinite(A0)), "A0 should be finite"
        assert np.all(np.isfinite(B0)), "B0 should be finite"

    # ------------------------------------------------------------------
    # Test 2: Known system, full observation
    # ------------------------------------------------------------------
    def test_known_system_full_obs(self) -> None:
        """H=eye(n), known A,B, noiseless -> A0~A_true, B0~B_true."""
        n, _q = 2, 1
        N, L = 500, 5
        A_true = np.array([[0.9, 0.1], [-0.1, 0.8]])
        B_true = np.array([[0.5], [0.3]])
        H_obs = np.eye(n)

        Y, U = self._make_data(A_true, B_true, H_obs, N, L, sigma=0.0, seed=200)

        A0, B0 = lti_freq_io(Y, U, H_obs)

        # With H=I, basis should match, so compare A directly
        A_err = np.linalg.norm(A0 - A_true, "fro") / np.linalg.norm(A_true, "fro")
        assert A_err < 0.15, f"Known system (H=I): A error {A_err:.4f} too large"

        # Also compare eigenvalues (similarity-invariant)
        eig_true = np.sort(np.abs(np.linalg.eigvals(A_true)))
        eig_est = np.sort(np.abs(np.linalg.eigvals(A0)))
        eig_err = np.linalg.norm(eig_true - eig_est) / np.linalg.norm(eig_true)
        assert eig_err < 0.1, f"Eigenvalue error {eig_err:.4f} too large"

    # ------------------------------------------------------------------
    # Test 3: Stability
    # ------------------------------------------------------------------
    def test_stability(self) -> None:
        """All eigenvalues of A0 have |lambda| < 1."""
        n, _q = 2, 1
        N, L = 200, 3
        A_true = np.array([[0.99, 0.05], [-0.05, 0.98]])
        B_true = np.array([[0.5], [0.3]])
        H_obs = np.eye(n)

        Y, U = self._make_data(A_true, B_true, H_obs, N, L, seed=300)

        A0, _ = lti_freq_io(Y, U, H_obs)

        max_eig = np.max(np.abs(np.linalg.eigvals(A0)))
        assert max_eig < 1.0 + 1e-10, (
            f"All eigenvalues should be inside unit circle, max |eig|={max_eig:.6f}"
        )

    # ------------------------------------------------------------------
    # Test 4: MaxStabilize enforces eigenvalue bound
    # ------------------------------------------------------------------
    def test_max_stabilize(self) -> None:
        """max_stabilize=0.95 -> eigenvalues <= 0.95."""
        n, _q = 2, 1
        N, L = 200, 3
        A_true = np.array([[0.99, 0.05], [-0.05, 0.98]])
        B_true = np.array([[0.5], [0.3]])
        H_obs = np.eye(n)

        Y, U = self._make_data(A_true, B_true, H_obs, N, L, seed=300)

        A0, _ = lti_freq_io(Y, U, H_obs, max_stabilize=0.95)

        max_eig = np.max(np.abs(np.linalg.eigvals(A0)))
        assert max_eig <= 0.95 + 1e-10, f"Eigenvalues should be <= 0.95, max |eig|={max_eig:.6f}"

    # ------------------------------------------------------------------
    # Test 5: Custom horizon runs without error
    # ------------------------------------------------------------------
    def test_custom_horizon(self) -> None:
        """horizon=20 runs without error."""
        n, _q = 2, 1
        N, L = 200, 3
        A_true = np.array([[0.9, 0.1], [-0.1, 0.8]])
        B_true = np.array([[0.5], [0.3]])
        H_obs = np.eye(n)

        Y, U = self._make_data(A_true, B_true, H_obs, N, L, seed=700)

        A0, B0 = lti_freq_io(Y, U, H_obs, horizon=20)

        assert np.all(np.isfinite(A0)), "Horizon=20 produced non-finite A"
        assert np.all(np.isfinite(B0)), "Horizon=20 produced non-finite B"
