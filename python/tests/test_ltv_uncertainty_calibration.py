# Copyright (c) 2026 Pedro Lourenco. All rights reserved.
# This code is released under the MIT License. See LICENSE file in the
# project root for full license information.
#
# This module is part of the Open Source System Identification Toolbox (SID).
# https://github.com/pdlourenco/sid

"""COSMIC Bayesian-uncertainty calibration tests (issue #137, SPEC.md §8.9).

The 1/sqrt(N) data scaling (§8.3.2) makes the returned MAP the minimiser of
``||unscaled residual||^2 + N*lambda*||dC||^2``, so the reported posterior must
use the effective prior weight ``N*lambda``: ``P(k) = P_scaled(k)/N``.  Two
tiers, following the seeded-Monte-Carlo convention:

* **Gate smoke** -- deterministic, fixed-seed, small ``M``.  Includes a
  machine-precision unit test that pins ``P`` to the estimator Hessian directly
  (fails the old ``V^T V + lambda F^T F`` formula by up to a factor ``N``), and
  a mid-``lambda`` reported-vs-empirical calibration that fails the old
  ~1.9x over-statement.
* **Campaign** -- larger ``M`` sweep, gated behind ``SID_MC_CAMPAIGN`` so it
  never runs in the normal suite.
"""

from __future__ import annotations

import os

import numpy as np
import pytest

from sid import ltv_disc
from sid._internal.ltv_uncertainty_backward_pass import uncertainty_backward_pass

# --- fixed ground-truth LTV (constant) system -----------------------------
_P, _Q, _N, _L = 2, 1, 80, 3
_A_TRUE = np.array([[0.7, 0.2], [-0.1, 0.6]])
_B_TRUE = np.array([[0.5], [0.3]])
_SIGMA = 0.05


def _make_trajectory(seed: int) -> tuple[np.ndarray, np.ndarray]:
    """One noisy realisation of the true system: states (N+1,p,L), inputs (N,q,L)."""
    rng = np.random.default_rng(seed)
    X = np.zeros((_N + 1, _P, _L))
    U = rng.standard_normal((_N, _Q, _L))
    for traj in range(_L):
        X[0, :, traj] = rng.standard_normal(_P)
        for k in range(_N):
            X[k + 1, :, traj] = (
                _A_TRUE @ X[k, :, traj] + _B_TRUE @ U[k, :, traj] + _SIGMA * rng.standard_normal(_P)
            )
    return X, U


class TestBackwardPassExact:
    """The recursion returns the estimator's posterior, not a mis-weighted one."""

    def test_equals_estimator_hessian_diagonal(self) -> None:
        """P(k) == [A_est^{-1}]_kk with A_est = V^T V + N*lambda*F^T F.

        Constructed brute-force from random SPD scaled blocks, so it is
        engine- and pipeline-independent.  The old reconstruction
        (V^T V + lambda F^T F) fails this by up to a factor N.
        """
        rng = np.random.default_rng(0)
        n, d = 6, 3
        lam = rng.uniform(0.5, 5.0, n - 1)
        S_scaled = np.zeros((d, d, n))
        for k in range(n):
            Ds = rng.standard_normal((5, d))
            reg = lam[0] if k == 0 else lam[n - 2] if k == n - 1 else lam[k - 1] + lam[k]
            S_scaled[:, :, k] = Ds.T @ Ds + reg * np.eye(d)

        # A_est = N * A_scaled: diag N*S_scaled, off-diag -N*lam_k I.
        A = np.zeros((n * d, n * d))
        for k in range(n):
            A[k * d : (k + 1) * d, k * d : (k + 1) * d] = n * S_scaled[:, :, k]
        for k in range(n - 1):
            off = -n * lam[k] * np.eye(d)
            A[k * d : (k + 1) * d, (k + 1) * d : (k + 2) * d] = off
            A[(k + 1) * d : (k + 2) * d, k * d : (k + 1) * d] = off
        A_inv = np.linalg.inv(A)
        P_brute = np.stack(
            [A_inv[k * d : (k + 1) * d, k * d : (k + 1) * d] for k in range(n)], axis=2
        )

        P = uncertainty_backward_pass(S_scaled, lam, n, d)
        np.testing.assert_allclose(P, P_brute, rtol=1e-10, atol=1e-12)

    def test_scales_as_p_scaled_over_n(self) -> None:
        """N * P(k) == [A_scaled^{-1}]_kk (the scaled-Hessian posterior)."""
        rng = np.random.default_rng(1)
        n, d = 5, 2
        lam = rng.uniform(0.5, 3.0, n - 1)
        S_scaled = np.zeros((d, d, n))
        for k in range(n):
            Ds = rng.standard_normal((4, d))
            reg = lam[0] if k == 0 else lam[n - 2] if k == n - 1 else lam[k - 1] + lam[k]
            S_scaled[:, :, k] = Ds.T @ Ds + reg * np.eye(d)

        # A_scaled: diag S_scaled, off-diag -lam_k I (un-inflated).
        A = np.zeros((n * d, n * d))
        for k in range(n):
            A[k * d : (k + 1) * d, k * d : (k + 1) * d] = S_scaled[:, :, k]
        for k in range(n - 1):
            off = -lam[k] * np.eye(d)
            A[k * d : (k + 1) * d, (k + 1) * d : (k + 2) * d] = off
            A[(k + 1) * d : (k + 2) * d, k * d : (k + 1) * d] = off
        A_inv = np.linalg.inv(A)
        P_scaled = np.stack(
            [A_inv[k * d : (k + 1) * d, k * d : (k + 1) * d] for k in range(n)], axis=2
        )
        P = uncertainty_backward_pass(S_scaled, lam, n, d)
        np.testing.assert_allclose(P * n, P_scaled, rtol=1e-10, atol=1e-12)


def _calibration_ratio(lam: float, n_trials: int, seed0: int) -> float:
    """Median over A-entries of reported std / empirical std at the middle step."""
    k0 = _N // 2
    a_hats = np.stack(
        [ltv_disc(*_make_trajectory(seed0 + m), lambda_=lam).a[:, :, k0] for m in range(n_trials)]
    )
    emp_std = a_hats.std(axis=0)
    rep_std = ltv_disc(*_make_trajectory(seed0 - 1), lambda_=lam, uncertainty=True).a_std[:, :, k0]
    return float(np.median(rep_std / emp_std))


class TestReportedStdCalibration:
    """Reported AStd tracks the empirical sampling spread (the user-facing check)."""

    def test_calibrated_at_lcurve_corner(self) -> None:
        """At mid-lambda, reported std ~= empirical std; the old formula over-states.

        Fixed seed, so deterministic.  New code: median ratio ~= 0.89.  The old
        V^T V + lambda F^T F reconstruction gives ~1.64 here (a 1.9x
        over-statement of sigma), which the < 1.2 bound rejects -- this test
        would be decoration if it passed both.
        """
        ratio = _calibration_ratio(lam=1e2, n_trials=120, seed0=1000)
        assert 0.6 < ratio < 1.2, (
            f"reported/empirical std median ratio {ratio:.3f} is out of the "
            "calibrated band [0.6, 1.2] (old #137 code gives ~1.64 here)"
        )


@pytest.mark.skipif(
    not os.environ.get("SID_MC_CAMPAIGN"),
    reason="Monte-Carlo campaign runs only when SID_MC_CAMPAIGN is set (out of the gate).",
)
class TestUncertaintyCampaign:
    """Larger seeded MC sweep across lambda -- run with SID_MC_CAMPAIGN=1."""

    @pytest.mark.parametrize("lam", [1e1, 1e2, 1e3, 1e4])
    def test_calibration_sweep(self, lam: float) -> None:
        ratio = _calibration_ratio(lam=lam, n_trials=500, seed0=50000)
        assert 0.6 < ratio < 1.4, f"lambda={lam:.0e}: ratio {ratio:.3f} out of band"
