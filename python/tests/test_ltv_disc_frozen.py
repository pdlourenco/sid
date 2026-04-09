# Copyright (c) 2026 Pedro Lourenco. All rights reserved.
# This code is released under the MIT License. See LICENSE file in the
# project root for full license information.
#
# This module is part of the Open Source System Identification Toolbox (SID).
# https://github.com/pdlourenco/sid-matlab

"""Tests for ltv_disc_frozen from sid.

Port of test_sidLTVdiscFrozen.m (8 tests).
"""

from __future__ import annotations

import numpy as np

from sid.ltv_disc import ltv_disc
from sid.ltv_disc_frozen import ltv_disc_frozen


class TestLTVDiscFrozen:
    """Unit tests for frozen transfer function computation."""

    # ------------------------------------------------------------------
    # Test 1: Result fields and dimensions
    # ------------------------------------------------------------------
    def test_result_fields_and_dims(self) -> None:
        """FrozenResult has all fields with correct shapes."""
        rng = np.random.default_rng(2001)
        p, q, N, L = 2, 1, 20, 5
        A_true = np.array([[0.9, 0.1], [-0.1, 0.8]])
        B_true = np.array([[0.5], [0.3]])
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

        ltv = ltv_disc(X, U, lambda_=1e4, uncertainty=True)
        frz = ltv_disc_frozen(ltv)

        # Check all fields exist
        assert hasattr(frz, "frequency")
        assert hasattr(frz, "frequency_hz")
        assert hasattr(frz, "time_steps")
        assert hasattr(frz, "response")
        assert hasattr(frz, "response_std")
        assert hasattr(frz, "sample_time")
        assert hasattr(frz, "method")
        assert frz.method == "ltv_disc_frozen"

        # Shapes
        nf = 128
        nk = N
        assert frz.frequency.shape == (nf,)
        assert frz.response.shape == (nf, p, q, nk)
        assert frz.response_std is not None
        assert frz.response_std.shape == (nf, p, q, nk)

    # ------------------------------------------------------------------
    # Test 2: Custom frequencies and time steps
    # ------------------------------------------------------------------
    def test_custom_freqs_and_timesteps(self) -> None:
        """Custom frequency vector and time_steps produce correct shapes."""
        rng = np.random.default_rng(2001)
        p, q, N, L = 2, 1, 20, 5
        A_true = np.array([[0.9, 0.1], [-0.1, 0.8]])
        B_true = np.array([[0.5], [0.3]])
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

        ltv = ltv_disc(X, U, lambda_=1e4, uncertainty=True)

        w = np.logspace(-2, np.log10(np.pi), 50)
        time_steps = np.array([0, 9, 19])  # 0-based
        frz = ltv_disc_frozen(ltv, frequencies=w, time_steps=time_steps)

        assert frz.response.shape == (50, p, q, 3)
        np.testing.assert_allclose(frz.time_steps, time_steps)
        np.testing.assert_allclose(frz.frequency, w, atol=1e-15)

    # ------------------------------------------------------------------
    # Test 3: LTI frozen TF matches analytic
    # ------------------------------------------------------------------
    def test_lti_frozen_matches_analytic(self) -> None:
        """Noiseless LTI with high lambda -> frozen TF matches analytic G(w)."""
        rng = np.random.default_rng(2002)
        p, q, N, L = 2, 1, 30, 10
        A_true = np.array([[0.9, 0.1], [-0.1, 0.8]])
        B_true = np.array([[0.5], [0.3]])
        X = np.zeros((N + 1, p, L))
        U = rng.standard_normal((N, q, L))
        for ll in range(L):
            X[0, :, ll] = rng.standard_normal(p)
            for k in range(N):
                # Noiseless
                X[k + 1, :, ll] = A_true @ X[k, :, ll] + B_true.ravel() * U[k, :, ll]

        ltv = ltv_disc(X, U, lambda_=1e8)
        w = np.linspace(0.01, np.pi, 64)
        frz = ltv_disc_frozen(ltv, frequencies=w)

        Ip = np.eye(p)
        max_err = 0.0
        mid_k = N // 2  # 0-based mid time step
        for iw in range(len(w)):
            z = np.exp(1j * w[iw])
            G_analytic = np.linalg.solve(z * Ip - A_true, B_true)
            G_frozen_mid = frz.response[iw, :, :, mid_k]
            err = np.linalg.norm(G_frozen_mid - G_analytic) / max(
                np.linalg.norm(G_analytic), np.finfo(float).eps
            )
            max_err = max(max_err, err)

        assert max_err < 0.02, f"LTI frozen TF maxRelErr={max_err:.4f}, expected < 0.02"

    # ------------------------------------------------------------------
    # Test 4: Time-varying DC gain
    # ------------------------------------------------------------------
    def test_time_varying_dc_gain(self) -> None:
        """DC gain |G(0.01)| = |B/(1-A)| should correlate with A(k)."""
        rng = np.random.default_rng(2003)
        p, q, N, L = 1, 1, 40, 15
        A_seq = 0.5 + 0.4 * np.arange(N) / (N - 1)  # ramp 0.5 -> 0.9
        B_true = 1.0
        X = np.zeros((N + 1, p, L))
        U = rng.standard_normal((N, q, L))
        for ll in range(L):
            X[0, :, ll] = rng.standard_normal(p)
            for k in range(N):
                X[k + 1, :, ll] = (
                    A_seq[k] * X[k, :, ll] + B_true * U[k, :, ll] + 0.01 * rng.standard_normal(p)
                )

        ltv = ltv_disc(X, U, lambda_=1e2)
        w_dc = np.array([0.01])
        frz = ltv_disc_frozen(ltv, frequencies=w_dc)

        dc_gain = np.abs(frz.response[0, 0, 0, :])
        corr_mat = np.corrcoef(dc_gain, A_seq)
        assert corr_mat[0, 1] > 0.8, (
            f"DC gain should correlate with A(k) ramp (corr={corr_mat[0, 1]:.3f})"
        )

    # ------------------------------------------------------------------
    # Test 5: Uncertainty finite positive
    # ------------------------------------------------------------------
    def test_uncertainty_finite_positive(self) -> None:
        """With uncertainty=True, response_std > 0 and finite everywhere."""
        rng = np.random.default_rng(2004)
        p, q, N, L = 2, 1, 20, 5
        A_true = np.array([[0.9, 0.1], [-0.1, 0.8]])
        B_true = np.array([[0.5], [0.3]])
        X = np.zeros((N + 1, p, L))
        U = rng.standard_normal((N, q, L))
        for ll in range(L):
            X[0, :, ll] = rng.standard_normal(p)
            for k in range(N):
                X[k + 1, :, ll] = (
                    A_true @ X[k, :, ll]
                    + B_true.ravel() * U[k, :, ll]
                    + 0.05 * rng.standard_normal(p)
                )

        ltv = ltv_disc(X, U, lambda_=1e4, uncertainty=True)
        frz = ltv_disc_frozen(ltv)

        assert frz.response_std is not None
        assert np.all(frz.response_std > 0), "ResponseStd should be positive"
        assert np.all(np.isfinite(frz.response_std)), "ResponseStd should be finite"

    # ------------------------------------------------------------------
    # Test 6: No uncertainty -> response_std is None
    # ------------------------------------------------------------------
    def test_no_uncertainty_empty_std(self) -> None:
        """Without uncertainty, response_std is None."""
        rng = np.random.default_rng(2004)
        p, q, N, L = 2, 1, 20, 5
        A_true = np.array([[0.9, 0.1], [-0.1, 0.8]])
        B_true = np.array([[0.5], [0.3]])
        X = np.zeros((N + 1, p, L))
        U = rng.standard_normal((N, q, L))
        for ll in range(L):
            X[0, :, ll] = rng.standard_normal(p)
            for k in range(N):
                X[k + 1, :, ll] = (
                    A_true @ X[k, :, ll]
                    + B_true.ravel() * U[k, :, ll]
                    + 0.05 * rng.standard_normal(p)
                )

        ltv = ltv_disc(X, U, lambda_=1e4)
        frz = ltv_disc_frozen(ltv)

        assert frz.response_std is None

    # ------------------------------------------------------------------
    # Test 7: SampleTime affects frequency_hz
    # ------------------------------------------------------------------
    def test_sample_time_affects_hz(self) -> None:
        """SampleTime=0.01 -> frequency_hz = frequency / (2*pi*0.01)."""
        rng = np.random.default_rng(2004)
        p, q, N, L = 2, 1, 20, 5
        A_true = np.array([[0.9, 0.1], [-0.1, 0.8]])
        B_true = np.array([[0.5], [0.3]])
        X = np.zeros((N + 1, p, L))
        U = rng.standard_normal((N, q, L))
        for ll in range(L):
            X[0, :, ll] = rng.standard_normal(p)
            for k in range(N):
                X[k + 1, :, ll] = (
                    A_true @ X[k, :, ll]
                    + B_true.ravel() * U[k, :, ll]
                    + 0.05 * rng.standard_normal(p)
                )

        ltv = ltv_disc(X, U, lambda_=1e4)
        frz = ltv_disc_frozen(ltv, sample_time=0.01)

        expected_hz = frz.frequency / (2.0 * np.pi * 0.01)
        np.testing.assert_allclose(
            frz.frequency_hz,
            expected_hz,
            atol=1e-12,
            err_msg="frequency_hz should be w / (2*pi*Ts)",
        )

    # ------------------------------------------------------------------
    # Test 8: Exact variance vs brute-force Kronecker
    # ------------------------------------------------------------------
    def test_exact_variance_vs_kronecker(self) -> None:
        """Compact rank-1 formula matches brute-force Kronecker product."""
        rng = np.random.default_rng(2008)
        p, q, N, L = 2, 1, 15, 8
        A_true = np.array([[0.9, 0.1], [-0.1, 0.8]])
        B_true = np.array([[0.5], [0.3]])
        X = np.zeros((N + 1, p, L))
        U = rng.standard_normal((N, q, L))
        for ll in range(L):
            X[0, :, ll] = rng.standard_normal(p)
            for k in range(N):
                noise = np.array([0.02, 0.08]) * rng.standard_normal(p)
                X[k + 1, :, ll] = A_true @ X[k, :, ll] + B_true.ravel() * U[k, :, ll] + noise

        ltv = ltv_disc(X, U, lambda_=1e3, uncertainty=True, covariance_mode="full")

        w = np.array([0.5, 1.5])
        ki_vec = np.array([4, 9])  # 0-based
        frz = ltv_disc_frozen(ltv, frequencies=w, time_steps=ki_vec)

        d = p + q
        Sigma = ltv.noise_cov
        Ip = np.eye(p)
        max_rel_err = 0.0

        for ik in range(2):
            ki = ki_vec[ik]
            Pk = ltv.p_cov[:, :, ki]
            Ak = ltv.a[:, :, ki]
            Bk = ltv.b[:, :, ki]

            for iw in range(2):
                z = np.exp(1j * w[iw])
                R = np.linalg.solve(z * Ip - Ak, Ip)
                Gk = R @ Bk

                for a in range(p):
                    for b in range(q):
                        var_bf = 0.0 + 0.0j
                        for r in range(d):
                            for j in range(p):
                                if r < p:
                                    dG_rj = R[a, j] * Gk[r, b]
                                elif r == p + b:
                                    dG_rj = R[a, j]
                                else:
                                    dG_rj = 0.0
                                for s in range(d):
                                    for l2 in range(p):
                                        if s < p:
                                            dG_sl = R[a, l2] * Gk[s, b]
                                        elif s == p + b:
                                            dG_sl = R[a, l2]
                                        else:
                                            dG_sl = 0.0
                                        var_bf += np.conj(dG_rj) * dG_sl * Sigma[j, l2] * Pk[r, s]

                        var_impl = frz.response_std[iw, a, b, ik] ** 2
                        rel_err = abs(var_bf.real - var_impl) / max(var_bf.real, 1e-30)
                        max_rel_err = max(max_rel_err, rel_err)

        assert max_rel_err < 1e-10, (
            f"Exact formula should match brute-force (maxRelErr={max_rel_err:.2e})"
        )
