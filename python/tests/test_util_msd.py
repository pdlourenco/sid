# Copyright (c) 2026 Pedro Lourenco. All rights reserved.
# This code is released under the MIT License. See LICENSE file in the
# project root for full license information.
#
# This module is part of the Open Source System Identification Toolbox (SID).
# https://github.com/pdlourenco/sid

"""Unit tests for the example-suite SMD plant helpers.

The module under test (:mod:`util_msd`) lives in ``python/examples/`` as a
sibling of the notebooks. It is made importable here by the test-suite
``conftest.py``, which prepends the examples directory to ``sys.path``.

The n=3 parity check against a legacy ``sid._internal.test_msd`` fixture
used to live in this module; that fixture has been removed as part of
the consolidation described in the commit history, and its regression
role is now carried by ``test_cross_validation.TestCrossValidationTestMSD``
which compares the helper's output against the MATLAB-generated
``testdata/reference_test_msd.json`` vectors.
"""

from __future__ import annotations

import numpy as np
import pytest
from scipy.linalg import expm

from util_msd import util_msd, util_msd_ltv, util_msd_nl  # noqa: E402


# ----------------------------------------------------------------------
# LTI chain — analytical and shape checks
# ----------------------------------------------------------------------


class TestUtilMsdLti:
    """LTI chain discretization."""

    def test_n1_sdof_against_analytical(self) -> None:
        """n=1 SDOF: matches the direct ZOH of the analytic 2x2 Ac."""
        m = np.array([1.0])
        k = np.array([100.0])
        c = np.array([0.5])
        F = np.array([[1.0]])
        Ts = 0.01

        Ad, Bd = util_msd(m, k, c, F, Ts)

        Ac = np.array([[0.0, 1.0], [-k[0] / m[0], -c[0] / m[0]]])
        Bc = np.array([[0.0], [1.0 / m[0]]])
        Ad_expected = expm(Ac * Ts)
        Bd_expected = np.linalg.solve(Ac, (Ad_expected - np.eye(2))) @ Bc

        np.testing.assert_allclose(Ad, Ad_expected, atol=1e-14)
        np.testing.assert_allclose(Bd, Bd_expected, atol=1e-14)

    def test_n2_eigenvalues_match_modal_analysis(self) -> None:
        """n=2 continuous eigenvalues should equal the analytic modes."""
        m = np.array([1.0, 1.0])
        k = np.array([100.0, 80.0])
        c = np.array([0.0, 0.0])  # undamped so eigenvalues are pure imaginary
        F = np.eye(2)
        Ts = 0.001  # small so discrete eigenvalues approach continuous

        Ad, _ = util_msd(m, k, c, F, Ts)
        # Continuous eigenvalues of Ac: solve det(K - omega^2 M) = 0.
        # For M=I and K = [[k0+k1, -k1],[-k1, k1]] with k0=100, k1=80:
        # eigenvalues of K: trace/2 +- sqrt((trace/2)^2 - det)
        K = np.array([[180.0, -80.0], [-80.0, 80.0]])
        omega_sq = np.sort(np.linalg.eigvalsh(K))
        omega = np.sqrt(omega_sq)

        # Discrete Ad eigenvalues should have magnitude 1 (undamped) and
        # arguments ~ +-omega_i * Ts.
        disc_eigs = np.linalg.eigvals(Ad)
        mags = np.abs(disc_eigs)
        args = np.sort(np.unique(np.abs(np.angle(disc_eigs))))

        np.testing.assert_allclose(mags, np.ones(4), atol=1e-6)
        np.testing.assert_allclose(args, omega * Ts, atol=1e-6)

    def test_n5_shapes_and_stability(self) -> None:
        """n=5 chain: shapes correct and all discrete eigenvalues inside unit circle."""
        n = 5
        m = np.ones(n)
        k = np.linspace(100.0, 60.0, n)
        c = 0.5 * np.ones(n)
        F = np.eye(n)[:, :1]  # single force on mass 0
        Ts = 0.01

        Ad, Bd = util_msd(m, k, c, F, Ts)

        assert Ad.shape == (2 * n, 2 * n)
        assert Bd.shape == (2 * n, 1)

        eigs = np.linalg.eigvals(Ad)
        assert np.all(np.abs(eigs) < 1.0), "damped chain must be strictly stable"

    def test_rejects_mismatched_shapes(self) -> None:
        m = np.array([1.0, 1.0])
        k = np.array([100.0])  # wrong length
        c = np.array([0.5, 0.5])
        F = np.eye(2)
        with pytest.raises(ValueError):
            util_msd(m, k, c, F, 0.01)


# ----------------------------------------------------------------------
# LTV chain — construction and collapse-to-LTI
# ----------------------------------------------------------------------


class TestUtilMsdLtv:
    """LTV chain construction from per-step parameter arrays."""

    def test_ltv_collapses_to_lti_when_all_constant(self) -> None:
        """All 1-D inputs + explicit N returns replicated LTI result."""
        m = np.array([1.0, 1.0])
        k = np.array([100.0, 80.0])
        c = np.array([0.5, 0.5])
        F = np.array([[1.0], [0.0]])
        Ts = 0.01
        N = 50

        Ad, Bd = util_msd_ltv(m, k, c, F, Ts, N=N)
        Ad_ref, Bd_ref = util_msd(m, k, c, F, Ts)

        assert Ad.shape == (4, 4, N)
        assert Bd.shape == (4, 1, N)
        for kk in range(N):
            np.testing.assert_allclose(Ad[:, :, kk], Ad_ref, atol=1e-14)
            np.testing.assert_allclose(Bd[:, :, kk], Bd_ref, atol=1e-14)

    def test_ltv_constant_inputs_without_N_raises(self) -> None:
        m = np.array([1.0])
        k = np.array([100.0])
        c = np.array([0.5])
        F = np.array([[1.0]])
        with pytest.raises(ValueError, match="time-invariant"):
            util_msd_ltv(m, k, c, F, 0.01)

    def test_ltv_tracks_time_varying_stiffness(self) -> None:
        """2-mass chain with k1 ramping: first-mode frequency should drift."""
        N = 100
        m = np.array([1.0, 1.0])
        c = np.array([0.5, 0.5])
        k_tv = np.zeros((2, N))
        k_tv[0, :] = np.linspace(100.0, 400.0, N)
        k_tv[1, :] = 80.0
        F = np.array([[1.0], [0.0]])
        Ts = 0.01

        Ad, _ = util_msd_ltv(m, k_tv, c, F, Ts)
        assert Ad.shape == (4, 4, N)

        # First and last slices should match a direct LTI call.
        Ad0_ref, _ = util_msd(m, k_tv[:, 0], c, F, Ts)
        AdN_ref, _ = util_msd(m, k_tv[:, -1], c, F, Ts)
        np.testing.assert_allclose(Ad[:, :, 0], Ad0_ref, atol=1e-14)
        np.testing.assert_allclose(Ad[:, :, -1], AdN_ref, atol=1e-14)

        # Middle slice should differ from both endpoints.
        assert not np.allclose(Ad[:, :, N // 2], Ad0_ref, atol=1e-6)
        assert not np.allclose(Ad[:, :, N // 2], AdN_ref, atol=1e-6)

    def test_ltv_step_change_discontinuity(self) -> None:
        """Discrete step change in k[0] at N/2 should be visible in Ad."""
        N = 40
        m = np.array([1.0])
        c = np.array([0.5])
        k_step = np.zeros((1, N))
        k_step[0, : N // 2] = 100.0
        k_step[0, N // 2 :] = 50.0
        F = np.array([[1.0]])
        Ts = 0.01

        Ad, _ = util_msd_ltv(m, k_step, c, F, Ts)
        # Steps on either side of the discontinuity should be constant.
        np.testing.assert_allclose(Ad[:, :, 0], Ad[:, :, N // 2 - 1], atol=1e-14)
        np.testing.assert_allclose(Ad[:, :, N // 2], Ad[:, :, N - 1], atol=1e-14)
        # Across the discontinuity they should differ.
        assert not np.allclose(Ad[:, :, N // 2 - 1], Ad[:, :, N // 2], atol=1e-6)

    def test_ltv_time_varying_F(self) -> None:
        """Time-varying input distribution is accepted."""
        N = 20
        m = np.array([1.0, 1.0])
        k = np.array([100.0, 80.0])
        c = np.array([0.5, 0.5])
        F_tv = np.zeros((2, 1, N))
        F_tv[:, 0, :] = np.linspace(0, 1, N)[None, :]  # fades in on both masses
        Ts = 0.01
        Ad, Bd = util_msd_ltv(m, k, c, F_tv, Ts)
        assert Ad.shape == (4, 4, N)
        assert Bd.shape == (4, 1, N)
        # Ad is LTI (m, k, c are fixed) — all slices equal.
        for kk in range(1, N):
            np.testing.assert_allclose(Ad[:, :, kk], Ad[:, :, 0], atol=1e-14)
        # Bd scales with F.
        ratio = Bd[3, 0, -1] / Bd[3, 0, 10] if Bd[3, 0, 10] != 0 else 0
        assert np.isfinite(ratio)


# ----------------------------------------------------------------------
# Nonlinear (Duffing) simulation
# ----------------------------------------------------------------------


class TestUtilMsdNl:
    """RK4 simulation of the chain with cubic stiffness."""

    def test_linear_case_matches_zoh_ad_bd(self) -> None:
        """k_cubic=0: RK4 matches the exact ZOH Ad/Bd propagation to O(Ts^5)."""
        m = np.array([1.0])
        k = np.array([100.0])
        kc = np.array([0.0])
        c = np.array([0.5])
        F = np.array([[1.0]])
        Ts = 0.005  # small enough that RK4 local error is ~1e-12
        N = 200
        rng = np.random.default_rng(7)
        u = rng.standard_normal((N, 1))

        # RK4 reference
        x_rk = util_msd_nl(m, k, kc, c, F, Ts, u, substeps=4)

        # Exact ZOH reference via Ad/Bd propagation
        Ad, Bd = util_msd(m, k, c, F, Ts)
        x_zoh = np.zeros((N + 1, 2))
        for kk in range(N):
            x_zoh[kk + 1] = Ad @ x_zoh[kk] + (Bd @ u[kk])

        np.testing.assert_allclose(x_rk, x_zoh, atol=1e-8, rtol=1e-6)

    def test_duffing_deviates_from_linear_at_large_amplitude(self) -> None:
        """Cubic term should produce measurable deviation at high drive."""
        m = np.array([1.0])
        k_lin = np.array([100.0])
        k_cub = np.array([1000.0])
        c = np.array([0.5])
        F = np.array([[1.0]])
        Ts = 0.01
        N = 1000
        rng = np.random.default_rng(42)
        u_big = 3.0 * rng.standard_normal((N, 1))  # drive to large amplitude

        x_lin = util_msd_nl(m, k_lin, np.zeros(1), c, F, Ts, u_big, substeps=4)
        x_duf = util_msd_nl(m, k_lin, k_cub, c, F, Ts, u_big, substeps=4)

        # The two trajectories must remain finite and differ materially.
        assert np.all(np.isfinite(x_lin))
        assert np.all(np.isfinite(x_duf))

        rms_diff = float(np.sqrt(np.mean((x_lin[:, 0] - x_duf[:, 0]) ** 2)))
        rms_lin = float(np.sqrt(np.mean(x_lin[:, 0] ** 2)))
        assert rms_diff / rms_lin > 0.05, "Duffing response should differ > 5 %"

    def test_zero_input_zero_initial_state(self) -> None:
        """No input, no initial displacement -> trajectory stays at origin."""
        m = np.array([1.0, 1.0])
        k = np.array([100.0, 80.0])
        kc = np.array([500.0, 0.0])
        c = np.array([0.5, 0.5])
        F = np.array([[1.0], [0.0]])
        Ts = 0.01
        N = 200
        u = np.zeros((N, 1))

        x = util_msd_nl(m, k, kc, c, F, Ts, u, substeps=1)
        np.testing.assert_allclose(x, 0.0, atol=1e-14)

    def test_stability_at_target_parameters(self) -> None:
        """Canonical notebook parameters must not blow up at default substeps."""
        m = np.array([1.0])
        k_lin = np.array([100.0])
        k_cub = np.array([1000.0])
        c = np.array([0.5])
        F = np.array([[1.0]])
        Ts = 0.01
        N = 2000
        rng = np.random.default_rng(0)
        u = 2.0 * rng.standard_normal((N, 1))

        x = util_msd_nl(m, k_lin, k_cub, c, F, Ts, u, substeps=4)
        assert np.all(np.isfinite(x))
        assert np.max(np.abs(x[:, 0])) < 10.0

    def test_rejects_mismatched_input_width(self) -> None:
        m = np.array([1.0])
        k = np.array([100.0])
        kc = np.array([0.0])
        c = np.array([0.5])
        F = np.array([[1.0]])  # q=1
        u = np.zeros((10, 2))  # q=2 mismatch
        with pytest.raises(ValueError, match="inputs"):
            util_msd_nl(m, k, kc, c, F, 0.01, u)
