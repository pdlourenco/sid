# Copyright (c) 2026 Pedro Lourenco. All rights reserved.
# This code is released under the MIT License. See LICENSE file in the
# project root for full license information.
#
# This module is part of the Open Source System Identification Toolbox (SID).
# https://github.com/pdlourenco/sid-matlab

"""Tests for sid_uncertainty from sid._internal.uncertainty."""

from __future__ import annotations

import numpy as np

from sid._internal.hann_win import hann_win
from sid._internal.uncertainty import sid_uncertainty


class TestSidUncertainty:
    """Unit tests for asymptotic standard deviations."""

    def test_phi_v_std_formula(self) -> None:
        """PhiVStd = sqrt(2*CW/N) * |PhiV| for known W, CW, PhiV, N.

        W = [1, 0.5, 0], CW = 1^2 + 2*(0.5^2 + 0^2) = 1.5.
        """
        W = np.array([1.0, 0.5, 0.0])
        N = 100
        PhiV = np.array([2.0, 3.0, 5.0])
        CW = 1.5
        expected_std = np.sqrt(2 * CW / N) * np.abs(PhiV)

        _, PhiVStd = sid_uncertainty(None, PhiV, None, N, W)
        np.testing.assert_allclose(PhiVStd, expected_std, atol=1e-12)

    def test_g_std_formula(self) -> None:
        """GStd matches the SISO formula with known G, PhiV, Coh.

        GVar = (CW / N) * |G|^2 * (1 - Coh) / Coh.
        """
        G = np.array([1 + 1j, 2 - 0.5j, 0.5 + 0.3j])
        PhiV = np.array([1.0, 1.0, 1.0])
        Coh = np.array([0.9, 0.5, 0.99])
        N = 1000
        W = hann_win(10)
        CW = W[0] ** 2 + 2 * np.sum(W[1:] ** 2)

        expected_var = (CW / N) * np.abs(G) ** 2 * (1 - Coh) / Coh
        expected_std = np.sqrt(expected_var)

        GStd, _ = sid_uncertainty(G, PhiV, Coh, N, W)
        np.testing.assert_allclose(GStd, expected_std, atol=1e-12)

    def test_time_series(self) -> None:
        """Time series mode: G=None gives g_std=None, phi_v_std computed."""
        W = hann_win(5)
        PhiV = np.array([1.0, 2.0])
        GStd, PhiVStd = sid_uncertainty(None, PhiV, None, 100, W)
        assert GStd is None
        assert len(PhiVStd) == 2

    def test_mimo_no_phi_u(self) -> None:
        """MIMO with Coh=None and PhiU=None gives GStd = NaN array."""
        rng = np.random.default_rng(42)
        G_mimo = rng.standard_normal((10, 2, 3)) + 1j * rng.standard_normal((10, 2, 3))
        PhiV_mimo = np.abs(rng.standard_normal((10, 2, 2)))
        GStd, _ = sid_uncertainty(G_mimo, PhiV_mimo, None, 200, hann_win(10))
        assert GStd.shape == G_mimo.shape
        assert np.all(np.isnan(GStd))

    def test_scaling_with_n(self) -> None:
        """Larger N gives smaller standard deviations."""
        W = hann_win(20)
        G = np.array([1 + 0.5j])
        PhiV = np.array([1.0])
        Coh = np.array([0.8])

        GStd1, PhiVStd1 = sid_uncertainty(G, PhiV, Coh, 100, W)
        GStd2, PhiVStd2 = sid_uncertainty(G, PhiV, Coh, 10000, W)
        assert GStd2[0] < GStd1[0]
        assert PhiVStd2[0] < PhiVStd1[0]

    def test_coherence_effect(self) -> None:
        """Higher coherence gives lower GStd."""
        W = hann_win(20)
        G = np.array([1 + 0.5j])
        PhiV = np.array([1.0])

        GStd_hi, _ = sid_uncertainty(G, PhiV, np.array([0.99]), 1000, W)
        GStd_lo, _ = sid_uncertainty(G, PhiV, np.array([0.3]), 1000, W)
        assert GStd_hi[0] < GStd_lo[0]

    def test_cw_m2(self) -> None:
        """CW for M=2: W=[1,0.5,0], CW=1+2*0.25=1.5."""
        W = hann_win(2)
        _, PhiVStd = sid_uncertainty(None, np.array([1.0]), None, 100, W)
        expected = np.sqrt(2 * 1.5 / 100) * 1.0
        np.testing.assert_allclose(PhiVStd[0], expected, atol=1e-12)

    def test_zero_coherence_clamp(self) -> None:
        """Zero coherence is clamped to 1e-10, result stays finite."""
        G = np.array([1.0 + 0j])
        PhiV = np.array([1.0])
        Coh = np.array([0.0])
        GStd, _ = sid_uncertainty(G, PhiV, Coh, 1000, hann_win(10))
        assert np.all(np.isfinite(GStd))
