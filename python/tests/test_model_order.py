# Copyright (c) 2026 Pedro Lourenco. All rights reserved.
# This code is released under the MIT License. See LICENSE file in the
# project root for full license information.
#
# This module is part of the Open Source System Identification Toolbox (SID).
# https://github.com/pdlourenco/sid

"""Tests for model_order from sid.

Port of test_sidModelOrder.m (8 tests -- simplified from 14).
"""

from __future__ import annotations

import numpy as np
import pytest
from scipy.signal import lfilter

from sid._exceptions import SidError
from sid.freq_bt import freq_bt
from sid.model_order import model_order


class TestModelOrder:
    """Unit tests for model order estimation via Hankel SVD."""

    # ------------------------------------------------------------------
    # Test 1: Known 2nd-order system
    # ------------------------------------------------------------------
    def test_known_second_order(self) -> None:
        """2nd-order system -> model_order returns n=2."""
        rng = np.random.default_rng(42)
        N = 1000
        u = rng.standard_normal(N)

        # 2nd-order system: poles at 0.9*exp(+/-j*pi/4)
        r_pole = 0.9
        theta = np.pi / 4
        a1 = -2 * r_pole * np.cos(theta)
        a2 = r_pole**2
        # y(k) + a1*y(k-1) + a2*y(k-2) = u(k-1)
        # Transfer function: z^{-1} / (1 + a1*z^{-1} + a2*z^{-2})
        b_poly = np.array([0.0, 1.0])
        a_poly = np.array([1.0, a1, a2])
        y = lfilter(b_poly, a_poly, u) + 0.01 * rng.standard_normal(N)

        G = freq_bt(y, u, window_size=60)
        n, sv = model_order(G)

        assert n == 2, f"Expected n=2 for 2nd-order system, got n={n}"

    # ------------------------------------------------------------------
    # Test 2: Known 4th-order system
    # ------------------------------------------------------------------
    def test_known_fourth_order(self) -> None:
        """4th-order system (two pairs of poles) -> n=4."""
        rng = np.random.default_rng(100)
        N = 2000
        u = rng.standard_normal(N)

        # Two resonant modes:
        #   Pair 1: 0.85*exp(+/-j*0.4)
        #   Pair 2: 0.85*exp(+/-j*2.0)
        r_pole = 0.85
        a1 = np.array([1.0, -2 * r_pole * np.cos(0.4), r_pole**2])
        a2 = np.array([1.0, -2 * r_pole * np.cos(2.0), r_pole**2])
        a_poly = np.convolve(a1, a2)
        y = lfilter([1.0], a_poly, u) + 0.005 * rng.standard_normal(N)

        G = freq_bt(y, u, window_size=100)
        n, sv = model_order(G)

        assert n == 4, f"Expected n=4 for 4th-order system, got n={n}"

    # ------------------------------------------------------------------
    # Test 3: Threshold method
    # ------------------------------------------------------------------
    def test_threshold_method(self) -> None:
        """Threshold method gives n in a reasonable range."""
        rng = np.random.default_rng(42)
        N = 1000
        u = rng.standard_normal(N)

        r_pole = 0.9
        theta = np.pi / 4
        a1 = -2 * r_pole * np.cos(theta)
        a2 = r_pole**2
        b_poly = np.array([0.0, 1.0])
        a_poly = np.array([1.0, a1, a2])
        y = lfilter(b_poly, a_poly, u) + 0.01 * rng.standard_normal(N)

        G = freq_bt(y, u, window_size=60)
        n, sv = model_order(G, threshold=0.01)

        assert n >= 1, f"Threshold method: n should be >= 1, got {n}"
        assert n <= 10, f"Threshold method: n should be <= 10, got {n}"

    # ------------------------------------------------------------------
    # Test 4: Output fields
    # ------------------------------------------------------------------
    def test_output_fields(self) -> None:
        """sv dict has 'singular_values' and 'horizon'."""
        rng = np.random.default_rng(200)
        N = 500
        u = rng.standard_normal(N)
        y = lfilter([1.0], [1.0, -0.9], u) + 0.01 * rng.standard_normal(N)

        G = freq_bt(y, u)
        n, sv = model_order(G)

        assert "singular_values" in sv, "Missing key: singular_values"
        assert "horizon" in sv, "Missing key: horizon"
        assert isinstance(sv["singular_values"], np.ndarray)
        assert sv["singular_values"].ndim == 1
        assert isinstance(sv["horizon"], (int, np.integer))
        assert np.all(sv["singular_values"] >= 0), "Singular values must be non-negative"

    # ------------------------------------------------------------------
    # Test 5: Custom horizon
    # ------------------------------------------------------------------
    def test_custom_horizon(self) -> None:
        """horizon=20 -> sv['horizon']==20."""
        rng = np.random.default_rng(300)
        N = 500
        u = rng.standard_normal(N)
        y = lfilter([0.5], [1.0, -0.8], u) + 0.01 * rng.standard_normal(N)

        G = freq_bt(y, u)
        n, sv = model_order(G, horizon=20)

        assert sv["horizon"] == 20, f"Expected horizon=20, got {sv['horizon']}"
        assert len(sv["singular_values"]) > 0, "SingularValues should not be empty"

    # ------------------------------------------------------------------
    # Test 6: SISO first-order system
    # ------------------------------------------------------------------
    def test_siso_first_order(self) -> None:
        """AR(1) a=0.9, N=2000 -> n=1."""
        rng = np.random.default_rng(400)
        N = 2000
        u = rng.standard_normal(N)
        y = lfilter([0.8], [1.0, -0.9], u) + 0.005 * rng.standard_normal(N)

        G = freq_bt(y, u, window_size=40)
        n, sv = model_order(G)

        assert n == 1, f"Expected n=1 for 1st-order system, got n={n}"

    # ------------------------------------------------------------------
    # Test 7: Singular values are decreasing
    # ------------------------------------------------------------------
    def test_singular_values_decreasing(self) -> None:
        """sigmas[0] >= sigmas[1] >= ..."""
        rng = np.random.default_rng(500)
        N = 1000
        u = rng.standard_normal(N)
        y = lfilter([1.0], [1.0, -0.85], u) + 0.01 * rng.standard_normal(N)

        G = freq_bt(y, u)
        n, sv = model_order(G)

        sigmas = sv["singular_values"]
        for i in range(len(sigmas) - 1):
            assert sigmas[i] >= sigmas[i + 1] - 1e-15, (
                f"Singular values not non-increasing at index {i}: {sigmas[i]} < {sigmas[i + 1]}"
            )

    # ------------------------------------------------------------------
    # Test 8: Error on bad input
    # ------------------------------------------------------------------
    def test_error_bad_input(self) -> None:
        """Non-FreqResult -> raises SidError."""
        with pytest.raises(SidError) as exc_info:
            model_order(42)
        assert exc_info.value.code == "bad_input"


class TestModelOrderStrictlyProper:
    """Regression for the lag-0 Hankel bug (issue #139).

    Uses exact analytic frequency responses on the (0, pi] grid so the order
    is unambiguous. The pre-fix code kept the lag-0 IFFT sample and built the
    Hankel of z^{-1}G(z), overestimating the order by one for a strictly
    proper plant with a nonzero finite zero.
    """

    @staticmethod
    def _result(G: np.ndarray, nf: int):
        from types import SimpleNamespace

        w = np.arange(1, nf + 1) * np.pi / nf
        return SimpleNamespace(frequency=w, response=G)

    def test_finite_zero_second_order(self) -> None:
        """G(z) = (z + 0.5) / (z^2 - 1.2728 z + 0.81), true n = 2.

        The lag-0 fix is verified directly on the Hankel singular values:
        pre-fix the third value was a structural ~0.4/6.1 (order 3); post-fix
        it collapses to numerical noise (~1e-4). Both the threshold and the
        (default) gap method return n=2. Pre-ADR-0004 the MATLAB gap method
        included the resolvable->floor cliff and returned n=3 here while Python
        returned n=2; the count-based rule makes both agree. See ADR-0004.
        """
        nf = 128
        w = np.arange(1, nf + 1) * np.pi / nf
        z = np.exp(1j * w)
        G = ((z + 0.5) / (z**2 - 1.2728 * z + 0.81)).reshape(nf, 1, 1)
        n, sv = model_order(self._result(G, nf), threshold=1e-2)
        assert n == 2, f"finite-zero strictly proper plant: expected n=2, got {n}"
        # The third singular value must be numerical noise, not structural.
        s = sv["singular_values"]
        assert s[2] / s[0] < 1e-3, f"spurious structural SV: {s[2] / s[0]:.2e}"
        # Gap method (default) now also returns n=2 (count-based, cliff excluded).
        n_gap, _ = model_order(self._result(G, nf))
        assert n_gap == 2, f"finite-zero gap method: expected n=2, got {n_gap}"

    def test_zero_at_origin_unchanged(self) -> None:
        """Masked family (numerator z) still returns n = 2."""
        nf = 128
        w = np.arange(1, nf + 1) * np.pi / nf
        z = np.exp(1j * w)
        G = (z / (z**2 - 1.2728 * z + 0.81)).reshape(nf, 1, 1)
        n, _ = model_order(self._result(G, nf), threshold=1e-2)
        assert n == 2

    def test_mimo_block_hankel(self) -> None:
        """2x2 MIMO map with four distinct first-order channels -> n = 4.

        Uses the threshold method: the four modes span two orders of
        magnitude in Hankel singular value, so the largest-gap heuristic
        latches onto the first (largest) gap; the threshold counts all four
        structural values. Exercises the MIMO block-Hankel assembly.
        """
        nf = 256
        w = np.arange(1, nf + 1) * np.pi / nf
        z = np.exp(1j * w)
        G = np.zeros((nf, 2, 2), dtype=complex)
        G[:, 0, 0] = 1.0 / (z - 0.7)
        G[:, 0, 1] = 0.3 / (z - 0.4)
        G[:, 1, 0] = 0.2 / (z - 0.6)
        G[:, 1, 1] = 1.0 / (z - 0.5)
        n, sv = model_order(self._result(G, nf), horizon=8, threshold=1e-3)
        assert n == 4, f"MIMO 4-pole plant: expected n=4, got {n}"
        # Block-Hankel is (r*ny) x (r*nu) = 16 x 16.
        assert len(sv["singular_values"]) == 16
