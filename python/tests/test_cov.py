# Copyright (c) 2026 Pedro Lourenco. All rights reserved.
# This code is released under the MIT License. See LICENSE file in the
# project root for full license information.
#
# This module is part of the Open Source System Identification Toolbox (SID).
# https://github.com/pdlourenco/sid-matlab

"""Tests for sid_cov from sid._internal.cov."""

from __future__ import annotations

import numpy as np
import pytest

from sid._internal.cov import sid_cov


class TestSidCov:
    """Unit tests for biased cross-covariance estimation."""

    def test_hand_computed_auto(self) -> None:
        """Known hand-computed auto-covariance values.

        x = [1, 2, 3, 4], N = 4.
        R(0) = (1+4+9+16)/4 = 7.5
        R(1) = (2*1+3*2+4*3)/4 = 20/4 = 5.0
        R(2) = (3*1+4*2)/4 = 11/4 = 2.75
        """
        x = np.array([[1.0], [2.0], [3.0], [4.0]])
        R = sid_cov(x, x, 2)
        np.testing.assert_allclose(R[0], 7.5, atol=1e-12)
        np.testing.assert_allclose(R[1], 5.0, atol=1e-12)
        np.testing.assert_allclose(R[2], 2.75, atol=1e-12)

    def test_hand_computed_cross(self) -> None:
        """Known hand-computed cross-covariance values.

        x = [1, 2, 3, 4], z = [4, 3, 2, 1], N = 4.
        R_xz(0) = (1*4+2*3+3*2+4*1)/4 = 20/4 = 5.0
        R_xz(1) = (2*4+3*3+4*2)/4 = 25/4 = 6.25
        R_xz(2) = (3*4+4*3)/4 = 24/4 = 6.0
        """
        x = np.array([[1.0], [2.0], [3.0], [4.0]])
        z = np.array([[4.0], [3.0], [2.0], [1.0]])
        R = sid_cov(x, z, 2)
        np.testing.assert_allclose(R[0], 5.0, atol=1e-12)
        np.testing.assert_allclose(R[1], 6.25, atol=1e-12)
        np.testing.assert_allclose(R[2], 6.0, atol=1e-12)

    def test_output_shape_scalar(self) -> None:
        """Scalar signals: x(10x1), z(10x1), maxlag=4 gives shape (5,)."""
        rng = np.random.default_rng(42)
        x = rng.standard_normal((10, 1))
        z = rng.standard_normal((10, 1))
        R = sid_cov(x, z, 4)
        assert R.shape == (5,)

    def test_output_shape_matrix(self) -> None:
        """Matrix signals: x(10x2), z(10x3), maxlag=4 gives shape (5,2,3)."""
        rng = np.random.default_rng(42)
        x = rng.standard_normal((10, 2))
        z = rng.standard_normal((10, 3))
        R = sid_cov(x, z, 4)
        assert R.shape == (5, 2, 3)

    def test_biased_estimator(self) -> None:
        """Biased estimator divides by N: x=[1,0,0,0], R(0)=1/4=0.25."""
        x = np.array([[1.0], [0.0], [0.0], [0.0]])
        R = sid_cov(x, x, 3)
        np.testing.assert_allclose(R[0], 0.25, atol=1e-12)

    def test_white_noise(self) -> None:
        """White noise: R(0) approx 1, R(tau>0) approx 0."""
        rng = np.random.default_rng(123)
        N = 100000
        x = rng.standard_normal((N, 1))
        R = sid_cov(x, x, 10)
        np.testing.assert_allclose(R[0], 1.0, atol=0.02)
        for tau in range(1, 11):
            np.testing.assert_allclose(R[tau], 0.0, atol=0.02)

    def test_lag0_biased_variance(self) -> None:
        """R(0) = (1/N)*sum(x^2) for data with non-zero mean."""
        rng = np.random.default_rng(99)
        x = rng.standard_normal((200, 1)) * 3 + 5
        R = sid_cov(x, x, 0)
        biased_var = np.sum(x**2) / len(x)
        np.testing.assert_allclose(R[0], biased_var, atol=1e-10)

    def test_single_lag(self) -> None:
        """maxlag=0 returns shape (1,) for scalar output."""
        x = np.array([[1.0], [2.0], [3.0]])
        R = sid_cov(x, x, 0)
        assert R.shape == (1,)
