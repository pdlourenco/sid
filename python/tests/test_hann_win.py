# Copyright (c) 2026 Pedro Lourenco. All rights reserved.
# This code is released under the MIT License. See LICENSE file in the
# project root for full license information.
#
# This module is part of the Open Source System Identification Toolbox (SID).
# https://github.com/pdlourenco/sid-matlab

"""Tests for hann_win from sid._internal.hann_win."""

from __future__ import annotations

import numpy as np
import pytest

from sid._internal.hann_win import hann_win


class TestHannWin:
    """Unit tests for the Hann lag window function."""

    @pytest.mark.parametrize("M", [2, 5, 10, 30, 100])
    def test_boundary_values(self, M: int) -> None:
        """W[0]==1 and W[M]==0 for various M."""
        W = hann_win(M)
        np.testing.assert_allclose(W[0], 1.0, atol=1e-15)
        np.testing.assert_allclose(W[M], 0.0, atol=1e-15)

    @pytest.mark.parametrize("M", [2, 7, 50])
    def test_output_shape(self, M: int) -> None:
        """Output shape is (M+1,)."""
        W = hann_win(M)
        assert W.shape == (M + 1,)

    def test_value_range(self) -> None:
        """All window values are in [0, 1]."""
        W = hann_win(50)
        assert np.all(W >= 0)
        assert np.all(W <= 1)

    def test_known_m2(self) -> None:
        """M=2 gives [1, 0.5, 0]."""
        W = hann_win(2)
        expected = np.array([1.0, 0.5 * (1 + np.cos(np.pi / 2)), 0.5 * (1 + np.cos(np.pi))])
        np.testing.assert_allclose(W, expected, atol=1e-15)

    def test_known_m4(self) -> None:
        """M=4 gives 5 values from the formula."""
        W = hann_win(4)
        tau = np.arange(5, dtype=np.float64)
        expected = 0.5 * (1 + np.cos(np.pi * tau / 4))
        np.testing.assert_allclose(W, expected, atol=1e-15)

    def test_monotonicity(self) -> None:
        """Window is monotonically decreasing."""
        W = hann_win(20)
        assert np.all(np.diff(W) <= 0)

    def test_midpoint_even(self) -> None:
        """For M=10, W[5]==0.5."""
        M = 10
        W = hann_win(M)
        np.testing.assert_allclose(W[5], 0.5, atol=1e-15)
