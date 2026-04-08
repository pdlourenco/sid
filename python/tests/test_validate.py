# Copyright (c) 2026 Pedro Lourenco. All rights reserved.
# This code is released under the MIT License. See LICENSE file in the
# project root for full license information.
#
# This module is part of the Open Source System Identification Toolbox (SID).
# https://github.com/pdlourenco/sid-matlab

"""Tests for validate_data from sid._internal.validate_data."""

from __future__ import annotations

import numpy as np
import pytest

from sid._exceptions import SidError
from sid._internal.validate_data import validate_data


class TestValidateData:
    """Unit tests for input validation and orientation."""

    def test_basic_siso(self) -> None:
        """Basic SISO: y=1000x1, u=None returns correct metadata."""
        rng = np.random.default_rng(42)
        y = rng.standard_normal(1000)
        y_out, u_out, N, ny, nu, is_ts, n_traj = validate_data(y, None)
        assert N == 1000
        assert ny == 1
        assert nu == 0
        assert is_ts is True
        assert n_traj == 1

    def test_row_to_column(self) -> None:
        """1-D vector y gets reshaped to (N, 1) column."""
        rng = np.random.default_rng(42)
        y_row = rng.standard_normal(30)  # 1-D array (Python row equivalent)
        y_out, _ = validate_data(y_row, None)[:2]
        assert y_out.shape == (30, 1)

    def test_error_too_short(self) -> None:
        """N < 2 raises SidError with code 'too_short'."""
        with pytest.raises(SidError) as exc:
            validate_data(np.array([1.0]), None)
        assert exc.value.code == "too_short"

    def test_error_complex_data(self) -> None:
        """Complex y raises SidError with code 'complex_data'."""
        with pytest.raises(SidError) as exc:
            validate_data(np.array([1 + 1j, 2, 3]), None)
        assert exc.value.code == "complex_data"

    def test_error_nan(self) -> None:
        """NaN in y raises SidError with code 'non_finite'."""
        with pytest.raises(SidError) as exc:
            validate_data(np.array([1.0, np.nan, 3.0]), None)
        assert exc.value.code == "non_finite"

    def test_error_inf(self) -> None:
        """Inf in y raises SidError with code 'non_finite'."""
        with pytest.raises(SidError) as exc:
            validate_data(np.array([1.0, np.inf, 3.0]), None)
        assert exc.value.code == "non_finite"

    def test_error_size_mismatch(self) -> None:
        """y(100) and u(50) raises SidError with code 'size_mismatch'."""
        rng = np.random.default_rng(42)
        with pytest.raises(SidError) as exc:
            validate_data(rng.standard_normal(100), rng.standard_normal(50))
        assert exc.value.code == "size_mismatch"

    def test_multi_output(self) -> None:
        """y(50x2), u(50x3) returns ny=2, nu=3."""
        rng = np.random.default_rng(42)
        y = rng.standard_normal((50, 2))
        u = rng.standard_normal((50, 3))
        _, _, N, ny, nu, is_ts, _ = validate_data(y, u)
        assert N == 50
        assert ny == 2
        assert nu == 3
        assert is_ts is False

    def test_time_series(self) -> None:
        """u=None means time series: is_ts=True, nu=0."""
        rng = np.random.default_rng(42)
        y = rng.standard_normal(50)
        _, _, _, _, nu, is_ts, _ = validate_data(y, None)
        assert nu == 0
        assert is_ts is True
