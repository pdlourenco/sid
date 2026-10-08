# Copyright (c) 2026 Pedro Lourenco. All rights reserved.
# This code is released under the MIT License. See LICENSE file in the
# project root for full license information.
#
# This module is part of the Open Source System Identification Toolbox (SID).
# https://github.com/pdlourenco/sid

"""Shape regressions: one trajectory given as 3-D, and MISO ``freq_map``.

A ``(N, n, 1)`` array is one trajectory (SPEC.md §1) and must give exactly
the 2-D result; MATLAB cannot represent the trailing singleton at all, so this
is the Python-only half of the contract. ``freq_map`` on MISO data (one output,
several inputs) must run under both algorithms (§6.1, §6.8). Before the fix,
``freq_etfe``, Welch ``freq_map`` and ``detrend`` raised on 3-D ``L = 1``
input, and ``freq_map`` raised on any MISO data (review-v2 analysis §5.1 F1,
F2; §5.3 U5).
"""

from __future__ import annotations

import numpy as np
import pytest

import sid


@pytest.fixture
def siso() -> tuple[np.ndarray, np.ndarray]:
    rng = np.random.default_rng(7)
    u = rng.standard_normal(600)
    y = np.convolve(u, [0.0, 0.8, 0.3])[:600] + 0.05 * rng.standard_normal(600)
    return y, u


def _as3d(x: np.ndarray) -> np.ndarray:
    return x.reshape(x.shape[0], -1, 1)


@pytest.mark.parametrize(
    "call",
    [
        lambda y, u: sid.freq_bt(y, u),
        lambda y, u: sid.freq_btfdr(y, u),
        lambda y, u: sid.freq_etfe(y, u),
        lambda y, u: sid.freq_map(y, u, segment_length=256),
        lambda y, u: sid.freq_map(y, u, segment_length=256, algorithm="welch"),
    ],
    ids=["freq_bt", "freq_btfdr", "freq_etfe", "freq_map_bt", "freq_map_welch"],
)
def test_single_trajectory_3d_equals_2d(siso, call) -> None:
    y, u = siso
    r2 = call(y, u)
    r3 = call(_as3d(y), _as3d(u))
    np.testing.assert_array_equal(r3.response, r2.response)
    np.testing.assert_array_equal(r3.noise_spectrum, r2.noise_spectrum)
    assert r3.response.shape == r2.response.shape


def test_single_trajectory_3d_time_series_equals_2d(siso) -> None:
    y, _ = siso
    r2 = sid.freq_etfe(y, None)
    r3 = sid.freq_etfe(_as3d(y), None)
    np.testing.assert_array_equal(r3.noise_spectrum, r2.noise_spectrum)


def test_detrend_single_trajectory_3d_keeps_shape(siso) -> None:
    y, _ = siso
    x2 = np.column_stack([y, 2.0 * y + np.arange(y.size)])
    d2, t2 = sid.detrend(x2)
    d3, t3 = sid.detrend(x2[:, :, np.newaxis])
    assert d3.shape == x2.shape + (1,) and t3.shape == x2.shape + (1,)
    np.testing.assert_array_equal(d3[:, :, 0], d2)
    np.testing.assert_array_equal(t3[:, :, 0], t2)


@pytest.mark.parametrize("algorithm", ["bt", "welch"])
def test_freq_map_miso(algorithm: str) -> None:
    rng = np.random.default_rng(11)
    N, nu = 1200, 2
    u = rng.standard_normal((N, nu))
    y = 0.7 * u[:, 0] - 0.4 * u[:, 1] + 0.05 * rng.standard_normal(N)
    res = sid.freq_map(y, u, segment_length=600, algorithm=algorithm)
    nf, K = res.response.shape[:2]
    assert res.response.shape == (nf, K, 1, nu)
    assert res.noise_spectrum.shape == (nf, K, 1, 1)
    # A static MISO gain: the estimate recovers it per input channel.
    np.testing.assert_allclose(res.response[:, :, 0, 0].real.mean(), 0.7, atol=0.05)
    np.testing.assert_allclose(res.response[:, :, 0, 1].real.mean(), -0.4, atol=0.05)
