# Copyright (c) 2026 Pedro Lourenco. All rights reserved.
# This code is released under the MIT License. See LICENSE file in the
# project root for full license information.
#
# This module is part of the Open Source System Identification Toolbox (SID).
# https://github.com/pdlourenco/sid

"""ZOH discretization of a 3-mass spring-damper chain."""

from __future__ import annotations

import numpy as np
from scipy.linalg import expm


def test_msd(
    m: np.ndarray,
    k_spring: np.ndarray,
    c_damp: np.ndarray,
    F: np.ndarray,
    Ts: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Discretize a 3-mass spring-damper chain using exact ZOH.

    This is the Python port of ``sidTestMSD.m``.

    Builds the continuous-time state-space model for a chain of three
    masses connected by springs and dampers (wall--k1--m1--k2--m2--k3--m3)
    and discretizes it using the matrix exponential (zero-order hold).

    State vector: ``[x1, x2, x3, v1, v2, v3]`` (positions, velocities).

    Parameters
    ----------
    m : ndarray, shape ``(3,)``
        Masses (kg).
    k_spring : ndarray, shape ``(3,)``
        Spring constants (N/m).
    c_damp : ndarray, shape ``(3,)``
        Damping coefficients (N s/m).
    F : ndarray, shape ``(3, q)``
        Force input distribution matrix.
    Ts : float
        Sample time (s).

    Returns
    -------
    Ad : ndarray, shape ``(6, 6)``
        Discrete dynamics matrix.
    Bd : ndarray, shape ``(6, q)``
        Discrete input matrix.

    Examples
    --------
    >>> m = np.array([1.0, 1.0, 1.0])
    >>> k = np.array([10.0, 10.0, 10.0])
    >>> c = np.array([0.1, 0.1, 0.1])
    >>> F = np.eye(3)
    >>> Ad, Bd = test_msd(m, k, c, F, 0.01)  # doctest: +SKIP

    Notes
    -----
    **Specification:** (Test utility -- not yet in SPEC.md)

    See Also
    --------
    sid.ltv_disc : LTV system identification that uses this test model.

    Changelog
    ---------
    2026-04-08 : First version by Pedro Lourenco.
    """

    m = np.asarray(m, dtype=np.float64)
    k_spring = np.asarray(k_spring, dtype=np.float64)
    c_damp = np.asarray(c_damp, dtype=np.float64)
    F = np.asarray(F, dtype=np.float64)

    M_mat = np.diag(m)

    K_mat = np.array(
        [
            [k_spring[0] + k_spring[1], -k_spring[1], 0.0],
            [-k_spring[1], k_spring[1] + k_spring[2], -k_spring[2]],
            [0.0, -k_spring[2], k_spring[2]],
        ]
    )

    C_mat = np.array(
        [
            [c_damp[0] + c_damp[1], -c_damp[1], 0.0],
            [-c_damp[1], c_damp[1] + c_damp[2], -c_damp[2]],
            [0.0, -c_damp[2], c_damp[2]],
        ]
    )

    n = len(m)
    Minv = np.linalg.solve(M_mat, np.eye(n))

    Ac = np.block(
        [
            [np.zeros((n, n)), np.eye(n)],
            [-Minv @ K_mat, -Minv @ C_mat],
        ]
    )

    Bc = np.block(
        [
            [np.zeros((n, F.shape[1]))],
            [Minv @ F],
        ]
    )

    ns = 2 * n
    Ad = expm(Ac * Ts)
    Bd = np.linalg.solve(Ac, (Ad - np.eye(ns))) @ Bc

    return Ad, Bd
