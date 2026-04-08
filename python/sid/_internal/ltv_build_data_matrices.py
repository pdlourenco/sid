# Copyright (c) 2026 Pedro Lourenco. All rights reserved.
# This code is released under the MIT License. See LICENSE file in the
# project root for full license information.
#
# This module is part of the Open Source System Identification Toolbox (SID).
# https://github.com/pdlourenco/sid-matlab

"""Construct COSMIC data matrices for uniform and variable-length trajectories."""

from __future__ import annotations

import numpy as np


def build_data_matrices(
    X: np.ndarray,
    U: np.ndarray,
    N: int,
    p: int,
    q: int,
    L: int,
) -> tuple[np.ndarray, np.ndarray]:
    """Build COSMIC data matrices for uniform-length trajectories.

    This is the Python port of ``sidLTVbuildDataMatrices.m``.

    Constructs the data matrices ``D(k)`` and next-state matrices ``Xl(k)``
    for the COSMIC algorithm. Each is normalized by ``1 / sqrt(N)``.

    .. math::

        D(k)  = [X(k)^T \\; U(k)^T] / \\sqrt{N}, \\quad
        Xl(k) = X(k+1)^T / \\sqrt{N}

    Parameters
    ----------
    X : ndarray, shape ``(N+1, p, L)``
        State data for all trajectories.
    U : ndarray, shape ``(N, q, L)``
        Input data for all trajectories.
    N : int
        Number of time steps.
    p : int
        State dimension.
    q : int
        Input dimension.
    L : int
        Number of trajectories.

    Returns
    -------
    D : ndarray, shape ``(L, p+q, N)``
        Data matrices, normalized by ``1 / sqrt(N)``.
    Xl : ndarray, shape ``(L, p, N)``
        Next-state matrices, normalized by ``1 / sqrt(N)``.

    Examples
    --------
    >>> D, Xl = build_data_matrices(X, U, N, p, q, L)  # doctest: +SKIP

    Notes
    -----
    **Specification:** SPEC.md §8.2 -- Inputs

    The ``1 / sqrt(N)`` scaling makes the normal equations independent of
    *N*.

    See Also
    --------
    build_data_matrices_var_len : Variant for variable-length trajectories.
    sid._internal.ltv_build_block_terms.build_block_terms :
        Consumes these matrices to form the block tridiagonal system.

    Changelog
    ---------
    2026-04-08 : First version by Pedro Lourenco.
    """

    sqrtN = np.sqrt(N)
    D = np.zeros((L, p + q, N))
    Xl = np.zeros((L, p, N))

    for k in range(N):
        # X[k, :, :] is (p, L) -> transpose to (L, p)
        # U[k, :, :] is (q, L) -> transpose to (L, q)
        D[:, :, k] = (
            np.hstack(
                [
                    X[k, :, :].T,
                    U[k, :, :].T,
                ]
            )
            / sqrtN
        )
        Xl[:, :, k] = X[k + 1, :, :].T / sqrtN

    return D, Xl


def build_data_matrices_var_len(
    X: list[np.ndarray],
    U: list[np.ndarray],
    N: int,
    p: int,
    q: int,
    L: int,
    horizons: np.ndarray,
) -> tuple[list[np.ndarray], list[np.ndarray]]:
    """Build COSMIC data matrices for variable-length trajectories.

    This is the Python port of ``sidLTVbuildDataMatricesVarLen.m``.

    Like :func:`build_data_matrices`, but handles trajectories with
    different lengths.  At each time step *k*, only trajectories with
    ``horizon > k`` contribute.  Returns lists instead of 3-D arrays.

    Parameters
    ----------
    X : list of ndarray
        Length-*L* list where ``X[i]`` has shape ``(N_i+1, p)``.
    U : list of ndarray
        Length-*L* list where ``U[i]`` has shape ``(N_i, q)``.
    N : int
        Maximum horizon across all trajectories.
    p : int
        State dimension.
    q : int
        Input dimension.
    L : int
        Number of trajectories.
    horizons : ndarray, shape ``(L,)``
        Horizon length of each trajectory.

    Returns
    -------
    D : list of ndarray
        Length-*N* list where ``D[k]`` has shape ``(L_k, p+q)``,
        normalized by ``1 / sqrt(N)``.
    Xl : list of ndarray
        Length-*N* list where ``Xl[k]`` has shape ``(L_k, p)``,
        normalized by ``1 / sqrt(N)``.

    Examples
    --------
    >>> D, Xl = build_data_matrices_var_len(  # doctest: +SKIP
    ...     X, U, N, p, q, L, horizons
    ... )

    Notes
    -----
    **Specification:** SPEC.md §8.8 -- Variable-Length Trajectories

    The normalization uses the *maximum* horizon ``N``, matching the
    uniform-trajectory convention so that the effective regularization
    strength ``lambda`` is consistent regardless of how many trajectories
    are active at each time step.

    See Also
    --------
    build_data_matrices : Variant for uniform-length trajectories.
    sid._internal.ltv_build_block_terms.build_block_terms :
        Consumes these matrices to form the block tridiagonal system.

    Changelog
    ---------
    2026-04-08 : First version by Pedro Lourenco.
    """

    D: list[np.ndarray] = []
    Xl: list[np.ndarray] = []
    sqrtN = np.sqrt(N)

    for k in range(N):
        active = [i for i in range(L) if horizons[i] > k]
        Lk = len(active)

        if Lk == 0:
            D.append(np.zeros((0, p + q)))
            Xl.append(np.zeros((0, p)))
            continue

        Dk = np.zeros((Lk, p + q))
        Xlk = np.zeros((Lk, p))

        for ii, traj_idx in enumerate(active):
            Dk[ii, :] = (
                np.hstack(
                    [
                        X[traj_idx][k, :],
                        U[traj_idx][k, :],
                    ]
                )
                / sqrtN
            )
            Xlk[ii, :] = X[traj_idx][k + 1, :] / sqrtN

        D.append(Dk)
        Xl.append(Xlk)

    return D, Xl
