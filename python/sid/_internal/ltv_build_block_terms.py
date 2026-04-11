# Copyright (c) 2026 Pedro Lourenco. All rights reserved.
# This code is released under the MIT License. See LICENSE file in the
# project root for full license information.
#
# This module is part of the Open Source System Identification Toolbox (SID).
# https://github.com/pdlourenco/sid

"""Compute block diagonal and right-hand side terms for COSMIC."""

from __future__ import annotations

import numpy as np


def build_block_terms(
    D: np.ndarray | list[np.ndarray],
    Xl: np.ndarray | list[np.ndarray],
    lambda_: np.ndarray,
    N: int,
    p: int,
    q: int,
) -> tuple[np.ndarray, np.ndarray]:
    """Build block diagonal terms and right-hand side for the COSMIC system.

    This is the Python port of ``sidLTVbuildBlockTerms.m``.

    Computes ``S(k) = D(k)^T D(k) + reg(k) I`` and
    ``T(k) = D(k)^T Xl(k)`` for the COSMIC block tridiagonal system.
    Handles both 3-D array *D* (uniform trajectories) and list *D*
    (variable-length trajectories).

    Parameters
    ----------
    D : ndarray, shape ``(L, d, N)`` or list of ndarray
        Data matrices.  Either a 3-D array (uniform trajectories) or a
        length-*N* list where ``D[k]`` has shape ``(L_k, d)``.
    Xl : ndarray, shape ``(L, p, N)`` or list of ndarray
        Next-state matrices.  Either a 3-D array or a length-*N* list
        where ``Xl[k]`` has shape ``(L_k, p)``.
    lambda_ : ndarray, shape ``(N-1,)``
        Regularization weights.
    N : int
        Number of time steps.
    p : int
        State dimension.
    q : int
        Input dimension.

    Returns
    -------
    S : ndarray, shape ``(d, d, N)``
        Block diagonal terms with regularization added,
        where ``d = p + q``.
    T : ndarray, shape ``(d, p, N)``
        Right-hand side terms.

    Examples
    --------
    >>> S, T = build_block_terms(D, Xl, lambda_, N, p, q)  # doctest: +SKIP

    Notes
    -----
    **Specification:** SPEC.md §8.3 -- COSMIC Algorithm

    The regularization is distributed as follows:

    - ``S[:, :, 0]  += lambda_[0] * I``
    - ``S[:, :, -1] += lambda_[-1] * I``
    - ``S[:, :, k]  += (lambda_[k-1] + lambda_[k]) * I`` for interior *k*

    See Also
    --------
    sid._internal.ltv_cosmic_solve.cosmic_solve :
        Solves the block tridiagonal system built by this function.
    sid._internal.ltv_build_data_matrices.build_data_matrices :
        Produces the *D* and *Xl* inputs.

    Changelog
    ---------
    2026-04-08 : First version by Pedro Lourenco.
    """

    d = p + q
    S = np.zeros((d, d, N))
    T = np.zeros((d, p, N))
    use_cell = isinstance(D, list)

    for k in range(N):
        Dk = D[k] if use_cell else D[:, :, k]
        Xlk = Xl[k] if use_cell else Xl[:, :, k]
        S[:, :, k] = Dk.T @ Dk
        T[:, :, k] = Dk.T @ Xlk

    # Add regularization to diagonal blocks
    eye_d = np.eye(d)
    S[:, :, 0] += lambda_[0] * eye_d
    S[:, :, N - 1] += lambda_[N - 2] * eye_d
    for k in range(1, N - 1):
        S[:, :, k] += (lambda_[k - 1] + lambda_[k]) * eye_d

    return S, T
