# Copyright (c) 2026 Pedro Lourenco. All rights reserved.
# This code is released under the MIT License. See LICENSE file in the
# project root for full license information.
#
# This module is part of the Open Source System Identification Toolbox (SID).
# https://github.com/pdlourenco/sid-matlab

"""Evaluate the COSMIC cost function."""

from __future__ import annotations

import numpy as np


def evaluate_cost(
    A: np.ndarray,
    B: np.ndarray,
    D: np.ndarray | list[np.ndarray],
    Xl: np.ndarray | list[np.ndarray],
    lambda_: np.ndarray,
    N: int,
    p: int,
    q: int,
) -> tuple[float, float, float, float]:
    """Compute the COSMIC cost function value.

    This is the Python port of ``sidLTVevaluateCost.m``.

    Evaluates the total COSMIC cost as the sum of data fidelity and
    temporal regularization terms.  Handles both 3-D array *D* (uniform
    trajectories) and list *D* (variable-length trajectories).

    Parameters
    ----------
    A : ndarray, shape ``(p, p, N)``
        Time-varying dynamics matrices.
    B : ndarray, shape ``(p, q, N)``
        Time-varying input matrices.
    D : ndarray, shape ``(L, d, N)`` or list of ndarray
        Data matrices.
    Xl : ndarray, shape ``(L, p, N)`` or list of ndarray
        Next-state matrices.
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
    cost : float
        Total cost ``fidelity + reg``.
    fidelity : float
        ``(1/2) sum_k ||D(k) C(k) - Xl(k)||_F^2``.
    reg : float
        ``(1/2) sum_k lambda_k ||C(k) - C(k+1)||_F^2`` (weighted).
    variation : float
        ``(1/2) sum_k ||C(k) - C(k+1)||_F^2`` (unweighted, for L-curve).

    Examples
    --------
    >>> cost, fid, reg, var = evaluate_cost(  # doctest: +SKIP
    ...     A, B, D, Xl, lambda_, N, p, q
    ... )

    Notes
    -----
    **Specification:** SPEC.md §8.3 -- COSMIC Algorithm

    ``C(k) = [A(k)^T; B(k)^T]`` is the stacked coefficient matrix at
    each time step.

    See Also
    --------
    sid._internal.ltv_cosmic_solve.cosmic_solve :
        Produces the *A* and *B* matrices.

    Changelog
    ---------
    2026-04-08 : First version by Pedro Lourenco.
    """

    fidelity = 0.0
    prior_vec = np.zeros(N - 1)
    use_cell = isinstance(D, list)

    for k in range(N):
        # C(k) = [A(k)'; B(k)']  shape (d, p)
        Ck = np.vstack([A[:, :, k].T, B[:, :, k].T])

        Dk = D[k] if use_cell else D[:, :, k]
        Xlk = Xl[k] if use_cell else Xl[:, :, k]

        residual = Dk @ Ck - Xlk
        fidelity += np.sum(residual**2)

        if k < N - 1:
            Ck1 = np.vstack([A[:, :, k + 1].T, B[:, :, k + 1].T])
            prior_vec[k] = np.sum((Ck - Ck1) ** 2)

    fidelity *= 0.5
    reg = 0.5 * lambda_.dot(prior_vec)
    variation = 0.5 * np.sum(prior_vec)
    cost = fidelity + reg

    return cost, fidelity, reg, variation
