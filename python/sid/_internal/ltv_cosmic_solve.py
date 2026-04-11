# Copyright (c) 2026 Pedro Lourenco. All rights reserved.
# This code is released under the MIT License. See LICENSE file in the
# project root for full license information.
#
# This module is part of the Open Source System Identification Toolbox (SID).
# https://github.com/pdlourenco/sid

"""COSMIC forward-backward block tridiagonal solver."""

from __future__ import annotations

import warnings

import numpy as np


def cosmic_solve(
    S: np.ndarray,
    T: np.ndarray,
    lambda_: np.ndarray,
    N: int,
    p: int,
    q: int,
) -> tuple[np.ndarray, np.ndarray]:
    """Solve the COSMIC block tridiagonal system.

    This is the Python port of ``sidLTVcosmicSolve.m``.

    Solves the block tridiagonal system arising from the regularized
    least-squares formulation in the COSMIC algorithm.  Returns the
    forward Schur complements ``Lbd`` for reuse in uncertainty
    computation.

    Parameters
    ----------
    S : ndarray, shape ``(d, d, N)``
        Block diagonal terms from :func:`build_block_terms`.
    T : ndarray, shape ``(d, p, N)``
        Right-hand side terms.
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
    C : ndarray, shape ``(d, p, N)``
        Solution ``[A(k)^T; B(k)^T]`` for each time step,
        where ``d = p + q``.
    Lbd : ndarray, shape ``(d, d, N)``
        Forward Schur complements.

    Examples
    --------
    >>> C, Lbd = cosmic_solve(S, T, lambda_, N, p, q)  # doctest: +SKIP

    Notes
    -----
    **Specification:** SPEC.md §8.3 -- COSMIC Algorithm

    **Algorithm:**

    Forward pass::

        Lbd(k) = S(k) - lambda(k-1)^2 * Lbd(k-1)^{-1}
        Y(k) = Lbd(k)^{-1} * (T(k) + lambda(k-1) * Y(k-1))

    Backward pass::

        C(k) = Y(k) + lambda(k) * Lbd(k)^{-1} * C(k+1)

    Complexity: O(N * d^3), d = p + q.

    **References:**

    Carvalho, Soares, Lourenco, Ventura. "COSMIC: fast closed-form
    identification from large-scale data for LTV systems."
    arXiv:2112.04355, 2022.

    See Also
    --------
    sid._internal.ltv_build_block_terms.build_block_terms :
        Builds the *S* and *T* inputs.
    sid._internal.ltv_uncertainty_backward_pass.uncertainty_backward_pass :
        Uses *Lbd* for posterior covariance estimation.

    Changelog
    ---------
    2026-04-08 : First version by Pedro Lourenco.
    """

    d = p + q
    Lbd = np.zeros((d, d, N))
    Y = np.zeros((d, p, N))
    C = np.zeros((d, p, N))
    eye_d = np.eye(d)

    # ---- Forward pass: Schur complement recursion -------------------------
    Lbd[:, :, 0] = S[:, :, 0]
    Y[:, :, 0] = np.linalg.solve(Lbd[:, :, 0], T[:, :, 0])

    eps_mach = np.finfo(np.float64).eps
    for k in range(1, N):
        Lbd_prev = Lbd[:, :, k - 1]
        # Reciprocal condition number estimate (matches MATLAB sidLTVcosmicSolve)
        rc = 1.0 / np.linalg.cond(Lbd_prev)
        if rc < eps_mach:
            warnings.warn(
                f"COSMIC forward pass: Lbd({k - 1}) is near-singular "
                f"(rcond={rc:.2e}). Results may be unreliable. "
                "Try adjusting lambda.",
                stacklevel=2,
            )
        Lbd[:, :, k] = S[:, :, k] - lambda_[k - 1] ** 2 * np.linalg.solve(Lbd_prev, eye_d)
        Y[:, :, k] = np.linalg.solve(
            Lbd[:, :, k],
            T[:, :, k] + lambda_[k - 1] * Y[:, :, k - 1],
        )

    # ---- Backward pass ----------------------------------------------------
    C[:, :, N - 1] = Y[:, :, N - 1]

    for k in range(N - 2, -1, -1):
        C[:, :, k] = Y[:, :, k] + lambda_[k] * np.linalg.solve(Lbd[:, :, k], C[:, :, k + 1])

    return C, Lbd
