# Copyright (c) 2026 Pedro Lourenco. All rights reserved.
# This code is released under the MIT License. See LICENSE file in the
# project root for full license information.
#
# This module is part of the Open Source System Identification Toolbox (SID).
# https://github.com/pdlourenco/sid

"""Bayesian posterior covariance diagonal blocks for COSMIC."""

from __future__ import annotations

import numpy as np


def uncertainty_backward_pass(
    S_scaled: np.ndarray,
    lambda_: np.ndarray,
    N: int,
    d: int,
) -> np.ndarray:
    """Compute diagonal blocks of the inverse COSMIC Hessian.

    This is the Python port of ``sidLTVuncertaintyBackwardPass.m``.

    Computes ``P(k) = [A_unscaled^{-1}]_{kk}``, the diagonal blocks of
    the inverse of the unscaled COSMIC Hessian.  These are the row
    covariance matrices needed for Bayesian uncertainty estimation.

    The COSMIC algorithm normalizes data by ``1 / sqrt(N)``, so
    *S_scaled* contains ``D_s^T D_s + regularization``.  This function
    reconstructs the unscaled diagonal blocks, then computes left and
    right Schur complements to obtain the diagonal blocks of the inverse:

    .. math::

        P(k) = (\\Lambda_k^L + \\Lambda_k^R - S_{kk})^{-1}

    Parameters
    ----------
    S_scaled : ndarray, shape ``(d, d, N)``
        Scaled block diagonal terms (from :func:`build_block_terms`,
        with ``1 / sqrt(N)`` normalization).
    lambda_ : ndarray, shape ``(N-1,)``
        Regularization weights.
    N : int
        Number of time steps.
    d : int
        Combined dimension, ``d = p + q``.

    Returns
    -------
    P : ndarray, shape ``(d, d, N)``
        Diagonal blocks of the inverse Hessian.
        ``Cov(vec(C(k))) = Sigma (x) P(k)``.

    Examples
    --------
    >>> P = uncertainty_backward_pass(S, lambda_, N, d)  # doctest: +SKIP

    Notes
    -----
    **Specification:** SPEC.md §8.9 -- Bayesian Uncertainty Estimation

    **Algorithm:**

    1. Forward pass: left Schur complements ``Lbd^L(k)`` of the scaled Hessian.
    2. Backward pass: right Schur complements ``Lbd^R(k)``.
    3. Combine: ``P_scaled(k) = (Lbd^L(k) + Lbd^R(k) - S_scaled(k))^{-1}``.
    4. Scale to the returned estimator: ``P(k) = P_scaled(k) / N`` (effective
       prior weight ``N * lambda``; see SPEC.md §8.9.2 and issue #137).

    Complexity: O(N * d^3).

    **References:**

    Carvalho, Soares, Lourenco, Ventura. "COSMIC: fast closed-form
    identification from large-scale data for LTV systems."
    arXiv:2112.04355, 2022.

    See Also
    --------
    sid._internal.ltv_cosmic_solve.cosmic_solve :
        Solves the COSMIC system whose Hessian is inverted here.
    sid._internal.ltv_build_block_terms.build_block_terms :
        Produces the *S_scaled* input.

    Changelog
    ---------
    2026-04-08 : First version by Pedro Lourenco.
    """

    eye_d = np.eye(d)

    # ---- Schur recursion on the SCALED Hessian blocks --------------------
    # SPEC.md §8.9.2: the returned MAP minimises ||unscaled residual||^2 +
    # N*lambda*||dC||^2 (effective prior weight N*lambda, from the 1/sqrt(N)
    # scaling of §8.3.2), so its Hessian is A_est = N * A_scaled and the
    # reported posterior is P(k) = [A_est^{-1}]_kk = [A_scaled^{-1}]_kk / N =
    # P_scaled(k) / N.  We therefore run the left/right Schur complements on the
    # scaled blocks S_scaled (with the un-inflated lambda_k couplings) to obtain
    # P_scaled, then divide by N -- exactly equivalent to inflating S -> N*S and
    # the couplings -> (N*lambda_k)^2, but minimal.  (The previous code rebuilt
    # V^T V + lambda F^T F -- unscaled data with an un-inflated prior -- the
    # posterior of a different estimator, overstating P by up to a factor N; see
    # issue #137.)
    S = S_scaled

    # ---- Left Schur complements -- forward pass ---------------------------
    LbdL = np.zeros((d, d, N))
    LbdL[:, :, 0] = S[:, :, 0]
    for k in range(1, N):
        LbdL[:, :, k] = S[:, :, k] - lambda_[k - 1] ** 2 * np.linalg.solve(LbdL[:, :, k - 1], eye_d)

    # ---- Right Schur complements -- backward pass -------------------------
    LbdR = np.zeros((d, d, N))
    LbdR[:, :, N - 1] = S[:, :, N - 1]
    for k in range(N - 2, -1, -1):
        LbdR[:, :, k] = S[:, :, k] - lambda_[k] ** 2 * np.linalg.solve(LbdR[:, :, k + 1], eye_d)

    # ---- Combine: P_scaled(k) = (LbdL + LbdR - S)^{-1}, then P = P_scaled / N
    P = np.zeros((d, d, N))
    for k in range(N):
        M = LbdL[:, :, k] + LbdR[:, :, k] - S[:, :, k]
        P[:, :, k] = np.linalg.solve(M, eye_d) / N

    return P
