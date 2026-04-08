# Copyright (c) 2026 Pedro Lourenco. All rights reserved.
# This code is released under the MIT License. See LICENSE file in the
# project root for full license information.
#
# This module is part of the Open Source System Identification Toolbox (SID).
# https://github.com/pdlourenco/sid-matlab

"""Estimate noise covariance from COSMIC residuals."""

from __future__ import annotations

import numpy as np


def estimate_noise_cov(
    C: np.ndarray,
    D: np.ndarray | list[np.ndarray],
    Xl: np.ndarray | list[np.ndarray],
    P: np.ndarray,
    cov_mode: str,
    N: int,
    p: int,
    q: int,
) -> tuple[np.ndarray, float]:
    """Estimate noise covariance from COSMIC residuals.

    This is the Python port of ``sidEstimateNoiseCov.m``.

    Estimates the ``(p, p)`` noise covariance matrix from the residuals
    of a COSMIC solution.  The data *D* and *Xl* are scaled by
    ``1 / sqrt(N)`` (COSMIC convention).  The scaled residuals have noise
    covariance ``Sigma / N``.  This function returns the *unscaled*
    noise covariance ``Sigma``.

    Parameters
    ----------
    C : ndarray, shape ``(d, p, N)``
        COSMIC solution matrices.
    D : ndarray, shape ``(L, d, N)`` or list of ndarray
        Data matrices.  Either a 3-D array (uniform trajectories) or a
        length-*N* list where ``D[k]`` has shape ``(L_k, d)``.
    Xl : ndarray, shape ``(L, p, N)`` or list of ndarray
        Next-state matrices.  Either a 3-D array or a length-*N* list
        where ``Xl[k]`` has shape ``(L_k, p)``.
    P : ndarray, shape ``(d, d, N)``
        Posterior covariance diagonal blocks from
        :func:`uncertainty_backward_pass`.
    cov_mode : str
        Covariance structure: ``'diagonal'``, ``'full'``, or
        ``'isotropic'``.
    N : int
        Number of time steps.
    p : int
        State dimension.
    q : int
        Input dimension.

    Returns
    -------
    Sigma : ndarray, shape ``(p, p)``
        Estimated noise covariance matrix.
    dof : float
        Effective degrees of freedom used.

    Examples
    --------
    >>> Sigma, dof = estimate_noise_cov(  # doctest: +SKIP
    ...     C, D, Xl, P, 'diagonal', N, p, q
    ... )

    Notes
    -----
    **Specification:** SPEC.md §8.9.3 -- Noise Covariance Estimation

    The effective degrees of freedom are computed using the unscaled
    hat-matrix trace:  ``dof = total_obs - N * trace(sum_k DtD(k) P(k))``.
    A conservative fallback is used if *dof* is non-positive.

    See Also
    --------
    sid._internal.ltv_uncertainty_backward_pass.uncertainty_backward_pass :
        Produces the *P* input.
    sid._internal.extract_std.extract_std :
        Uses *Sigma* to compute element-wise standard deviations.

    Changelog
    ---------
    2026-04-08 : First version by Pedro Lourenco.
    """

    d = p + q
    use_cell = isinstance(D, list)

    # ---- Accumulate scaled residual scatter matrix and count observations --
    SSR_scaled = np.zeros((p, p))
    total_obs = 0

    for k in range(N):
        Ck = C[:, :, k]
        Dk = D[k] if use_cell else D[:, :, k]
        Xlk = Xl[k] if use_cell else Xl[:, :, k]
        Lk = Dk.shape[0]

        if Lk == 0:
            continue

        Ek = Xlk - Dk @ Ck
        SSR_scaled += Ek.T @ Ek
        total_obs += Lk

    # ---- Effective degrees of freedom -------------------------------------
    trace_sum = 0.0
    for k in range(N):
        Dk = D[k] if use_cell else D[:, :, k]
        if Dk.shape[0] > 0:
            DtD = Dk.T @ Dk
            trace_sum += np.sum(DtD * P[:, :, k])  # trace of product

    dof = total_obs - N * trace_sum

    # Conservative fallback if exact dof is non-positive
    if dof <= 0:
        dof = total_obs - N * d
        if dof <= 0:
            dof = max(total_obs, 1)

    # Unscaled noise covariance
    Sigma = N * SSR_scaled / dof

    # Apply covariance mode restriction
    if cov_mode == "diagonal":
        Sigma = np.diag(np.diag(Sigma))
    elif cov_mode == "isotropic":
        Sigma = (np.trace(Sigma) / p) * np.eye(p)

    return Sigma, float(dof)
