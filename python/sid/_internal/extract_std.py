# Copyright (c) 2026 Pedro Lourenco. All rights reserved.
# This code is released under the MIT License. See LICENSE file in the
# project root for full license information.
#
# This module is part of the Open Source System Identification Toolbox (SID).
# https://github.com/pdlourenco/sid

"""Extract element-wise standard deviations of A(k) and B(k)."""

from __future__ import annotations

import numpy as np


def extract_std(
    P: np.ndarray,
    Sigma: np.ndarray,
    N: int,
    p: int,
    q: int,
) -> tuple[np.ndarray, np.ndarray]:
    """Standard deviations of A(k) and B(k) entries.

    This is the Python port of ``sidExtractStd.m``.

    Computes element-wise standard deviations of the time-varying system
    matrices from the Kronecker posterior structure
    ``Cov(vec(C(k))) = Sigma (x) P(k)``.

    Parameters
    ----------
    P : ndarray, shape ``(d, d, N)``
        Posterior covariance diagonal blocks, ``d = p + q``.
    Sigma : ndarray, shape ``(p, p)``
        Noise covariance matrix.
    N : int
        Number of time steps.
    p : int
        State dimension.
    q : int
        Input dimension.

    Returns
    -------
    AStd : ndarray, shape ``(p, p, N)``
        Standard deviations of ``A(k)`` entries.
    BStd : ndarray, shape ``(p, q, N)``
        Standard deviations of ``B(k)`` entries.

    Examples
    --------
    >>> AStd, BStd = extract_std(P, Sigma, N, p, q)  # doctest: +SKIP

    Notes
    -----
    **Specification:** SPEC.md §8.9.4 -- Standard Deviations

    Since ``C(k) = [A(k)^T; B(k)^T]``, row *a* of ``C(k)`` corresponds
    to column *a* of ``A(k)`` for ``a = 0, ..., p-1`` and column
    ``a - p`` of ``B(k)`` for ``a = p, ..., d-1``.  The variance of each
    entry is:

    - ``Var(A(k)_{b,a}) = Sigma_{bb} * P(k)_{aa}``
    - ``Var(B(k)_{b,a}) = Sigma_{bb} * P(k)_{p+a, p+a}``

    See Also
    --------
    sid._internal.estimate_noise_cov.estimate_noise_cov :
        Produces the *Sigma* input.
    sid._internal.ltv_uncertainty_backward_pass.uncertainty_backward_pass :
        Produces the *P* input.

    Changelog
    ---------
    2026-04-08 : First version by Pedro Lourenco.
    """

    AStd = np.zeros((p, p, N))
    BStd = np.zeros((p, q, N))
    sig_diag = np.diag(Sigma)  # (p,)

    for k in range(N):
        p_diag = np.diag(P[:, :, k])  # (d,)

        for a in range(p):
            AStd[:, a, k] = np.sqrt(sig_diag * p_diag[a])

        for a in range(q):
            BStd[:, a, k] = np.sqrt(sig_diag * p_diag[p + a])

    return AStd, BStd
