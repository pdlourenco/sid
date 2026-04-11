# Copyright (c) 2026 Pedro Lourenco. All rights reserved.
# This code is released under the MIT License. See LICENSE file in the
# project root for full license information.
#
# This module is part of the Open Source System Identification Toolbox (SID).
# https://github.com/pdlourenco/sid

"""Biased cross-covariance estimate for spectral analysis."""

from __future__ import annotations

import numpy as np


def sid_cov(x: np.ndarray, z: np.ndarray, max_lag: int) -> np.ndarray:
    """Biased cross-covariance estimate for lags ``0 .. max_lag``.

    This is the Python port of ``sidCov.m``.

    Computes the biased sample cross-covariance:

    .. math::

        \\hat R_{xz}(\\tau) = \\frac{1}{N}
            \\sum_{t=\\tau+1}^{N} x(t)\\,z(t-\\tau)^T

    for :math:`\\tau = 0, 1, \\dots, \\text{max\\_lag}`.

    For multi-trajectory data (3-D arrays with shape ``(N, p, L)``),
    the covariance is ensemble-averaged across *L* trajectories.

    Parameters
    ----------
    x : ndarray, shape ``(N, p)`` or ``(N, p, L)``
        First signal.
    z : ndarray, shape ``(N, q)`` or ``(N, q, L)``
        Second signal (may equal *x* for auto-covariance).
    max_lag : int
        Maximum lag *M* (non-negative integer).

    Returns
    -------
    R : ndarray
        Covariance estimates.  Shape ``(max_lag + 1, p, q)`` for matrix
        signals, or ``(max_lag + 1,)`` for scalar signals (``p == q == 1``).

    Examples
    --------
    >>> R = sid_cov(x, x, 30)  # doctest: +SKIP

    Notes
    -----
    **Specification:** SPEC.md §2.3 -- Covariance Estimation

    See Also
    --------
    sid._internal.hann_win.hann_win : Hann window applied to covariances.
    sid._internal.windowed_dft.windowed_dft : DFT of windowed covariances.

    Changelog
    ---------
    2026-04-08 : First version by Pedro Lourenco.
    """

    N: int = x.shape[0]
    p: int = x.shape[1]
    q: int = z.shape[1]

    # Detect multi-trajectory (3-D arrays)
    if x.ndim == 3:
        L: int = x.shape[2]
    else:
        L = 1

    R = np.zeros((max_lag + 1, p, q), dtype=np.float64)

    if L == 1:
        # Single trajectory
        for tau in range(max_lag + 1):
            # R(tau) = (1/N) * x[tau:N, :]^T @ z[0:N-tau, :]
            R[tau, :, :] = x[tau:N, :].T @ z[: N - tau, :] / N
    else:
        # Ensemble average across L trajectories
        for tau in range(max_lag + 1):
            R_sum = np.zeros((p, q), dtype=np.float64)
            for traj in range(L):
                R_sum += x[tau:N, :, traj].T @ z[: N - tau, :, traj]
            R[tau, :, :] = R_sum / (L * N)

    # Squeeze trailing singleton dimensions for scalar signals
    if p == 1 and q == 1:
        R = R.ravel()

    return R
