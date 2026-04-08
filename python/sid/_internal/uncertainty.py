# Copyright (c) 2026 Pedro Lourenco. All rights reserved.
# This code is released under the MIT License. See LICENSE file in the
# project root for full license information.
#
# This module is part of the Open Source System Identification Toolbox (SID).
# https://github.com/pdlourenco/sid-matlab

"""Asymptotic standard deviations for spectral estimates."""

from __future__ import annotations

import numpy as np


def sid_uncertainty(
    G: np.ndarray | None,
    phi_v: np.ndarray,
    coh: np.ndarray | None,
    N: int,
    W: np.ndarray,
    n_traj: int = 1,
    phi_u: np.ndarray | None = None,
) -> tuple[np.ndarray | None, np.ndarray]:
    """Asymptotic standard deviations for spectral estimates.

    This is the Python port of ``sidUncertainty.m``.

    Computes the asymptotic standard deviations of the frequency response
    and noise spectrum estimates, based on Ljung (1999), pp. 184 and 188.

    For multi-trajectory ensemble averaging, pass ``n_traj > 1``.  The
    variance scales as ``1 / (N * n_traj)``.

    For MIMO systems, pass *phi_u* (the input auto-spectrum) to enable a
    diagonal approximation of the per-element variance.

    Parameters
    ----------
    G : ndarray or None
        Complex frequency response estimate, shape ``(nf,)``,
        ``(nf, ny, nu)``, or ``None`` for time-series mode.
    phi_v : ndarray
        Noise spectrum estimate, shape ``(nf,)`` or ``(nf, ny, ny)``.
    coh : ndarray or None
        Squared coherence, shape ``(nf,)``.  ``None`` for time-series or
        MIMO.
    N : int
        Number of data samples per trajectory.
    W : ndarray, shape ``(M+1,)``
        Hann window values for lags ``0 .. M``.
    n_traj : int, optional
        Number of trajectories (default 1).
    phi_u : ndarray or None, optional
        Input auto-spectrum, shape ``(nf, nu, nu)``.  Required for MIMO
        uncertainty; ignored for SISO.

    Returns
    -------
    g_std : ndarray or None
        Standard deviation of *G*, same shape as *G*, or ``None`` when
        *G* is ``None``.
    phi_v_std : ndarray
        Standard deviation of *phi_v*, same shape as *phi_v*.

    Examples
    --------
    >>> g_std, phi_v_std = sid_uncertainty(G, phi_v, coh, N, W)  # doctest: +SKIP

    Notes
    -----
    **Specification:** SPEC.md §3 -- Uncertainty Estimation

    See Also
    --------
    sid.freq_bt : Main function that calls this for uncertainty estimation.
    sid._internal.hann_win.hann_win : Window used in the C_W computation.

    Changelog
    ---------
    2026-04-08 : First version by Pedro Lourenco.
    """

    # Effective sample size
    Neff: int = N * n_traj

    # Window norm: C_W = W[0]^2 + 2 * sum(W[1:]^2)
    CW: float = float(W[0] ** 2 + 2.0 * np.sum(W[1:] ** 2))

    # ---- Noise spectrum variance ----------------------------------------
    phi_v_std = np.sqrt(2.0 * CW / Neff) * np.abs(phi_v)

    # ---- Transfer function variance -------------------------------------
    if G is None:
        return None, phi_v_std

    if coh is not None:
        # SISO case
        eps_floor = 1e-10
        coh_safe = np.maximum(coh, eps_floor)
        g_var = (CW / Neff) * np.abs(G) ** 2 * (1.0 - coh_safe) / coh_safe
        g_std = np.sqrt(g_var)

    elif phi_u is not None:
        # MIMO case: diagonal approximation
        nf: int = G.shape[0]
        ny: int = G.shape[1]
        nu: int = G.shape[2]
        g_std = np.zeros((nf, ny, nu), dtype=np.float64)
        eps_floor = 1e-10

        for k in range(nf):
            for ii in range(ny):
                # Diagonal of noise spectrum at this frequency
                if phi_v.ndim == 3:
                    phiV_ii = float(np.real(phi_v[k, ii, ii]))
                else:
                    phiV_ii = float(np.real(phi_v[k]))

                for jj in range(nu):
                    phiU_jj = float(np.real(phi_u[k, jj, jj]))
                    if phiU_jj > eps_floor:
                        g_std[k, ii, jj] = np.sqrt(CW / Neff * phiV_ii / phiU_jj)
                    else:
                        g_std[k, ii, jj] = np.inf
    else:
        # MIMO case without phi_u: cannot compute uncertainty
        g_std = np.full_like(G, np.nan, dtype=np.float64)

    return g_std, phi_v_std
