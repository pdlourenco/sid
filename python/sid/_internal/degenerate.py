# Copyright (c) 2026 Pedro Lourenco. All rights reserved.
# This code is released under the MIT License. See LICENSE file in the
# project root for full license information.
#
# This module is part of the Open Source System Identification Toolbox (SID).
# https://github.com/pdlourenco/sid

"""Shared degenerate-input handling for the frequency-domain estimators.

Factored out of ``freq_bt`` (the reference implementation) so that
``freq_bt``, ``freq_btfdr``, ``freq_etfe``, and the Welch path of ``freq_map``
apply the SPEC §2.6 / §2.7 / §3.3 / §10.2 / §10.3 rules identically and cannot
drift apart (issues #143, #141).

Two entry points:

- :func:`input_excitation_degenerate` -- the whole-signal absolute
  input-excitation check of §10.3 (constant or zero input).
- :func:`regularize_response` -- the per-frequency ``Phi_u`` guard, ``NaN``
  substitution, noise-spectrum PSD clamp, and coherence handling of
  §2.6 / §2.7, for both SISO and MIMO.

The ``sigma_G = Inf`` sentinel of §3.3 is applied downstream in
``sid_uncertainty`` (which sees ``coherence < eps`` at degenerate
frequencies); it is not this module's job.
"""

from __future__ import annotations

import warnings

import numpy as np

# Relative floor / conditioning threshold shared by every estimator (§2.6).
EPS_REG: float = 1e-10

# Canonical near-singular warning message (matches the historical freq_bt text
# so downstream ``pytest.warns(match=...)`` and users see one wording).
_SINGULAR_MSG = (
    "Input spectrum Phi_u is near-singular at some frequencies. G set to NaN at those points."
)
_CONSTANT_INPUT_MSG = (
    "Input u has negligible variation (constant or zero). The frequency "
    "response is unidentifiable; G set to NaN and sigma_G to Inf everywhere."
)


def input_excitation_degenerate(u: np.ndarray | None, eps: float = EPS_REG) -> bool:
    """Return ``True`` when the input carries no usable excitation (SPEC.md §10.3).

    A constant input (zero variance about its own mean) or an identically-zero
    input cannot identify any dynamics on the ``(0, pi]`` grid. This is an
    *absolute* check on the sample variance relative to the mean square, so it
    fires even though a constant input has a nonzero ``Phi_u`` under the
    biased-covariance convention (which the relative per-frequency floor of
    §10.2 can never detect).

    Parameters
    ----------
    u : ndarray or None
        Input signal, shape ``(N,)``, ``(N, nu)``, or ``(N, nu, L)``. ``None``
        (time-series mode) is never degenerate in this sense.
    eps : float
        Relative tolerance; ``max_ch(var) <= eps * max_ch(mean_square)`` counts
        as constant. Default ``1e-10``.

    Returns
    -------
    bool
        ``True`` when the input is constant or identically zero.

    Examples
    --------
    >>> import numpy as np
    >>> input_excitation_degenerate(np.ones(100))  # doctest: +SKIP
    True

    See Also
    --------
    regularize_response : Per-frequency degenerate handling this pairs with.

    Changelog
    ---------
    2026-07-23 : First version by Pedro Lourenco.
    """
    if u is None:
        return False
    u = np.asarray(u, dtype=np.float64)
    if u.ndim == 1:
        u = u[:, np.newaxis]
    elif u.ndim == 3:
        # (N, nu, L): pool trajectories per channel before the variance check.
        u = u.transpose(0, 2, 1).reshape(-1, u.shape[1])
    var_ac = u.var(axis=0)  # per-channel variance about the mean
    mean_sq = (u**2).mean(axis=0)  # per-channel mean square
    max_ms = float(np.max(mean_sq)) if mean_sq.size else 0.0
    return bool(np.max(var_ac) <= eps * max_ms or max_ms <= np.finfo(np.float64).tiny)


def dead_input_channels(u: np.ndarray | None, eps: float = EPS_REG) -> np.ndarray:
    """Flag individual (near-)constant input channels (SPEC.md §10.3).

    A single constant channel among otherwise active channels does *not* make
    the whole input degenerate (:func:`input_excitation_degenerate` stays
    ``False`` because a healthy channel dominates), and it may not push
    ``cond(Phi_u)`` past ``1/eps`` at every frequency -- yet that channel's
    transfer-function column is unidentifiable. Callers use this to warn
    without NaNing the healthy channels.

    Parameters
    ----------
    u : ndarray or None
        Input signal, shape ``(N,)``, ``(N, nu)``, or ``(N, nu, L)``.
    eps : float
        Relative tolerance, as in :func:`input_excitation_degenerate`.

    Returns
    -------
    ndarray
        Boolean array ``(nu,)``; ``True`` where the channel is (near-)constant.

    Examples
    --------
    >>> import numpy as np
    >>> dead_input_channels(np.c_[np.ones(50), np.arange(50.0)])  # doctest: +SKIP
    array([ True, False])

    See Also
    --------
    input_excitation_degenerate : Whole-signal version.

    Changelog
    ---------
    2026-07-23 : First version by Pedro Lourenco.
    """
    if u is None:
        return np.zeros(0, dtype=bool)
    u = np.asarray(u, dtype=np.float64)
    if u.ndim == 1:
        u = u[:, np.newaxis]
    elif u.ndim == 3:
        u = u.transpose(0, 2, 1).reshape(-1, u.shape[1])
    var_ac = u.var(axis=0)
    mean_sq = (u**2).mean(axis=0)
    scale = float(np.max(mean_sq)) if mean_sq.size else 0.0
    tiny = np.finfo(np.float64).tiny
    return (var_ac <= eps * scale) | (mean_sq <= tiny)


def clamp_psd_scalar(phi_v: np.ndarray) -> np.ndarray:
    """Clamp a SISO/scalar noise spectrum to be non-negative (SPEC.md §2.7).

    Only *finite* negative values are clamped to 0; ``NaN`` entries (from a
    degenerate frequency) are preserved. ``np.maximum(NaN, 0)`` already returns
    ``NaN`` in NumPy, but this wrapper makes the intent explicit and mirrors
    the MATLAB port, where ``max(NaN, 0) == 0`` would otherwise erase the NaN.

    Parameters
    ----------
    phi_v : ndarray
        Scalar noise spectrum ``(nf,)``.

    Returns
    -------
    ndarray
        ``phi_v`` with finite negatives set to 0 and ``NaN`` preserved.

    Examples
    --------
    >>> import numpy as np
    >>> clamp_psd_scalar(np.array([1.0, -1e-9, np.nan]))  # doctest: +SKIP
    array([ 1.,  0., nan])

    See Also
    --------
    clamp_psd_matrix : MIMO (matrix) counterpart.

    Changelog
    ---------
    2026-07-23 : First version by Pedro Lourenco.
    """
    out = np.array(phi_v, dtype=np.float64, copy=True)
    neg = np.isfinite(out) & (out < 0.0)
    out[neg] = 0.0
    return out


def clamp_psd_matrix(V: np.ndarray) -> np.ndarray:
    """Clamp a single ``(ny, ny)`` noise-spectrum matrix to PSD (§2.7).

    The **real part is taken first**, then the matrix is symmetrised and any
    negative eigenvalues are zeroed. Taking the real part before ``eigh`` is
    the pre-refactor convention for the MIMO noise spectrum ``Φ̂_v`` and avoids
    handing a complex-symmetric (not Hermitian) matrix to ``eigh``, which would
    silently reinterpret one triangle. A matrix containing ``NaN`` (degenerate
    frequency) is returned unchanged so the NaN survives.

    Parameters
    ----------
    V : ndarray
        A single ``(ny, ny)`` noise-spectrum matrix (may be complex).

    Returns
    -------
    ndarray
        The real, symmetric, PSD-clamped ``(ny, ny)`` matrix (or the real part
        unchanged when it contains ``NaN``).

    Examples
    --------
    >>> import numpy as np
    >>> clamp_psd_matrix(np.array([[1.0, 0.0], [0.0, -1e-9]]))  # doctest: +SKIP
    array([[1., 0.], [0., 0.]])

    Notes
    -----
    **Specification:** SPEC.md §2.7 -- Noise Spectrum Estimate (PSD clamp).

    See Also
    --------
    clamp_psd_scalar : Scalar (SISO) counterpart.

    Changelog
    ---------
    2026-07-23 : First version by Pedro Lourenco.
    """
    V = np.real(np.asarray(V))  # real part first (parity with the old path)
    if not np.all(np.isfinite(V)):
        return V
    Vk = (V + V.T) / 2.0
    eigvals, eigvecs = np.linalg.eigh(Vk)
    if np.any(eigvals < 0):
        eigvals = np.maximum(eigvals, 0.0)
        Vk = eigvecs @ np.diag(eigvals) @ eigvecs.T
    return np.real(Vk)


def regularize_response(
    phi_yu: np.ndarray,
    phi_u: np.ndarray,
    phi_y: np.ndarray,
    *,
    force_degenerate: bool = False,
    eps_reg: float = EPS_REG,
    warn: bool = True,
) -> tuple[np.ndarray, np.ndarray, np.ndarray | None]:
    """Form ``G``, ``Phi_v`` and coherence with the shared degenerate handling.

    Implements SPEC.md §2.6 (regularization) and §2.7 (PSD clamp).

    Parameters
    ----------
    phi_yu, phi_u, phi_y : ndarray
        Cross- and auto-spectra. SISO: 1-D ``(nf,)`` each. MIMO: ``phi_yu``
        ``(nf, ny, nu)``, ``phi_u`` ``(nf, nu, nu)``, ``phi_y`` ``(nf, ny, ny)``.
    force_degenerate : bool
        When ``True`` (the §10.3 input-excitation check fired), *every*
        frequency is treated as degenerate: ``G = NaN`` everywhere,
        ``Phi_v = Phi_y`` (all output is unexplained), coherence ``0``.
    eps_reg : float
        Relative floor (SISO) / conditioning threshold (MIMO), §2.6.
    warn : bool
        Emit the near-singular warning when any frequency is degenerate.

    Returns
    -------
    G : ndarray
        Frequency response, ``(nf,)`` (SISO, complex) or ``(nf, ny, nu)``.
    Phi_v : ndarray
        Noise spectrum, ``(nf,)`` (SISO) or ``(nf, ny, ny)``.
    coherence : ndarray or None
        Squared coherence ``(nf,)`` for SISO; ``None`` for MIMO.

    Examples
    --------
    >>> import numpy as np
    >>> G, Phi_v, coh = regularize_response(  # doctest: +SKIP
    ...     phi_yu, phi_u, phi_y
    ... )

    See Also
    --------
    input_excitation_degenerate : Whole-signal check that sets ``force_degenerate``.
    sid._internal.uncertainty.sid_uncertainty : Applies the sigma_G = Inf sentinel.

    Changelog
    ---------
    2026-07-23 : First version by Pedro Lourenco.
    """
    is_siso = phi_yu.ndim == 1

    if is_siso:
        PhiYU = phi_yu
        PhiU = np.real(phi_u)
        PhiY = np.real(phi_y)

        PhiU_abs = np.abs(PhiU)
        PhiU_max = float(np.max(PhiU_abs)) if PhiU_abs.size > 0 else 0.0
        # Defense-in-depth: an all-zero Phi_u makes the relative mask vacuous
        # (0 < 0 is False), so treat it as fully degenerate even if a caller
        # forgot the input-excitation check upstream -- otherwise the §143 bug
        # shape (silent NaN, no warning) returns by omission.
        if force_degenerate or PhiU_max <= np.finfo(np.float64).tiny:
            singular_mask = np.ones_like(PhiU_abs, dtype=bool)
            force_degenerate = True
        else:
            singular_mask = PhiU_abs < eps_reg * PhiU_max

        G = np.empty_like(PhiYU)
        with np.errstate(divide="ignore", invalid="ignore"):
            G[:] = PhiYU / phi_u
        if np.any(singular_mask):
            G[singular_mask] = np.nan + 1j * np.nan

        with np.errstate(divide="ignore", invalid="ignore"):
            PhiV = PhiY - np.abs(PhiYU) ** 2 / PhiU
        if np.any(singular_mask):
            PhiV[singular_mask] = PhiY[singular_mask]
        PhiV = clamp_psd_scalar(PhiV)

        with np.errstate(divide="ignore", invalid="ignore"):
            Coh = np.abs(PhiYU) ** 2 / (PhiY * PhiU)
        if np.any(singular_mask):
            Coh[singular_mask] = 0.0
        Coh = np.clip(Coh, 0.0, 1.0)

        if warn and np.any(singular_mask):
            warnings.warn(
                _CONSTANT_INPUT_MSG if force_degenerate else _SINGULAR_MSG,
                stacklevel=2,
            )
        return G, PhiV, Coh

    # ---- MIMO ----
    nf = phi_yu.shape[0]
    ny = phi_yu.shape[1]
    nu = phi_yu.shape[2]
    PhiY = phi_y if phi_y.ndim == 3 else phi_y[:, np.newaxis, np.newaxis]
    PhiU = phi_u if phi_u.ndim == 3 else phi_u[:, np.newaxis, np.newaxis]
    PhiYU = phi_yu

    G = np.zeros((nf, ny, nu), dtype=complex)
    PhiV = np.zeros((nf, ny, ny))
    any_singular = False

    for k in range(nf):
        PhiU_k = PhiU[k, :, :].reshape(nu, nu)
        PhiYU_k = PhiYU[k, :, :].reshape(ny, nu)
        PhiY_k = np.real(PhiY[k, :, :].reshape(ny, ny))

        rc = 1.0 / np.linalg.cond(PhiU_k) if np.all(np.isfinite(PhiU_k)) else 0.0
        if force_degenerate or rc < eps_reg:
            G[k, :, :] = np.nan
            PhiV[k, :, :] = clamp_psd_matrix(PhiY_k)  # all output is noise
            any_singular = True
        else:
            Gk = np.linalg.solve(PhiU_k.T, PhiYU_k.T).T
            G[k, :, :] = Gk
            PhiV[k, :, :] = clamp_psd_matrix(PhiY_k - Gk @ PhiYU_k.conj().T)

    if warn and any_singular:
        warnings.warn(
            _CONSTANT_INPUT_MSG if force_degenerate else _SINGULAR_MSG,
            stacklevel=2,
        )
    return G, PhiV, None
