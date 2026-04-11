# Copyright (c) 2026 Pedro Lourenco. All rights reserved.
# This code is released under the MIT License. See LICENSE file in the
# project root for full license information.
#
# This module is part of the Open Source System Identification Toolbox (SID).
# https://github.com/pdlourenco/sid

"""Estimate model order from frequency response via Hankel SVD."""

from __future__ import annotations

import warnings

import numpy as np

from sid._exceptions import SidError


def model_order(
    result: object,
    *,
    horizon: int | None = None,
    threshold: float | None = None,
) -> tuple[int, dict]:
    """Estimate model order from frequency response data.

    This is the Python port of ``sidModelOrder.m``.

    Estimates the state dimension *n* of a linear system from the
    singular value decomposition of a block Hankel matrix built from
    impulse response coefficients.  The impulse response is obtained via
    IFFT of the frequency response estimate.

    Parameters
    ----------
    result : FreqResult
        Result from any ``sid.freq_*`` function.  Must expose
        ``.response``, ``.noise_spectrum``, and ``.frequency``
        attributes.  For time-series data (no input), ``response`` may be
        ``None``; in that case ``.noise_spectrum`` is used.
    horizon : int or None, optional
        Block Hankel prediction horizon *r*.
        Default: ``min(N_imp // 3, 50)``.
    threshold : float or None, optional
        If specified, count singular values with
        ``sigma_k / sigma_1 > threshold`` instead of gap detection.
        Default: ``None`` (use gap method).

    Returns
    -------
    n : int
        Estimated model order (state dimension).
    sv_dict : dict
        Dictionary with keys:

        - ``'singular_values'`` -- 1-D ndarray of singular values.
        - ``'horizon'`` -- int, prediction horizon used.

    Raises
    ------
    SidError
        If *result* is missing required attributes (code:
        ``'bad_input'``).
    SidError
        If *horizon* is not a positive integer (code:
        ``'bad_horizon'``).
    SidError
        If *threshold* is not a positive scalar (code:
        ``'bad_threshold'``).
    SidError
        If the impulse response is too short for a Hankel matrix (code:
        ``'too_short'``).

    Examples
    --------
    Automated model order detection:

    >>> import sid  # doctest: +SKIP
    >>> G = sid.freq_bt(y, u)  # doctest: +SKIP
    >>> n, sv = sid.model_order(G)  # doctest: +SKIP

    With a threshold:

    >>> n, sv = sid.model_order(G, threshold=0.01)  # doctest: +SKIP

    Notes
    -----
    **Specification:** SPEC.md section 8.12 -- Output-COSMIC: Partial
    State Observation

    **Algorithm:**

    1. Compute impulse response ``g(k)`` via IFFT of the frequency
       response, using a conjugate-symmetric extension to obtain a
       real-valued sequence.
    2. Build block Hankel matrix *H* with *r* block-rows and *r*
       block-columns.  For MIMO systems each entry ``g(k)`` is an
       ``(ny, nu)`` block.
    3. Compute SVD of *H*.
    4. Detect model order as ``argmax_k (sigma_k / sigma_{k+1})``
       (gap method) or count singular values above *threshold*.

    **References:**

    Kung, S.Y. "A new identification and model reduction algorithm via
    singular value decomposition." Proc. 12th Asilomar Conference, 1978.

    See Also
    --------
    sid.freq_bt : Blackman-Tukey frequency response estimation.
    sid.freq_etfe : Empirical transfer function estimate.

    Changelog
    ---------
    2026-04-09 : First version (Python port) by Pedro Lourenco.
    """

    # ------------------------------------------------------------------
    # 1. Validate input
    # ------------------------------------------------------------------
    if not hasattr(result, "frequency"):
        raise SidError("bad_input", "Result must have a 'frequency' attribute.")

    is_time_series = (
        not hasattr(result, "response") or result.response is None  # type: ignore[union-attr]
    )

    if is_time_series:
        if not hasattr(result, "noise_spectrum") or result.noise_spectrum is None:  # type: ignore[union-attr]
            raise SidError(
                "bad_input",
                "Result must have a 'response' or 'noise_spectrum' attribute.",
            )

    if horizon is not None:
        if not isinstance(horizon, (int, np.integer)) or horizon < 1:
            raise SidError("bad_horizon", "Horizon must be a positive integer.")

    if threshold is not None:
        if not np.isscalar(threshold) or threshold <= 0:
            raise SidError("bad_threshold", "Threshold must be a positive scalar.")

    # ------------------------------------------------------------------
    # 2. Extract frequency response
    # ------------------------------------------------------------------
    if is_time_series:
        G = np.asarray(result.noise_spectrum)  # type: ignore[union-attr]
    else:
        G = np.asarray(result.response)  # type: ignore[union-attr]

    # Determine dimensions
    if G.ndim == 1:
        nf = G.shape[0]
        ny = 1
        nu = 1
    elif G.ndim == 2:
        nf = G.shape[0]
        ny = G.shape[1]
        nu = 1
    else:
        nf = G.shape[0]
        ny = G.shape[1]
        nu = G.shape[2]

    # ------------------------------------------------------------------
    # 3. Compute impulse response via IFFT (SPEC.md section 8.13)
    # ------------------------------------------------------------------
    Nfft = 2 * nf
    g_all = np.zeros((Nfft, ny, nu))

    for iy in range(ny):
        for iu in range(nu):
            if G.ndim == 3:
                Gvec = G[:, iy, iu]
            elif G.ndim == 2:
                Gvec = G[:, iy]
            else:
                Gvec = G[:]

            # Build conjugate-symmetric full-circle spectrum
            Gfull = np.zeros(Nfft, dtype=complex)
            Gfull[0] = np.real(Gvec[0])  # DC approximation
            Gfull[1:nf] = Gvec[0 : nf - 1]  # w1 to w_{nf-1}
            Gfull[nf] = np.real(Gvec[nf - 1])  # Nyquist (real)
            Gfull[nf + 1 :] = np.conj(Gvec[nf - 2 :: -1])  # mirror

            g_all[:, iy, iu] = np.real(np.fft.ifft(Gfull))

    # Use causal part (first half) as impulse response coefficients
    N_imp = nf
    g = g_all[:N_imp, :, :]

    # ------------------------------------------------------------------
    # 4. Determine horizon
    # ------------------------------------------------------------------
    if horizon is None:
        horizon = min(N_imp // 3, 50)

    if horizon < 2:
        raise SidError(
            "too_short",
            f"Impulse response too short for Hankel matrix (need N_imp >= 6, got {N_imp}).",
        )

    # Need at least 2*horizon - 1 impulse response coefficients
    if 2 * horizon - 1 > N_imp:
        horizon = (N_imp + 1) // 2
        if horizon < 2:
            raise SidError(
                "too_short",
                "Impulse response too short for Hankel matrix.",
            )

    r = horizon

    # ------------------------------------------------------------------
    # 5. Build block Hankel matrix (SPEC.md section 8.13)
    # ------------------------------------------------------------------
    H_hankel = np.zeros((r * ny, r * nu))

    for bi in range(r):
        for bj in range(r):
            idx = bi + bj  # 0-based: bi + bj corresponds to g(bi+bj+1-1)
            if idx < N_imp:
                H_hankel[bi * ny : (bi + 1) * ny, bj * nu : (bj + 1) * nu] = g[idx, :, :].reshape(
                    ny, nu
                )

    # ------------------------------------------------------------------
    # 6. SVD
    # ------------------------------------------------------------------
    _, Sigma, _ = np.linalg.svd(H_hankel, full_matrices=False)
    sigmas = Sigma  # 1-D array of singular values

    n_sigma = len(sigmas)

    # ------------------------------------------------------------------
    # 7. Detect model order (SPEC.md section 8.13)
    # ------------------------------------------------------------------
    if threshold is not None:
        # Threshold method: count sigma_k / sigma_1 > tau
        if sigmas[0] < np.finfo(float).eps:
            n = 1
            warnings.warn(
                "All singular values near zero. Returning n = 1.",
                stacklevel=2,
            )
        else:
            n = int(np.sum(sigmas / sigmas[0] > threshold))
            n = max(n, 1)
    else:
        # Gap method: argmax_k (sigma_k / sigma_{k+1})
        if n_sigma < 2 or sigmas[0] < np.finfo(float).eps:
            n = 1
            if sigmas[0] < np.finfo(float).eps:
                warnings.warn(
                    "All singular values near zero. Returning n = 1.",
                    stacklevel=2,
                )
        else:
            # Only consider ratios among singular values above a noise
            # floor to avoid spurious gaps in the numerical tail.
            noise_floor = sigmas[0] * np.sqrt(n_sigma) * np.finfo(float).eps
            last_sig = np.where(sigmas > noise_floor)[0]
            last_sig = last_sig[-1] if len(last_sig) > 0 else 0
            max_k = min(last_sig, n_sigma // 2)
            max_k = max(max_k, 1)

            ratios = sigmas[:max_k] / sigmas[1 : max_k + 1]

            # Largest gap
            n = int(np.argmax(ratios)) + 1  # 1-based order

    # ------------------------------------------------------------------
    # 8. Pack output
    # ------------------------------------------------------------------
    sv_dict: dict = {
        "singular_values": sigmas,
        "horizon": r,
    }

    return n, sv_dict
