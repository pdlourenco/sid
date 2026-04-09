# Copyright (c) 2026 Pedro Lourenco. All rights reserved.
# This code is released under the MIT License. See LICENSE file in the
# project root for full license information.
#
# This module is part of the Open Source System Identification Toolbox (SID).
# https://github.com/pdlourenco/sid-matlab

"""Frozen transfer function from LTV state-space model."""

from __future__ import annotations

import numpy as np
from numpy.linalg import solve

from sid._exceptions import SidError
from sid._results import FrozenResult, LTVResult


def ltv_disc_frozen(
    ltv_result: LTVResult,
    *,
    frequencies: np.ndarray | None = None,
    time_steps: np.ndarray | None = None,
    sample_time: float = 1.0,
) -> FrozenResult:
    """Compute the frozen (instantaneous) transfer function from an LTV model.

    For each time step *k* and frequency *w*, evaluates the frozen
    transfer function

    .. math::

        G(w, k) = (e^{jw} I - A(k))^{-1} B(k)

    If the ``ltv_result`` includes Bayesian uncertainty (from
    :func:`sid.ltv_disc` with ``uncertainty=True``), the standard
    deviation of *G* is propagated via first-order (Jacobian)
    linearization using the rank-1 factorization described in
    SPEC.md S8.11.1.

    Parameters
    ----------
    ltv_result : LTVResult
        Result struct from :func:`sid.ltv_disc`.
    frequencies : ndarray of shape (nf,) or None, optional
        Frequency vector in rad/sample.  Default is 128 linearly
        spaced points in (0, pi].
    time_steps : ndarray of shape (nk,) or None, optional
        0-based indices of time steps to evaluate.  Default is all
        time steps ``np.arange(N)``.
    sample_time : float, optional
        Sample time in seconds, used only for the ``frequency_hz``
        output.  Default is ``1.0``.

    Returns
    -------
    FrozenResult
        Frozen dataclass with fields:

        - **frequency** (*ndarray, shape (nf,)*) -- Frequency in rad/sample.
        - **frequency_hz** (*ndarray, shape (nf,)*) -- Frequency in Hz.
        - **time_steps** (*ndarray, shape (nk,)*) -- Selected 0-based indices.
        - **response** (*ndarray, shape (nf, p, q, nk)*) -- Complex
          frozen transfer function.
        - **response_std** (*ndarray or None, shape (nf, p, q, nk)*) --
          Standard deviation of ``response``.  ``None`` when the input
          ``LTVResult`` has no uncertainty.
        - **sample_time** (*float*) -- Sample time.
        - **method** (*str*) -- ``'ltv_disc_frozen'``.

    Raises
    ------
    SidError
        If any element of ``time_steps`` is outside ``[0, N-1]``
        (code: ``'bad_time_steps'``).

    Examples
    --------
    Basic usage with default frequencies and all time steps:

    >>> import sid  # doctest: +SKIP
    >>> ltv = sid.ltv_disc(X, U, lambda_=1e5, uncertainty=True)  # doctest: +SKIP
    >>> frz = sid.ltv_disc_frozen(ltv)  # doctest: +SKIP

    Custom frequencies and selected time steps:

    >>> import numpy as np  # doctest: +SKIP
    >>> w = np.logspace(-2, np.log10(np.pi), 200)  # doctest: +SKIP
    >>> frz = sid.ltv_disc_frozen(ltv, frequencies=w, time_steps=np.array([0, 49, 99]))  # doctest: +SKIP

    Notes
    -----
    **Algorithm:**

    1. Build a frequency grid (default 128 points in (0, pi]).
    2. For each selected time step *k* and frequency *w*, compute the
       resolvent R = (z I - A(k))^{-1} via ``numpy.linalg.solve`` and
       the frozen response G(k) = R B(k).
    3. If the ``LTVResult`` carries uncertainty, propagate it using the
       rank-1 Jacobian factorization:

       .. math::

           \\operatorname{Var}(G_{ab}) =
               (v^H P(k) v) \\cdot (r_a \\Sigma r_a^H)

       where v = [G_k(:,b); e_b], r_a = R(a,:), P(k) is the row-wise
       posterior covariance, and Sigma is the noise covariance.

    **Specification:** SPEC.md S8.9

    References
    ----------
    .. [1] Carvalho, Soares, Lourenco, Ventura. "COSMIC: fast closed-form
       identification from large-scale data for LTV systems."
       arXiv:2112.04355, 2022.

    See Also
    --------
    sid.ltv_disc : LTV state-space identification.
    sid.freq_map : Time-varying frequency response via short-time windows.

    Changelog
    ---------
    2026-04-08 : First version (Python port) by Pedro Lourenco.
    """
    # ------------------------------------------------------------------
    # 1. Defaults
    # ------------------------------------------------------------------
    nf_default = 128
    if frequencies is None:
        w = np.arange(1, nf_default + 1) * (np.pi / nf_default)
    else:
        w = np.asarray(frequencies, dtype=np.float64).ravel()

    # ------------------------------------------------------------------
    # 2. Extract from ltv_result
    # ------------------------------------------------------------------
    A = ltv_result.a  # (p, p, N)
    B = ltv_result.b  # (p, q, N)
    p = ltv_result.state_dim
    q = ltv_result.input_dim
    N = ltv_result.data_length

    has_uncertainty = ltv_result.p_cov is not None

    if time_steps is None:
        k_vec = np.arange(N)
    else:
        k_vec = np.asarray(time_steps, dtype=np.intp).ravel()

    nk = len(k_vec)
    nf = len(w)

    # ------------------------------------------------------------------
    # 3. Validate time step indices (0-based)
    # ------------------------------------------------------------------
    if np.any(k_vec < 0) or np.any(k_vec > N - 1):
        raise SidError(
            "bad_time_steps",
            f"time_steps must be in range [0, {N - 1}].",
        )

    # ------------------------------------------------------------------
    # 4. Compute frozen transfer function
    # ------------------------------------------------------------------
    G = np.zeros((nf, p, q, nk), dtype=np.complex128)
    G_std: np.ndarray | None = None
    if has_uncertainty:
        G_std = np.zeros((nf, p, q, nk))

    Ip = np.eye(p)

    for ik in range(nk):
        ki = k_vec[ik]
        Ak = A[:, :, ki]  # (p, p)
        Bk = B[:, :, ki]  # (p, q)

        for iw in range(nf):
            z = np.exp(1j * w[iw])
            # R = (z I - A(k))^{-1}  via solve: (z I - Ak) R = I
            R = solve(z * Ip - Ak, Ip)  # (p, p)
            Gk = R @ Bk  # (p, q)
            G[iw, :, :, ik] = Gk

        # ---- Uncertainty propagation ----
        if has_uncertainty:
            Pk = ltv_result.p_cov[:, :, ki]  # (d, d), d = p + q
            Sigma = ltv_result.noise_cov  # (p, p)
            d = p + q

            for iw in range(nf):
                z = np.exp(1j * w[iw])
                R = solve(z * Ip - Ak, Ip)  # (p, p) resolvent
                Gk = R @ Bk  # (p, q)

                # Sigma quadratic form for each output row
                sig_quad = np.zeros(p)
                for a in range(p):
                    ra = R[a, :]  # (p,) complex
                    sig_quad[a] = np.real(ra @ Sigma @ ra.conj())

                # P quadratic form for each input column
                var_G = np.zeros((p, q))
                for b in range(q):
                    v = np.zeros(d, dtype=np.complex128)
                    v[:p] = Gk[:, b]
                    v[p + b] = 1.0
                    p_quad = np.real(v.conj() @ Pk @ v)  # scalar
                    var_G[:, b] = p_quad * sig_quad

                G_std[iw, :, :, ik] = np.sqrt(var_G)

    # ------------------------------------------------------------------
    # 5. Pack result
    # ------------------------------------------------------------------
    return FrozenResult(
        frequency=w,
        frequency_hz=w / (2.0 * np.pi * sample_time),
        time_steps=k_vec,
        response=G,
        response_std=G_std,
        sample_time=sample_time,
        method="ltv_disc_frozen",
    )
