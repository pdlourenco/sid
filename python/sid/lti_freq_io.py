# Copyright (c) 2026 Pedro Lourenco. All rights reserved.
# This code is released under the MIT License. See LICENSE file in the
# project root for full license information.
#
# This module is part of the Open Source System Identification Toolbox (SID).
# https://github.com/pdlourenco/sid

"""LTI state-space realization from input-output data via Ho-Kalman."""

from __future__ import annotations

import warnings

import numpy as np

from sid._exceptions import SidError
from sid.freq_bt import freq_bt


def lti_freq_io(
    Y: np.ndarray | list,
    U: np.ndarray | list,
    H: np.ndarray,
    *,
    horizon: int | None = None,
    max_stabilize: float = 0.999,
) -> tuple[np.ndarray, np.ndarray]:
    """Identify an LTI state-space model from input-output data.

    Estimates a constant (LTI) state-space realization from input-output
    data, given a known observation matrix *H*:

    .. math::

        x(k+1) = A_0\\,x(k) + B_0\\,u(k)

        y(k)   = H\\,x(k)

    Uses the Ho-Kalman realization algorithm applied to the frequency
    response estimated via :func:`sid.freq_bt`.  The realization is
    transformed to the H-basis so that the observation equation
    ``y = H x`` holds exactly with the returned ``(A0, B0)``.

    This is the Python port of ``sidLTIfreqIO.m``.

    Parameters
    ----------
    Y : ndarray or list of ndarray
        Output data.  Accepted formats:

        - ``(N+1, py)`` -- single trajectory
        - ``(N+1, py, L)`` -- *L* trajectories, same horizon *N*
        - ``[Y_1, Y_2, ...]`` -- list of arrays ``(N_l+1, py)``.
          Trajectories shorter than ``ceil(2/3 * max_horizon)`` are
          discarded; the rest are trimmed to a common length.
    U : ndarray or list of ndarray
        Input data, matching format of *Y*:

        - ``(N, q)`` -- single trajectory
        - ``(N, q, L)`` -- *L* trajectories
        - ``[U_1, U_2, ...]`` -- list of arrays ``(N_l, q)``
    H : ndarray, shape ``(py, n)``
        Observation matrix.
    horizon : int or None, optional
        Hankel matrix depth *r*.  Default is ``min(N_imp // 3, 50)``
        where *N_imp* is the number of impulse response coefficients.
    max_stabilize : float, optional
        Maximum eigenvalue magnitude after stabilization.  Unstable
        eigenvalues are reflected inside the unit circle, then clamped
        to this radius.  Default is ``0.999``.

    Returns
    -------
    A0 : ndarray, shape ``(n, n)``
        Estimated LTI dynamics matrix.
    B0 : ndarray, shape ``(n, q)``
        Estimated LTI input matrix.

    Raises
    ------
    SidError
        If data is too short for Hankel construction
        (code: ``'too_short'``).
    SidError
        If data dimensions are inconsistent
        (code: ``'dim_mismatch'``).

    Examples
    --------
    >>> import numpy as np
    >>> import sid  # doctest: +SKIP
    >>> H = np.array([[1, 0]])
    >>> A0, B0 = sid.lti_freq_io(Y, U, H)  # doctest: +SKIP

    Notes
    -----
    **Algorithm:**

    1. Estimate transfer function G(e^{jw}) via :func:`sid.freq_bt`.
    2. Compute Markov parameters g(k) = H A^{k-1} B via IFFT.
    3. Build block Hankel matrices H_0 and H_1 (shifted).
    4. SVD of H_0, truncate to order *n* (Ho-Kalman realization).
    5. Transform realization to H-basis: find T s.t. C_r T^{-1} = H.
    6. Stabilize eigenvalues if needed.

    **Specification:** SPEC.md section 8.12 -- Output-COSMIC (LTI
    initialization)

    References
    ----------
    .. [1] Ho, B.L. and Kalman, R.E. "Effective construction of linear
       state-variable models from input/output functions."
       Regelungstechnik, 14(12):545-548, 1966.

    See Also
    --------
    sid.freq_bt : Frequency response estimation used in step 1.
    sid.ltv_disc_io : Output-COSMIC identification that uses this as
        initializer.

    Changelog
    ---------
    2026-04-09 : First version (Python port) by Pedro Lourenco.
    """
    # ------------------------------------------------------------------
    # 1. Parse inputs
    # ------------------------------------------------------------------
    Y_trim, U_out, py, n, q = _parse_inputs(Y, U, H)

    # ------------------------------------------------------------------
    # 2. Frequency response estimation (SPEC.md S8.13)
    # ------------------------------------------------------------------
    G_result = freq_bt(Y_trim, U_out)
    G = G_result.response  # (nf, py, q) or (nf,) complex
    nf = G.shape[0]

    # Ensure G is 3-D: (nf, py, q)
    if G.ndim == 1:
        G = G[:, np.newaxis, np.newaxis]
    elif G.ndim == 2:
        G = G[:, :, np.newaxis]

    # ------------------------------------------------------------------
    # 3. Impulse response via IFFT (SPEC.md S8.13)
    # ------------------------------------------------------------------
    g = _freq_to_impulse(G, nf, py, q)
    N_imp = g.shape[0]

    # ------------------------------------------------------------------
    # 4. Determine horizon
    # ------------------------------------------------------------------
    if horizon is None:
        horizon = min(N_imp // 3, 50)

    if horizon < 2:
        raise SidError(
            "too_short",
            f"Data too short for Hankel matrix (need N_imp >= 6, got {N_imp}).",
        )

    # Need at least 2*horizon impulse response coefficients for H_1
    if 2 * horizon > N_imp:
        horizon = N_imp // 2
        if horizon < 2:
            raise SidError(
                "too_short",
                "Data too short for shifted Hankel matrix.",
            )

    r = horizon

    # Check that Hankel matrix can support order n
    if r * py < n or r * q < n:
        raise SidError(
            "too_short",
            f"Hankel size ({r}*{py} x {r}*{q}) too small for "
            f"order n={n}. Increase data length or horizon.",
        )

    # ------------------------------------------------------------------
    # 5. Build block Hankel matrices (SPEC.md S8.13)
    # ------------------------------------------------------------------
    H0, H1 = _build_hankel(g, r, py, q, N_imp)

    # ------------------------------------------------------------------
    # 6. Ho-Kalman SVD realization (SPEC.md S8.13)
    # ------------------------------------------------------------------
    A_r, B_r, C_r = _ho_kalman(H0, H1, n, py, q)

    # ------------------------------------------------------------------
    # 7. Transform to H-basis
    # ------------------------------------------------------------------
    A0, B0 = _transform_to_h_basis(A_r, B_r, C_r, H, n)

    # ------------------------------------------------------------------
    # 8. Stabilize
    # ------------------------------------------------------------------
    A0 = _stabilize(A0, max_stabilize)

    return A0, B0


# ======================================================================
# Private helpers
# ======================================================================


def _parse_inputs(
    Y: np.ndarray | list,
    U: np.ndarray | list,
    H: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, int, int, int]:
    """Validate and parse inputs for :func:`lti_freq_io`.

    For cell array (list) inputs, trajectories shorter than
    ``ceil(2/3 * max_horizon)`` are discarded and the rest are trimmed
    to a common length.

    Returns
    -------
    tuple
        ``(Y_trim, U_out, py, n, q)``
    """
    H = np.asarray(H, dtype=np.float64)
    if H.ndim != 2:
        raise SidError("dim_mismatch", "H must be a 2-D matrix (py x n).")
    py = H.shape[0]
    n = H.shape[1]

    if isinstance(Y, list):
        # -- Variable-length: trim and stack into 3D --
        if not isinstance(U, list):
            raise SidError(
                "bad_input",
                "When Y is a list, U must also be a list.",
            )
        L = len(Y)
        if len(U) != L:
            raise SidError(
                "dim_mismatch",
                f"Y has {L} trajectories but U has {len(U)}.",
            )
        if L == 0:
            raise SidError("bad_input", "Trajectory lists must not be empty.")

        Y_list = [np.asarray(y, dtype=np.float64) for y in Y]
        U_list = [np.asarray(u, dtype=np.float64) for u in U]

        # Ensure 2-D
        for i in range(L):
            if Y_list[i].ndim == 1:
                Y_list[i] = Y_list[i][:, np.newaxis]
            if U_list[i].ndim == 1:
                U_list[i] = U_list[i][:, np.newaxis]

        q = U_list[0].shape[1]

        horizons = np.array([u.shape[0] for u in U_list])

        for i in range(L):
            if Y_list[i].shape[1] != py:
                raise SidError(
                    "dim_mismatch",
                    f"Y[{i}] has {Y_list[i].shape[1]} columns but H has {py} rows.",
                )
            if U_list[i].shape[1] != q:
                raise SidError(
                    "dim_mismatch",
                    f"U[{i}] has {U_list[i].shape[1]} columns, expected {q}.",
                )

        # Keep trajectories with horizon >= ceil(2/3 * max)
        max_h = int(np.max(horizons))
        threshold = int(np.ceil(2 * max_h / 3))
        keep = [i for i in range(L) if horizons[i] >= threshold]
        if len(keep) == 0:
            raise SidError(
                "too_short",
                "No trajectories long enough for spectral estimation.",
            )

        N_common = int(min(horizons[i] for i in keep))
        L_kept = len(keep)

        # Stack into 3D arrays, trimmed to common length
        # Y_trim: (N_common, py, L_kept), U_out: (N_common, q, L_kept)
        Y_trim = np.zeros((N_common, py, L_kept))
        U_out = np.zeros((N_common, q, L_kept))
        for ii, idx in enumerate(keep):
            Y_trim[:, :, ii] = Y_list[idx][:N_common, :]
            U_out[:, :, ii] = U_list[idx][:N_common, :]

    else:
        # -- Uniform-horizon mode --
        Y_arr = np.asarray(Y, dtype=np.float64)
        U_arr = np.asarray(U, dtype=np.float64)

        # Ensure 3D
        if Y_arr.ndim == 1:
            Y_arr = Y_arr[:, np.newaxis, np.newaxis]
        elif Y_arr.ndim == 2:
            Y_arr = Y_arr[:, :, np.newaxis]

        if U_arr.ndim == 1:
            U_arr = U_arr[:, np.newaxis, np.newaxis]
        elif U_arr.ndim == 2:
            U_arr = U_arr[:, :, np.newaxis]

        N_u = U_arr.shape[0]
        q = U_arr.shape[1]

        if Y_arr.shape[0] < N_u:
            raise SidError(
                "dim_mismatch",
                f"Y must have at least N={N_u} rows, got {Y_arr.shape[0]}.",
            )
        if Y_arr.shape[1] != py:
            raise SidError(
                "dim_mismatch",
                f"Y has {Y_arr.shape[1]} columns but H has {py} rows.",
            )
        if U_arr.shape[2] != Y_arr.shape[2]:
            raise SidError(
                "dim_mismatch",
                f"U has {U_arr.shape[2]} trajectories but Y has {Y_arr.shape[2]}.",
            )

        # Trim Y to first N rows to match U
        Y_trim = Y_arr[:N_u, :, :]
        U_out = U_arr

        # Squeeze single-trajectory dimension so freq_bt sees 2D input
        if Y_trim.shape[2] == 1:
            Y_trim = Y_trim[:, :, 0]
            U_out = U_out[:, :, 0]

    return Y_trim, U_out, py, n, q


def _freq_to_impulse(
    G: np.ndarray,
    nf: int,
    ny: int,
    nu: int,
) -> np.ndarray:
    """Convert frequency response to impulse response via IFFT.

    Uses DC extrapolation and skips lag 0 (no feedthrough).

    Parameters
    ----------
    G : ndarray, shape ``(nf, ny, nu)``
        Complex frequency response.
    nf : int
        Number of frequency bins.
    ny : int
        Number of output channels.
    nu : int
        Number of input channels.

    Returns
    -------
    g : ndarray, shape ``(N_imp, ny, nu)``
        Impulse response coefficients for lags 1..N_imp, where
        ``N_imp = nf - 1``.
    """
    Nfft = 2 * nf
    g_all = np.zeros((Nfft, ny, nu))

    for iy in range(ny):
        for iu in range(nu):
            Gvec = G[:, iy, iu]  # (nf,) complex

            # Build full-circle: DC, positive freqs, Nyquist, mirror
            Gfull = np.zeros(Nfft, dtype=complex)

            # DC: extrapolate from first two bins (DC must be real)
            if nf >= 2:
                Gfull[0] = np.real(2 * Gvec[0] - Gvec[1])
            else:
                Gfull[0] = np.real(Gvec[0])

            Gfull[1:nf] = Gvec[: nf - 1]  # w_1 to w_{nf-1}
            Gfull[nf] = np.real(Gvec[nf - 1])  # Nyquist (real)
            Gfull[nf + 1 :] = np.conj(Gvec[nf - 2 :: -1])  # mirror

            g_all[:, iy, iu] = np.real(np.fft.ifft(Gfull))

    # Use causal part starting from lag 1 (skip direct feedthrough at
    # lag 0).  g_all[0] is lag 0, g_all[1] is lag 1 = H*B, etc.
    nf - 1
    g = g_all[1:nf, :, :]
    return g


def _build_hankel(
    g: np.ndarray,
    r: int,
    py: int,
    q: int,
    N_imp: int,
) -> tuple[np.ndarray, np.ndarray]:
    """Build block Hankel matrix H_0 and shifted Hankel H_1.

    H_0{i,j} = g(i+j-1), H_1{i,j} = g(i+j), using 1-based block
    indexing (translated to 0-based array indexing).

    Parameters
    ----------
    g : ndarray, shape ``(N_imp, py, q)``
        Impulse response coefficients (1-indexed in MATLAB).
    r : int
        Hankel matrix depth.
    py : int
        Output dimension.
    q : int
        Input dimension.
    N_imp : int
        Number of impulse response coefficients.

    Returns
    -------
    H0 : ndarray, shape ``(r*py, r*q)``
        Block Hankel matrix.
    H1 : ndarray, shape ``(r*py, r*q)``
        Shifted block Hankel matrix.
    """
    H0 = np.zeros((r * py, r * q))
    H1 = np.zeros((r * py, r * q))

    for bi in range(1, r + 1):
        for bj in range(1, r + 1):
            idx0 = bi + bj - 1  # 1-based lag index
            idx1 = bi + bj  # 1-based lag index

            row_s = (bi - 1) * py
            row_e = bi * py
            col_s = (bj - 1) * q
            col_e = bj * q

            if idx0 <= N_imp:
                H0[row_s:row_e, col_s:col_e] = g[idx0 - 1, :, :].reshape(py, q)
            if idx1 <= N_imp:
                H1[row_s:row_e, col_s:col_e] = g[idx1 - 1, :, :].reshape(py, q)

    return H0, H1


def _ho_kalman(
    H0: np.ndarray,
    H1: np.ndarray,
    n: int,
    py: int,
    q: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Ho-Kalman realization from Hankel matrices.

    Given H_0 and H_1, compute the minimal realization (A_r, B_r, C_r)
    of order *n* such that G(z) = C_r (zI - A_r)^{-1} B_r.

    Parameters
    ----------
    H0 : ndarray, shape ``(r*py, r*q)``
        Block Hankel matrix.
    H1 : ndarray, shape ``(r*py, r*q)``
        Shifted block Hankel matrix.
    n : int
        System order.
    py : int
        Output dimension.
    q : int
        Input dimension.

    Returns
    -------
    A_r : ndarray, shape ``(n, n)``
        Dynamics matrix in arbitrary basis.
    B_r : ndarray, shape ``(n, q)``
        Input matrix in arbitrary basis.
    C_r : ndarray, shape ``(py, n)``
        Output matrix in arbitrary basis.
    """
    U_svd, sigmas_full, Vt_svd = np.linalg.svd(H0, full_matrices=False)

    # sigmas_full has length min(r*py, r*q)
    if len(sigmas_full) < n:
        raise SidError(
            "too_short",
            f"Hankel SVD has only {len(sigmas_full)} singular values but n={n} requested.",
        )

    # Truncate to order n
    U_n = U_svd[:, :n]
    V_n = Vt_svd[:n, :].T  # numpy svd returns V^T; V_n is (r*q, n)

    S_sqrt = np.diag(np.sqrt(sigmas_full[:n]))
    S_isqrt = np.diag(1.0 / np.sqrt(sigmas_full[:n]))

    # Realization in arbitrary basis
    A_r = S_isqrt @ U_n.T @ H1 @ V_n @ S_isqrt  # (n, n)
    C_r = U_n[:py, :] @ S_sqrt  # (py, n)
    B_r = S_sqrt @ V_n[:q, :].T  # (n, q)

    return A_r, B_r, C_r


def _transform_to_h_basis(
    A_r: np.ndarray,
    B_r: np.ndarray,
    C_r: np.ndarray,
    H: np.ndarray,
    n: int,
) -> tuple[np.ndarray, np.ndarray]:
    """Transform realization to the H-basis.

    Find T such that C_r T^{-1} = H, then A0 = T A_r T^{-1},
    B0 = T B_r.

    Parameters
    ----------
    A_r : ndarray, shape ``(n, n)``
        Dynamics in arbitrary basis.
    B_r : ndarray, shape ``(n, q)``
        Input in arbitrary basis.
    C_r : ndarray, shape ``(py, n)``
        Output in arbitrary basis.
    H : ndarray, shape ``(py, n)``
        Desired observation matrix.
    n : int
        System order.

    Returns
    -------
    A0 : ndarray, shape ``(n, n)``
        Dynamics in H-basis.
    B0 : ndarray, shape ``(n, q)``
        Input in H-basis.
    """
    Cr_pinv = np.linalg.pinv(C_r)
    Tinv = Cr_pinv @ H + np.eye(n) - Cr_pinv @ C_r

    # Check conditioning
    rc = 1.0 / np.linalg.cond(Tinv)
    if rc < np.finfo(float).eps * 1e3:
        warnings.warn(
            f"Basis transform Tinv is near-singular (rcond={rc:.2e}). Using raw realization.",
            stacklevel=2,
        )
        return np.real(A_r), np.real(B_r)

    T = np.linalg.solve(Tinv, np.eye(n))
    A0 = np.real(T @ A_r @ Tinv)
    B0 = np.real(T @ B_r)
    return A0, B0


def _stabilize(A: np.ndarray, max_stab: float) -> np.ndarray:
    """Reflect unstable eigenvalues inside the unit circle.

    Eigenvalues with |lambda| > 1 are reflected: lambda <- 1/conj(lambda).
    The result is then clamped so that |lambda| <= max_stab.

    Parameters
    ----------
    A : ndarray, shape ``(n, n)``
        Dynamics matrix.
    max_stab : float
        Maximum permitted eigenvalue magnitude.

    Returns
    -------
    A : ndarray, shape ``(n, n)``
        Stabilized dynamics matrix.
    """
    n = A.shape[0]
    eigvals, V = np.linalg.eig(A)
    mags = np.abs(eigvals)

    if np.max(mags) <= max_stab:
        return A

    # Step 1: Reflect eigenvalues outside the unit circle
    outside = mags > 1
    eigvals[outside] = 1.0 / np.conj(eigvals[outside])
    mags[outside] = np.abs(eigvals[outside])

    # Step 2: Clamp to max_stab radius
    toolarge = mags > max_stab
    eigvals[toolarge] = max_stab * eigvals[toolarge] / mags[toolarge]

    A = np.real(V @ np.diag(eigvals) @ np.linalg.solve(V, np.eye(n)))
    return A
