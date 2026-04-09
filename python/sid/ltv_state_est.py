# Copyright (c) 2026 Pedro Lourenco. All rights reserved.
# This code is released under the MIT License. See LICENSE file in the
# project root for full license information.
#
# This module is part of the Open Source System Identification Toolbox (SID).
# https://github.com/pdlourenco/sid-matlab

"""Batch LTV state estimation (RTS smoother)."""

from __future__ import annotations

import numpy as np

from sid._exceptions import SidError
from sid._internal.ltv_blk_tri_solve import blk_tri_solve


def ltv_state_est(
    Y: np.ndarray | list[np.ndarray],
    U: np.ndarray | list[np.ndarray],
    A: np.ndarray,
    B: np.ndarray,
    H: np.ndarray,
    *,
    R: np.ndarray | None = None,
    Q: np.ndarray | None = None,
) -> np.ndarray | list[np.ndarray]:
    """Estimate state trajectories for a partially observed LTV system.

    This is the Python port of ``sidLTVStateEst.m``.

    Estimates state trajectories for a discrete-time LTV system with
    partial observations by minimising

    .. math::

        J = \\sum_k \\|y(k) - H\\,x(k)\\|^2_{R^{-1}}
          + \\sum_k \\|x(k{+}1) - A(k)\\,x(k) - B(k)\\,u(k)\\|^2_{Q^{-1}}

    This is equivalent to a Rauch-Tung-Striebel (RTS) fixed-interval
    smoother.  Solved via a block tridiagonal forward-backward pass in
    O(N n^3) per trajectory.

    Parameters
    ----------
    Y : ndarray or list of ndarray
        Output data.  Accepted formats:

        - ``(N+1, py)`` -- single trajectory with *py* outputs
        - ``(N+1, py, L)`` -- *L* trajectories, same horizon *N*
        - ``[Y_1, Y_2, ...]`` -- list of *L* arrays ``(N_l+1, py)``
          for variable-length trajectories
    U : ndarray or list of ndarray
        Input data, matching format of *Y*:

        - ``(N, q)`` -- single trajectory with *q* inputs
        - ``(N, q, L)`` -- *L* trajectories
        - ``[U_1, U_2, ...]`` -- list of *L* arrays ``(N_l, q)``
    A : ndarray, shape ``(n, n, N)``
        Time-varying dynamics matrices.
    B : ndarray, shape ``(n, q, N)``
        Time-varying input matrices.
    H : ndarray, shape ``(py, n)``
        Observation matrix.
    R : ndarray, shape ``(py, py)`` or None, optional
        Measurement noise covariance (symmetric positive definite).
        Default: ``eye(py)``.
    Q : ndarray, shape ``(n, n)`` or None, optional
        Process noise covariance (symmetric positive definite).
        Default: ``eye(n)``.

    Returns
    -------
    X_hat : ndarray or list of ndarray
        Estimated states.

        - ``(N+1, n, L)`` when *Y* is an ndarray.
        - ``[X_1, X_2, ...]`` list of ``(N_l+1, n)`` when *Y* is a
          list.

    Raises
    ------
    SidError
        If *Y* is a list but *U* is not (code: ``'bad_input'``).
    SidError
        If data dimensions are inconsistent (code: ``'dim_mismatch'``).

    Examples
    --------
    State estimation with known dynamics:

    >>> import numpy as np
    >>> import sid  # doctest: +SKIP
    >>> X_hat = sid.ltv_state_est(Y, U, A, B, H)  # doctest: +SKIP

    With known measurement and process noise:

    >>> X_hat = sid.ltv_state_est(  # doctest: +SKIP
    ...     Y, U, A, B, H, R=R_meas, Q=Q_proc
    ... )

    Notes
    -----
    **Specification:** SPEC.md section 8.12 -- Output-COSMIC state step;
    docs/cosmic_output.md -- Appendix A

    **Algorithm:**

    1. Precompute ``R^{-1}``, ``Q^{-1}``, ``H^T R^{-1} H``,
       ``H^T R^{-1}``.
    2. Build shared block tridiagonal system:

       - Diagonal blocks ``S[k]`` encoding observation and dynamics
         information.
       - Off-diagonal blocks ``U_c[k] = -A(k)^T Q^{-1}``.

    3. Per trajectory, construct the right-hand side ``Theta[k]`` from
       the output data and input data.
    4. Solve the block tridiagonal system via
       :func:`~sid._internal.ltv_blk_tri_solve.blk_tri_solve`.
    5. Extract state estimates from the solution.

    Complexity: O(L * N * n^3).

    See Also
    --------
    sid._internal.ltv_blk_tri_solve.blk_tri_solve :
        Block tridiagonal solver used internally.
    sid.ltv_disc : LTV state-space identification via COSMIC.

    Changelog
    ---------
    2026-04-09 : First version (Python port) by Pedro Lourenco.
    """

    # ------------------------------------------------------------------
    # 1. Parse and validate inputs
    # ------------------------------------------------------------------
    Y_arr, U_arr, A, B, H, R_mat, Q_mat, N, n, py, q, L, is_var_len, horizons = _parse_inputs(
        Y, U, A, B, H, R, Q
    )

    # ------------------------------------------------------------------
    # 2. Precompute shared quantities
    # ------------------------------------------------------------------
    I_py = np.eye(py)
    I_n = np.eye(n)
    Rinv = np.linalg.solve(R_mat, I_py)
    Qinv = np.linalg.solve(Q_mat, I_n)
    HtRinvH = H.T @ Rinv @ H  # (n, n)
    HtRinv = H.T @ Rinv  # (n, py)

    # Number of blocks: K = N+1 (spec indices 0..N)
    K = N + 1

    # ------------------------------------------------------------------
    # 3. Build block tridiagonal system (SPEC.md section 8.14)
    # ------------------------------------------------------------------
    AtQinv = np.zeros((n, n, N))
    for j in range(N):
        AtQinv[:, :, j] = A[:, :, j].T @ Qinv

    # Diagonal blocks S_blk (length K)
    S_blk: list[np.ndarray] = [None] * K  # type: ignore[list-item]

    S_blk[0] = HtRinvH + AtQinv[:, :, 0] @ A[:, :, 0]
    for j in range(1, N):
        S_blk[j] = HtRinvH + Qinv + AtQinv[:, :, j] @ A[:, :, j]
    S_blk[K - 1] = HtRinvH + Qinv

    # Off-diagonal blocks Uc_blk (length N = K-1)
    Uc_blk: list[np.ndarray] = [None] * N  # type: ignore[list-item]
    for j in range(N):
        Uc_blk[j] = -AtQinv[:, :, j]

    # ------------------------------------------------------------------
    # 4. Solve per trajectory
    # ------------------------------------------------------------------
    if is_var_len:
        X_hat_list: list[np.ndarray] = [None] * L  # type: ignore[list-item]
    else:
        X_hat_arr = np.zeros((K, n, L))

    for l_idx in range(L):
        if is_var_len:
            Nl = horizons[l_idx]
            Kl = Nl + 1
            Yl = Y_arr[l_idx]  # type: ignore[index]
            Ul = U_arr[l_idx]  # type: ignore[index]
        else:
            Nl = N
            Kl = K
            Yl = Y_arr[:, :, l_idx]  # type: ignore[index]
            Ul = U_arr[:, :, l_idx]  # type: ignore[index]

        # Build RHS Theta for this trajectory
        Theta: list[np.ndarray] = [None] * Kl  # type: ignore[list-item]

        # j=0 (spec k=0): Theta_0 = H'R^{-1}y(0) - A(0)'Q^{-1}b(0)
        b_curr = B[:, :, 0] @ Ul[0, :]
        Theta[0] = HtRinv @ Yl[0, :] - AtQinv[:, :, 0] @ b_curr

        # j=1..Nl-1 (spec k=1..Nl-1)
        for j in range(1, Nl):
            b_prev = b_curr
            b_curr = B[:, :, j] @ Ul[j, :]
            Theta[j] = HtRinv @ Yl[j, :] + Qinv @ b_prev - AtQinv[:, :, j] @ b_curr

        # j=Kl-1 (spec k=Nl): Theta_Nl = H'R^{-1}y(Nl) + Q^{-1}b(Nl-1)
        Theta[Kl - 1] = HtRinv @ Yl[Nl, :] + Qinv @ b_curr

        # Solve (slice shared blocks to trajectory horizon)
        w, _ = blk_tri_solve(S_blk[:Kl], Uc_blk[:Nl], Theta)

        # Extract states
        if is_var_len:
            X_hat_l = np.zeros((Kl, n))
            for j in range(Kl):
                X_hat_l[j, :] = w[j]
            X_hat_list[l_idx] = X_hat_l
        else:
            for j in range(Kl):
                X_hat_arr[j, :, l_idx] = w[j]

    if is_var_len:
        return X_hat_list
    return X_hat_arr


# ======================================================================
# Private helpers
# ======================================================================


def _parse_inputs(
    Y: np.ndarray | list[np.ndarray],
    U: np.ndarray | list[np.ndarray],
    A: np.ndarray,
    B: np.ndarray,
    H: np.ndarray,
    R: np.ndarray | None,
    Q: np.ndarray | None,
) -> tuple:
    """Validate and parse inputs for :func:`ltv_state_est`.

    Returns
    -------
    tuple
        ``(Y, U, A, B, H, R, Q, N, n, py, q, L, is_var_len, horizons)``
    """

    py = H.shape[0]
    n = H.shape[1]
    q = B.shape[1]
    N = A.shape[2]

    if A.shape[0] != n or A.shape[1] != n:
        raise SidError(
            "dim_mismatch",
            f"A must be ({n}, {n}, N), got {A.shape}.",
        )
    if B.shape[0] != n:
        raise SidError(
            "dim_mismatch",
            f"B must have {n} rows (state dim), got {B.shape[0]}.",
        )

    is_var_len = isinstance(Y, list)

    if is_var_len:
        # -- Variable-length trajectory mode --
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

        horizons = np.empty(L, dtype=np.intp)
        for i in range(L):
            Yi = np.asarray(Y[i], dtype=np.float64)
            Ui = np.asarray(U[i], dtype=np.float64)

            if Yi.ndim == 1:
                Yi = Yi[:, np.newaxis]
            if Ui.ndim == 1:
                Ui = Ui[:, np.newaxis]

            Y[i] = Yi  # type: ignore[index]
            U[i] = Ui  # type: ignore[index]

            if Yi.shape[1] != py:
                raise SidError(
                    "dim_mismatch",
                    f"Y[{i}] has {Yi.shape[1]} columns but H has {py} rows.",
                )
            if Ui.shape[1] != q:
                raise SidError(
                    "dim_mismatch",
                    f"U[{i}] has {Ui.shape[1]} columns but B has {q} columns.",
                )
            Nl = Ui.shape[0]
            if Yi.shape[0] != Nl + 1:
                raise SidError(
                    "dim_mismatch",
                    f"Y[{i}] has {Yi.shape[0]} rows but U[{i}] has {Nl} (need N_l+1 and N_l).",
                )
            if Nl > N:
                raise SidError(
                    "dim_mismatch",
                    f"Trajectory {i} has horizon {Nl} > size(A,2)={N}.",
                )
            horizons[i] = Nl

    else:
        # -- Uniform-horizon mode --
        horizons = None

        Y_nd = np.asarray(Y, dtype=np.float64)
        U_nd = np.asarray(U, dtype=np.float64)

        if Y_nd.ndim == 2:
            Y_nd = Y_nd[:, :, np.newaxis]
        if U_nd.ndim == 2:
            U_nd = U_nd[:, :, np.newaxis]

        L = Y_nd.shape[2]

        if Y_nd.shape[0] != N + 1:
            raise SidError(
                "dim_mismatch",
                f"Y must have N+1={N + 1} rows, got {Y_nd.shape[0]}.",
            )
        if Y_nd.shape[1] != py:
            raise SidError(
                "dim_mismatch",
                f"Y has {Y_nd.shape[1]} columns but H has {py} rows.",
            )
        if U_nd.shape[0] != N:
            raise SidError(
                "dim_mismatch",
                f"U must have N={N} rows, got {U_nd.shape[0]}.",
            )
        if U_nd.shape[1] != q:
            raise SidError(
                "dim_mismatch",
                f"U has {U_nd.shape[1]} columns but B has {q} columns.",
            )
        if U_nd.shape[2] != L:
            raise SidError(
                "dim_mismatch",
                f"U has {U_nd.shape[2]} trajectories but Y has {L}.",
            )

        Y = Y_nd  # type: ignore[assignment]
        U = U_nd  # type: ignore[assignment]

    # -- Defaults for R and Q --
    if R is None:
        R = np.eye(py)
    else:
        R = np.asarray(R, dtype=np.float64)
    if Q is None:
        Q = np.eye(n)
    else:
        Q = np.asarray(Q, dtype=np.float64)

    if R.shape != (py, py):
        raise SidError(
            "dim_mismatch",
            f"R must be ({py}, {py}), got {R.shape}.",
        )
    if Q.shape != (n, n):
        raise SidError(
            "dim_mismatch",
            f"Q must be ({n}, {n}), got {Q.shape}.",
        )

    return Y, U, A, B, H, R, Q, N, n, py, q, L, is_var_len, horizons
