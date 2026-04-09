# Copyright (c) 2026 Pedro Lourenco. All rights reserved.
# This code is released under the MIT License. See LICENSE file in the
# project root for full license information.
#
# This module is part of the Open Source System Identification Toolbox (SID).
# https://github.com/pdlourenco/sid-matlab

"""Discrete-time LTV state-space identification from partial observations."""

from __future__ import annotations


import numpy as np

from sid._exceptions import SidError
from sid._internal.ltv_build_block_terms import build_block_terms
from sid._internal.ltv_build_data_matrices import (
    build_data_matrices,
    build_data_matrices_var_len,
)
from sid._internal.ltv_cosmic_solve import cosmic_solve
from sid._internal.ltv_uncertainty_backward_pass import uncertainty_backward_pass
from sid._internal.estimate_noise_cov import estimate_noise_cov
from sid._internal.extract_std import extract_std
from sid._results import LTVIOResult
from sid.lti_freq_io import lti_freq_io
from sid.ltv_state_est import ltv_state_est


def ltv_disc_io(
    Y: np.ndarray | list,
    U: np.ndarray | list,
    H: np.ndarray,
    *,
    lambda_: float | np.ndarray,
    R: np.ndarray | None = None,
    max_iter: int = 50,
    tolerance: float = 1e-6,
    covariance_mode: str = "diagonal",
    trust_region: float | str = "off",
    trust_region_tol: float = 1e-6,
    uncertainty: bool = True,
) -> LTVIOResult:
    """Identify a discrete-time LTV system from partial observations.

    Identifies time-varying system matrices A(k), B(k) and estimates
    state trajectories from input-output data when only partial state
    observations are available:

    .. math::

        x(k+1) = A(k)\\,x(k) + B(k)\\,u(k), \\quad k = 0, \\ldots, N-1

        y(k) = H\\,x(k)

    Uses the Output-COSMIC algorithm: alternating minimisation between
    a COSMIC step (dynamics estimation) and an RTS smoother (state
    estimation), initialised via LTI realization from the I/O transfer
    function (:func:`sid.lti_freq_io`).

    When ``rank(H) == n`` (full state observation), this reduces to the
    standard COSMIC algorithm with no EM iterations.

    This is the Python port of ``sidLTVdiscIO.m``.

    Parameters
    ----------
    Y : ndarray or list of ndarray
        Output data.  Accepted formats:

        - ``(N+1, py)`` -- single trajectory
        - ``(N+1, py, L)`` -- *L* trajectories, same horizon
        - ``[Y_1, Y_2, ...]`` -- list of arrays ``(N_l+1, py)``
    U : ndarray or list of ndarray
        Input data, matching format of *Y*:

        - ``(N, q)`` -- single trajectory
        - ``(N, q, L)`` -- *L* trajectories
        - ``[U_1, U_2, ...]`` -- list of arrays ``(N_l, q)``
    H : ndarray, shape ``(py, n)``
        Observation matrix.
    lambda_ : float or ndarray, shape ``(N-1,)``
        Regularization strength.  Scalar (uniform) or per-step vector.
        Must be positive.
    R : ndarray, shape ``(py, py)`` or None, optional
        Measurement noise covariance (SPD).  Default: ``eye(py)``.
    max_iter : int, optional
        Maximum alternating iterations.  Default: ``50``.
    tolerance : float, optional
        Convergence tolerance on relative cost change.
        Default: ``1e-6``.
    covariance_mode : str, optional
        Noise covariance estimation mode: ``'diagonal'``, ``'full'``,
        or ``'isotropic'``.  Default: ``'diagonal'``.
    trust_region : float or ``'off'``, optional
        Trust-region parameter mu_0 in ``[0, 1]``, or ``'off'``.
        Default: ``'off'``.
    trust_region_tol : float, optional
        Minimum mu before final pass.  Default: ``1e-6``.
    uncertainty : bool, optional
        Compute Bayesian posterior uncertainty for A(k), B(k).
        Default: ``True``.  Uncertainty is always computed in
        Output-COSMIC; this parameter is accepted for API
        consistency.

    Returns
    -------
    LTVIOResult
        Frozen dataclass with fields:

        - **a** (*ndarray, shape (n, n, N)*) -- Time-varying dynamics.
        - **b** (*ndarray, shape (n, q, N)*) -- Time-varying input.
        - **x** (*ndarray or list*) -- Estimated state trajectories.
        - **h** (*ndarray, shape (py, n)*) -- Observation matrix (copy).
        - **r** (*ndarray, shape (py, py)*) -- Noise covariance used.
        - **cost** (*ndarray, shape (n_iter,)*) -- Cost at each iteration.
        - **iterations** (*int*) -- Number of alternating iterations.
        - **lambda_** (*ndarray, shape (N-1,)*) -- Lambda vector used.
        - **data_length** (*int*) -- *N*.
        - **state_dim** (*int*) -- *n*.
        - **output_dim** (*int*) -- *py*.
        - **input_dim** (*int*) -- *q*.
        - **num_trajectories** (*int*) -- *L*.
        - **a_std** (*ndarray or None*) -- Std dev of ``a``.
        - **b_std** (*ndarray or None*) -- Std dev of ``b``.
        - **p_cov** (*ndarray or None*) -- Row covariance blocks.
        - **noise_cov** (*ndarray or None*) -- Noise covariance.
        - **noise_cov_estimated** (*bool or None*) -- Whether estimated.
        - **noise_variance** (*float or None*) -- ``trace(Sigma)/n``.
        - **degrees_of_freedom** (*float or None*) -- Effective DOF.
        - **algorithm** (*str*) -- ``'cosmic'``.
        - **method** (*str*) -- ``'ltv_disc_io'``.

    Raises
    ------
    SidError
        If data dimensions are inconsistent (code: ``'dim_mismatch'``).
    SidError
        If lambda is invalid (code: ``'bad_lambda'``).
    SidError
        If covariance_mode is invalid (code: ``'bad_cov_mode'``).

    Examples
    --------
    Basic usage:

    >>> import numpy as np
    >>> import sid  # doctest: +SKIP
    >>> H = np.array([[1, 0]])
    >>> result = sid.ltv_disc_io(Y, U, H, lambda_=1e5)  # doctest: +SKIP

    With known measurement noise:

    >>> result = sid.ltv_disc_io(  # doctest: +SKIP
    ...     Y, U, H, lambda_=1e5, R=R_meas
    ... )

    With trust-region for difficult convergence:

    >>> result = sid.ltv_disc_io(  # doctest: +SKIP
    ...     Y, U, H, lambda_=1e5, trust_region=1.0
    ... )

    Notes
    -----
    **Algorithm (Output-COSMIC):**

    1. LTI initialisation: estimate A0, B0 via Ho-Kalman realization
       of the I/O transfer function (:func:`sid.lti_freq_io`).
    2. State step: fix dynamics, RTS smoother for states.
    3. COSMIC step: fix states, solve for A(k), B(k).
    4. Repeat steps 2-3 until convergence.

    Complexity: O(T * (L*N*n^3 + N*(n+q)^3)) where T is iterations.

    **Specification:** SPEC.md section 8.12 -- Output-COSMIC

    References
    ----------
    .. [1] Carvalho, Soares, Lourenco, Ventura. "COSMIC: fast
       closed-form identification from large-scale data for LTV
       systems." arXiv:2112.04355, 2022.

    See Also
    --------
    sid.lti_freq_io : LTI initializer used internally.
    sid.ltv_disc : LTV identification with full state observation.
    sid.ltv_state_est : State estimation via RTS smoother.

    Changelog
    ---------
    2026-04-09 : First version (Python port) by Pedro Lourenco.
    """
    # ------------------------------------------------------------------
    # 1. Parse and validate inputs
    # ------------------------------------------------------------------
    (
        Y_parsed,
        U_parsed,
        H,
        lambda_vec,
        R_mat,
        do_trust_region,
        mu,
        mu_tol,
        cov_mode,
        N,
        n,
        py,
        q,
        L,
        is_var_len,
        horizons,
    ) = _parse_inputs(
        Y,
        U,
        H,
        lambda_=lambda_,
        R=R,
        max_iter=max_iter,
        tolerance=tolerance,
        covariance_mode=covariance_mode,
        trust_region=trust_region,
        trust_region_tol=trust_region_tol,
    )

    # Precompute observation precision
    Rinv = np.linalg.solve(R_mat, np.eye(py))

    # ------------------------------------------------------------------
    # 2. Full-rank fast path (SPEC.md S8.12.2)
    # ------------------------------------------------------------------
    if np.linalg.matrix_rank(H) == n:
        Hpinv = np.linalg.solve(H.T @ Rinv @ H, H.T @ Rinv)  # (n, py)

        if is_var_len:
            X_hat: np.ndarray | list = [Y_parsed[traj] @ Hpinv.T for traj in range(L)]
        else:
            X_hat = np.zeros((N + 1, n, L))
            for traj in range(L):
                X_hat[:, :, traj] = Y_parsed[:, :, traj].reshape(N + 1, py) @ Hpinv.T

        A, B, S_c, D_c, Xl_c, C_c = _cosmic_step(
            X_hat, U_parsed, lambda_vec, N, n, q, L, is_var_len, horizons
        )

        J = _evaluate_full_cost(
            X_hat,
            A,
            B,
            Y_parsed,
            U_parsed,
            H,
            Rinv,
            lambda_vec,
            N,
            n,
            q,
            L,
            is_var_len,
            horizons,
        )

        # Add uncertainty
        A_std, B_std, P, Sigma, noise_var, dof = _add_uncertainty(
            S_c, D_c, Xl_c, C_c, lambda_vec, N, n, q, cov_mode
        )

        return LTVIOResult(
            a=A,
            b=B,
            x=X_hat,
            h=H.copy(),
            r=R_mat.copy(),
            cost=np.array([J]),
            iterations=0,
            lambda_=lambda_vec,
            data_length=N,
            state_dim=n,
            output_dim=py,
            input_dim=q,
            num_trajectories=L,
            a_std=A_std,
            b_std=B_std,
            p_cov=P,
            noise_cov=Sigma,
            noise_cov_estimated=True,
            noise_variance=noise_var,
            degrees_of_freedom=dof,
            algorithm="cosmic",
            method="ltv_disc_io",
        )

    # ------------------------------------------------------------------
    # 3. LTI Initialisation (SPEC.md S8.12.4)
    # ------------------------------------------------------------------
    A0, B0 = lti_freq_io(Y, U, H)
    A = np.tile(A0[:, :, np.newaxis], (1, 1, N))
    B = np.tile(B0[:, :, np.newaxis], (1, 1, N))
    A0_rep = np.tile(A0[:, :, np.newaxis], (1, 1, N))  # trust-region target

    # ------------------------------------------------------------------
    # 4. Alternating state-COSMIC loop (SPEC.md S8.12.3)
    # ------------------------------------------------------------------
    mu_current = float(do_trust_region) * mu
    cost_history: list[float] = []
    n_iter = 0

    # For trust-region accept/reject
    J_converged_mu = np.inf
    X_best: np.ndarray | list | None = None
    A_best = A.copy()
    B_best = B.copy()

    # Hold intermediates from last COSMIC step (for uncertainty)
    S_c = D_c = Xl_c = C_c = None

    for _iter in range(max_iter):
        # -- E-step: state estimation --
        if mu_current > 0:
            A_use = (1 - mu_current) * A + mu_current * A0_rep
        else:
            A_use = A

        X_hat = ltv_state_est(Y_parsed, U_parsed, A_use, B, H, R=R_mat)

        # -- M-step: COSMIC solve --
        A, B, S_c, D_c, Xl_c, C_c = _cosmic_step(
            X_hat, U_parsed, lambda_vec, N, n, q, L, is_var_len, horizons
        )

        # -- Evaluate cost and check convergence --
        J = _evaluate_full_cost(
            X_hat,
            A,
            B,
            Y_parsed,
            U_parsed,
            H,
            Rinv,
            lambda_vec,
            N,
            n,
            q,
            L,
            is_var_len,
            horizons,
        )
        cost_history.append(J)
        n_iter += 1

        if n_iter >= 2:
            J_prev = cost_history[-2]
            rel_change = abs(J - J_prev) / max(abs(J_prev), 1.0)

            if rel_change < tolerance:
                if do_trust_region and mu_current > mu_tol:
                    if J <= J_converged_mu:
                        J_converged_mu = J
                        X_best = X_hat
                        A_best = A.copy()
                        B_best = B.copy()
                        mu_current = mu_current / 2
                    else:
                        X_hat = X_best
                        A = A_best.copy()
                        B = B_best.copy()
                        mu_current = 0
                elif do_trust_region and mu_current > 0:
                    J_converged_mu = J
                    X_best = X_hat
                    A_best = A.copy()
                    B_best = B.copy()
                    mu_current = 0
                else:
                    break

    # ------------------------------------------------------------------
    # 5. Add uncertainty from final COSMIC step (SPEC.md S8.12.9)
    # ------------------------------------------------------------------
    A_std, B_std, P, Sigma, noise_var, dof = _add_uncertainty(
        S_c, D_c, Xl_c, C_c, lambda_vec, N, n, q, cov_mode
    )

    return LTVIOResult(
        a=A,
        b=B,
        x=X_hat,
        h=H.copy(),
        r=R_mat.copy(),
        cost=np.array(cost_history),
        iterations=n_iter,
        lambda_=lambda_vec,
        data_length=N,
        state_dim=n,
        output_dim=py,
        input_dim=q,
        num_trajectories=L,
        a_std=A_std,
        b_std=B_std,
        p_cov=P,
        noise_cov=Sigma,
        noise_cov_estimated=True,
        noise_variance=noise_var,
        degrees_of_freedom=dof,
        algorithm="cosmic",
        method="ltv_disc_io",
    )


# ======================================================================
# Private helpers
# ======================================================================


def _parse_inputs(
    Y: np.ndarray | list,
    U: np.ndarray | list,
    H: np.ndarray,
    *,
    lambda_: float | np.ndarray,
    R: np.ndarray | None,
    max_iter: int,
    tolerance: float,
    covariance_mode: str,
    trust_region: float | str,
    trust_region_tol: float,
) -> tuple:
    """Validate and parse all inputs for :func:`ltv_disc_io`.

    Returns
    -------
    tuple
        ``(Y, U, H, lambda_vec, R, do_trust_region, mu, mu_tol,
        cov_mode, N, n, py, q, L, is_var_len, horizons)``
    """
    H = np.asarray(H, dtype=np.float64)
    if H.ndim != 2:
        raise SidError("dim_mismatch", "H must be a 2-D matrix (py x n).")
    py = H.shape[0]
    n = H.shape[1]

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

        Y_list = [np.asarray(y, dtype=np.float64) for y in Y]
        U_list = [np.asarray(u, dtype=np.float64) for u in U]

        for i in range(L):
            if Y_list[i].ndim == 1:
                Y_list[i] = Y_list[i][:, np.newaxis]
            if U_list[i].ndim == 1:
                U_list[i] = U_list[i][:, np.newaxis]

        q = U_list[0].shape[1]
        horizons = np.empty(L, dtype=np.intp)

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
            Nl = U_list[i].shape[0]
            if Y_list[i].shape[0] != Nl + 1:
                raise SidError(
                    "dim_mismatch",
                    f"Y[{i}] has {Y_list[i].shape[0]} rows but U[{i}] "
                    f"has {Nl} (need N_l+1 and N_l).",
                )
            if Nl < 2:
                raise SidError(
                    "too_short",
                    f"Trajectory {i} has fewer than 3 measurements.",
                )
            horizons[i] = Nl

        N = int(np.max(horizons))
        Y_out: np.ndarray | list = Y_list
        U_out: np.ndarray | list = U_list

    else:
        # -- Uniform-horizon mode --
        horizons = None

        Y_arr = np.asarray(Y, dtype=np.float64)
        U_arr = np.asarray(U, dtype=np.float64)

        if Y_arr.ndim == 2:
            Y_arr = Y_arr[:, :, np.newaxis]
        if U_arr.ndim == 2:
            U_arr = U_arr[:, :, np.newaxis]

        N = U_arr.shape[0]
        q = U_arr.shape[1]
        L = Y_arr.shape[2]

        if Y_arr.shape[0] != N + 1:
            raise SidError(
                "dim_mismatch",
                f"Y must have N+1={N + 1} rows, got {Y_arr.shape[0]}.",
            )
        if Y_arr.shape[1] != py:
            raise SidError(
                "dim_mismatch",
                f"Y has {Y_arr.shape[1]} columns but H has {py} rows.",
            )
        if U_arr.shape[2] != L:
            raise SidError(
                "dim_mismatch",
                f"U has {U_arr.shape[2]} trajectories but Y has {L}.",
            )

        Y_out = Y_arr
        U_out = U_arr

    # -- Validate lambda --
    lam = np.asarray(lambda_, dtype=np.float64)
    if lam.ndim == 0:
        # Scalar -> expand to (N-1,) vector
        if lam <= 0:
            raise SidError("bad_lambda", "Lambda must be positive.")
        lambda_vec = np.full(N - 1, float(lam))
    else:
        lambda_vec = lam.ravel()
        if lambda_vec.size != N - 1:
            raise SidError(
                "bad_lambda",
                f"Lambda vector must have N-1 = {N - 1} elements, got {lambda_vec.size}.",
            )
        if np.any(lambda_vec <= 0):
            raise SidError(
                "bad_lambda",
                "All lambda values must be positive.",
            )

    # -- Validate R --
    if R is None:
        R_mat = np.eye(py)
    else:
        R_mat = np.asarray(R, dtype=np.float64)
        if R_mat.shape != (py, py):
            raise SidError(
                "dim_mismatch",
                f"R must be ({py}, {py}), got {R_mat.shape}.",
            )

    # -- Validate covariance_mode --
    cov_mode = covariance_mode.lower()
    if cov_mode not in ("full", "diagonal", "isotropic"):
        raise SidError(
            "bad_cov_mode",
            f"covariance_mode must be 'full', 'diagonal', or 'isotropic'. Got '{covariance_mode}'.",
        )

    # -- Trust region --
    if isinstance(trust_region, str) and trust_region.lower() == "off":
        do_trust_region = False
        mu = 1.0
    else:
        do_trust_region = True
        mu = float(trust_region)

    mu_tol = float(trust_region_tol)

    return (
        Y_out,
        U_out,
        H,
        lambda_vec,
        R_mat,
        do_trust_region,
        mu,
        mu_tol,
        cov_mode,
        N,
        n,
        py,
        q,
        L,
        is_var_len,
        horizons,
    )


def _cosmic_step(
    X_hat: np.ndarray | list,
    U: np.ndarray | list,
    lambda_vec: np.ndarray,
    N: int,
    n: int,
    q: int,
    L: int,
    is_var_len: bool,
    horizons: np.ndarray | None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, ...]:
    """Standard COSMIC solve on estimated states.

    Treats X_hat as observed states and solves for
    C(k) = [A(k)'; B(k)'].

    Returns A, B, S, D, Xl, C for optional uncertainty computation.
    """
    if is_var_len:
        D, Xl = build_data_matrices_var_len(X_hat, U, N, n, q, L, horizons)
    else:
        D, Xl = build_data_matrices(X_hat, U, N, n, q, L)

    S, T = build_block_terms(D, Xl, lambda_vec, N, n, q)
    C, _ = cosmic_solve(S, T, lambda_vec, N, n, q)

    # C has shape (n+q, n, N): rows 0..n-1 are A', rows n..n+q-1 are B'
    A = C[:n, :, :].transpose(1, 0, 2)  # (n, n, N)
    B = C[n:, :, :].transpose(1, 0, 2)  # (n, q, N)

    return A, B, S, D, Xl, C


def _evaluate_full_cost(
    X_hat: np.ndarray | list,
    A: np.ndarray,
    B: np.ndarray,
    Y: np.ndarray | list,
    U: np.ndarray | list,
    H: np.ndarray,
    Rinv: np.ndarray,
    lambda_vec: np.ndarray,
    N: int,
    n: int,
    q: int,
    L: int,
    is_var_len: bool,
    horizons: np.ndarray | None,
) -> float:
    """Compute full Output-COSMIC objective.

    J = obs_fidelity + dyn_fidelity + smoothness
    """
    obs_fidelity = 0.0
    dyn_fidelity = 0.0
    smoothness = 0.0

    for traj in range(L):
        if is_var_len:
            Nl = int(horizons[traj])
            Yl = Y[traj]
            Xl = X_hat[traj]
            Ul = U[traj]
        else:
            Nl = N
            Yl = Y[:, :, traj]
            Xl = X_hat[:, :, traj]
            Ul = U[:, :, traj]

        for k in range(Nl + 1):
            res_obs = Yl[k, :] - H @ Xl[k, :]
            obs_fidelity += res_obs @ Rinv @ res_obs

        for k in range(Nl):
            res_dyn = Xl[k + 1, :] - A[:, :, k] @ Xl[k, :] - B[:, :, k] @ Ul[k, :]
            dyn_fidelity += res_dyn @ res_dyn

    # Smoothness: lambda(k) * ||C(k+1) - C(k)||^2_F
    for k in range(N - 1):
        Ck = np.vstack([A[:, :, k].T, B[:, :, k].T])
        Ck1 = np.vstack([A[:, :, k + 1].T, B[:, :, k + 1].T])
        smoothness += lambda_vec[k] * np.sum((Ck1 - Ck) ** 2)

    return obs_fidelity + dyn_fidelity + smoothness


def _add_uncertainty(
    S: np.ndarray | None,
    D: np.ndarray | list | None,
    Xl: np.ndarray | list | None,
    C: np.ndarray | None,
    lambda_vec: np.ndarray,
    N: int,
    n: int,
    q: int,
    cov_mode: str,
) -> tuple[
    np.ndarray | None,
    np.ndarray | None,
    np.ndarray | None,
    np.ndarray | None,
    float | None,
    float | None,
]:
    """Compute Bayesian uncertainty from the final COSMIC step.

    Returns (A_std, B_std, P, Sigma, noise_var, dof).
    """
    if S is None or D is None or Xl is None or C is None:
        return None, None, None, None, None, None

    d = n + q

    # Diagonal blocks of the Hessian inverse
    P = uncertainty_backward_pass(S, lambda_vec, N, d)

    # Noise covariance from COSMIC residuals
    Sigma, dof = estimate_noise_cov(C, D, Xl, P, cov_mode, N, n, q)

    # Standard deviations of A(k) and B(k) entries
    A_std, B_std = extract_std(P, Sigma, N, n, q)

    noise_var = float(np.trace(Sigma) / n)

    return A_std, B_std, P, Sigma, noise_var, dof
