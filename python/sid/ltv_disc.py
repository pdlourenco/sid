# Copyright (c) 2026 Pedro Lourenco. All rights reserved.
# This code is released under the MIT License. See LICENSE file in the
# project root for full license information.
#
# This module is part of the Open Source System Identification Toolbox (SID).
# https://github.com/pdlourenco/sid-matlab

"""Discrete-time LTV state-space identification via the COSMIC algorithm."""

from __future__ import annotations

import warnings

import numpy as np

from sid._exceptions import SidError
from sid._internal.ltv_build_block_terms import build_block_terms
from sid._internal.ltv_build_data_matrices import (
    build_data_matrices,
    build_data_matrices_var_len,
)
from sid._internal.ltv_cosmic_solve import cosmic_solve
from sid._internal.ltv_evaluate_cost import evaluate_cost
from sid._internal.ltv_uncertainty_backward_pass import uncertainty_backward_pass
from sid._internal.estimate_noise_cov import estimate_noise_cov
from sid._internal.extract_std import extract_std
from sid._results import LTVResult


def ltv_disc(
    X: np.ndarray | list,
    U: np.ndarray | list,
    *,
    lambda_: float | np.ndarray | str = "auto",
    lambda_grid: np.ndarray | None = None,
    precondition: bool = False,
    algorithm: str = "cosmic",
    uncertainty: bool = False,
    noise_cov: np.ndarray | str = "estimate",
    covariance_mode: str = "diagonal",
) -> LTVResult:
    """Identify a discrete-time LTV state-space model from trajectory data.

    Estimates time-varying system matrices A(k) and B(k) for the model

    .. math::

        x(k+1) = A(k)\\,x(k) + B(k)\\,u(k), \\quad k = 0, \\ldots, N-1

    using the COSMIC algorithm (Carvalho *et al.*, 2022), which solves a
    regularized least-squares problem balancing data fidelity against
    temporal smoothness of the system matrices.  The closed-form
    block-tridiagonal solver has O(N (p+q)^3) complexity.

    This is an open-source replacement for proprietary LTV identification
    routines.

    Parameters
    ----------
    X : ndarray or list of ndarray
        State trajectory data.  Accepted formats:

        - ``(N+1, p)`` -- single trajectory with *p* states
        - ``(N+1, p, L)`` -- *L* trajectories, same horizon *N*
        - ``[X_1, X_2, ...]`` -- list of *L* arrays ``(N_l+1, p)``
          with variable horizons (``N = max(N_l)``)
    U : ndarray or list of ndarray
        Input trajectory data, matching format of *X*:

        - ``(N, q)`` -- single trajectory with *q* inputs
        - ``(N, q, L)`` -- *L* trajectories
        - ``[U_1, U_2, ...]`` -- list of *L* arrays ``(N_l, q)``
    lambda_ : float, ndarray, or ``'auto'``, optional
        Regularization strength.  Options:

        - *scalar* -- uniform regularization at all time steps
        - *ndarray of shape (N-1,)* -- per-step regularization
        - ``'auto'`` -- automatic selection via L-curve (default)

        Must be positive (scalar or all elements).
    lambda_grid : ndarray or None, optional
        Candidate lambda values for L-curve auto-selection.  Only used
        when ``lambda_='auto'``.  Default is ``logspace(-3, 15, 50)``.
    precondition : bool, optional
        Apply block-diagonal preconditioning.  Currently disabled in
        v1.0.  Default is ``False``.
    algorithm : str, optional
        Identification algorithm.  Only ``'cosmic'`` is supported in
        v1.0.  Default is ``'cosmic'``.
    uncertainty : bool, optional
        Compute Bayesian posterior uncertainty for A(k), B(k).  Doubles
        computation cost.  Default is ``False``.
    noise_cov : ndarray or ``'estimate'``, optional
        Measurement noise covariance.  Options:

        - ``(p, p)`` ndarray -- known noise covariance (automatically
          enables uncertainty)
        - ``'estimate'`` -- estimate from residuals (default)
    covariance_mode : str, optional
        How to estimate noise covariance when *noise_cov* is
        ``'estimate'``.  One of ``'diagonal'``, ``'full'``, or
        ``'isotropic'``.  Default is ``'diagonal'``.

    Returns
    -------
    LTVResult
        Frozen dataclass with fields:

        - **a** (*ndarray, shape (p, p, N)*) -- Time-varying dynamics.
        - **b** (*ndarray, shape (p, q, N)*) -- Time-varying input.
        - **a_std** (*ndarray or None*) -- Std dev of ``a`` entries.
        - **b_std** (*ndarray or None*) -- Std dev of ``b`` entries.
        - **p_cov** (*ndarray or None*) -- Row covariance blocks.
        - **noise_cov** (*ndarray or None*) -- Noise covariance.
        - **noise_cov_estimated** (*bool or None*) -- Whether estimated.
        - **noise_variance** (*float or None*) -- ``trace(Sigma)/p``.
        - **degrees_of_freedom** (*float or None*) -- Effective DOF.
        - **lambda_** (*ndarray, shape (N-1,)*) -- Lambda vector used.
        - **cost** (*ndarray, shape (3,)*) -- ``[total, fidelity, reg]``.
        - **data_length** (*int*) -- *N*.
        - **state_dim** (*int*) -- *p*.
        - **input_dim** (*int*) -- *q*.
        - **num_trajectories** (*int*) -- *L*.
        - **algorithm** (*str*) -- ``'cosmic'``.
        - **preconditioned** (*bool*) -- Preconditioning flag.
        - **method** (*str*) -- ``'ltv_disc'``.

    Raises
    ------
    SidError
        If data contains NaN/Inf (code: ``'non_finite'``).
    SidError
        If data dimensions are inconsistent (code: ``'dim_mismatch'``).
    SidError
        If trajectories are too short, N < 2 (code: ``'too_short'``).
    SidError
        If lambda is invalid (code: ``'bad_lambda'``).
    SidError
        If algorithm is unsupported (code: ``'bad_algorithm'``).
    SidError
        If noise_cov is invalid (code: ``'bad_noise_cov'``).
    SidError
        If covariance_mode is invalid (code: ``'bad_cov_mode'``).

    Examples
    --------
    Basic identification with automatic lambda:

    >>> import numpy as np
    >>> import sid  # doctest: +SKIP
    >>> N = 100; p = 2; q = 1
    >>> rng = np.random.default_rng(0)
    >>> X = rng.standard_normal((N + 1, p))
    >>> U = rng.standard_normal((N, q))
    >>> result = sid.ltv_disc(X, U)  # doctest: +SKIP
    >>> result.a.shape  # doctest: +SKIP
    (2, 2, 100)

    Manual uniform lambda:

    >>> result = sid.ltv_disc(X, U, lambda_=1e5)  # doctest: +SKIP

    With uncertainty estimation:

    >>> result = sid.ltv_disc(X, U, lambda_=1e5, uncertainty=True)  # doctest: +SKIP
    >>> result.a_std.shape  # doctest: +SKIP
    (2, 2, 100)

    Variable-length trajectories:

    >>> X_list = [rng.standard_normal((51, 2)), rng.standard_normal((31, 2))]
    >>> U_list = [rng.standard_normal((50, 1)), rng.standard_normal((30, 1))]
    >>> result = sid.ltv_disc(X_list, U_list, lambda_=1e3)  # doctest: +SKIP

    Notes
    -----
    **Algorithm:**

    1. Validate and orient input data; detect variable-length mode.
    2. Build per-step data matrices D(k), X_lead(k) from trajectories
       (SPEC.md S8.3.2).
    3. If ``lambda_='auto'``, run L-curve selection over a grid of
       candidate values.
    4. Build block-tridiagonal terms S(k), Theta(k) (SPEC.md S8.3.3).
    5. COSMIC forward-backward pass: forward pass computes Lbd(k) and
       Y(k), backward pass recovers C(k) = [A(k)'; B(k)']
       (SPEC.md S8.3.4).
    6. Extract A(k), B(k) from C(k).
    7. Evaluate total cost = fidelity + regularization.
    8. If uncertainty requested: backward recursion on P(k), noise
       covariance estimation, and standard deviation extraction
       (SPEC.md S8.9).

    **Specification:** SPEC.md S8 -- Discrete-Time LTV State-Space
    Identification

    References
    ----------
    .. [1] Carvalho, Soares, Lourenco, Ventura. "COSMIC: fast closed-form
       identification from large-scale data for LTV systems."
       arXiv:2112.04355, 2022.

    See Also
    --------
    sid.ltv_disc_tune : Tune lambda via cross-validation.
    sid.ltv_disc_frozen : Identify with frozen (constant) segments.
    sid.freq_map : Time-varying frequency response estimation.

    Changelog
    ---------
    2026-04-08 : First version (Python port) by Pedro Lourenco.
    """
    # ------------------------------------------------------------------
    # 1. Parse and validate inputs
    # ------------------------------------------------------------------
    (
        X_arr,
        U_arr,
        lambda_val,
        do_precondition,
        algo,
        do_uncertainty,
        noise_cov_val,
        cov_mode,
        lam_grid,
        N,
        p,
        q,
        L,
        is_var_len,
        horizons,
    ) = _parse_inputs(
        X,
        U,
        lambda_=lambda_,
        lambda_grid=lambda_grid,
        precondition=precondition,
        algorithm=algorithm,
        uncertainty=uncertainty,
        noise_cov=noise_cov,
        covariance_mode=covariance_mode,
    )

    # ------------------------------------------------------------------
    # 2. Build data matrices (SPEC.md S8.3.2)
    # ------------------------------------------------------------------
    if is_var_len:
        D, Xl = build_data_matrices_var_len(X_arr, U_arr, N, p, q, L, horizons)
    else:
        D, Xl = build_data_matrices(X_arr, U_arr, N, p, q, L)

    # ------------------------------------------------------------------
    # 3. Lambda selection (SPEC.md S8.4)
    # ------------------------------------------------------------------
    if isinstance(lambda_val, str) and lambda_val.lower() == "auto":
        lambda_val = _lcurve_lambda(D, Xl, N, p, q, lam_grid)

    if np.ndim(lambda_val) == 0:
        # Scalar -> expand to (N-1,) vector
        lambda_vec = np.full(N - 1, float(lambda_val))
    else:
        lambda_vec = np.asarray(lambda_val, dtype=np.float64).ravel()

    # ------------------------------------------------------------------
    # 4. Build block tridiagonal terms S_kk, Theta_k (SPEC.md S8.3.3)
    # ------------------------------------------------------------------
    S, T = build_block_terms(D, Xl, lambda_vec, N, p, q)

    # ------------------------------------------------------------------
    # Preconditioning
    # ------------------------------------------------------------------
    precondition_requested = bool(do_precondition)
    if do_precondition:
        warnings.warn(
            "Preconditioning is disabled in v1.0 due to a known issue "
            "with off-diagonal block handling. Using unpreconditioned solver.",
            stacklevel=2,
        )
        do_precondition = False

    # ------------------------------------------------------------------
    # 5. COSMIC forward-backward pass (SPEC.md S8.3.4)
    # ------------------------------------------------------------------
    C, Lbd = cosmic_solve(S, T, lambda_vec, N, p, q)

    # ------------------------------------------------------------------
    # 6. Extract A(k), B(k) from C(k) = [A(k)'; B(k)']
    # ------------------------------------------------------------------
    # C has shape (p+q, p, N): rows 0..p-1 are A', rows p..p+q-1 are B'
    # Transpose first two dims to get (p, p, N) and (p, q, N)
    A = C[:p, :, :].transpose(1, 0, 2)  # (p, p, N)
    B = C[p:, :, :].transpose(1, 0, 2)  # (p, q, N)

    # ------------------------------------------------------------------
    # 7. Evaluate cost
    # ------------------------------------------------------------------
    cost_total, fidelity, reg, _ = evaluate_cost(A, B, D, Xl, lambda_vec, N, p, q)

    # ------------------------------------------------------------------
    # 8. Uncertainty estimation (SPEC.md S8.9)
    # ------------------------------------------------------------------
    A_std = None
    B_std = None
    P = None
    Sigma = None
    noise_cov_estimated = None
    noise_var = None
    dof = None

    if do_uncertainty:
        d = p + q

        # P(k) = diagonal block of H^{-1}
        P = uncertainty_backward_pass(S, lambda_vec, N, d)

        # Noise covariance
        noise_cov_provided = not isinstance(noise_cov_val, str)
        if noise_cov_provided:
            Sigma = noise_cov_val
            dof = np.nan
        else:
            Sigma, dof = estimate_noise_cov(C, D, Xl, P, cov_mode, N, p, q)

        # Standard deviations of A(k), B(k) entries
        A_std, B_std = extract_std(P, Sigma, N, p, q)

        noise_cov_estimated = not noise_cov_provided
        noise_var = float(np.trace(Sigma) / p)

    # ------------------------------------------------------------------
    # 9. Pack result
    # ------------------------------------------------------------------
    return LTVResult(
        a=A,
        b=B,
        a_std=A_std,
        b_std=B_std,
        p_cov=P,
        noise_cov=Sigma,
        noise_cov_estimated=noise_cov_estimated,
        noise_variance=noise_var,
        degrees_of_freedom=dof,
        lambda_=lambda_vec,
        cost=np.array([cost_total, fidelity, reg]),
        data_length=N,
        state_dim=p,
        input_dim=q,
        num_trajectories=L,
        algorithm=algo,
        # v1.0: preconditioning is disabled. Signal this distinctly from
        # "not requested" when the user asked for it, matching MATLAB's
        # sidLTVdisc behaviour (Preconditioned = 'not_implemented').
        preconditioned="not_implemented" if precondition_requested else False,
        method="ltv_disc",
    )


# ======================================================================
# Private helpers
# ======================================================================


def _parse_inputs(
    X: np.ndarray | list,
    U: np.ndarray | list,
    *,
    lambda_: float | np.ndarray | str,
    lambda_grid: np.ndarray | None,
    precondition: bool,
    algorithm: str,
    uncertainty: bool,
    noise_cov: np.ndarray | str,
    covariance_mode: str,
) -> tuple:
    """Validate and parse all inputs for :func:`ltv_disc`.

    Returns
    -------
    tuple
        ``(X, U, lambda_, precondition, algorithm, uncertainty,
        noise_cov, covariance_mode, lambda_grid, N, p, q, L,
        is_var_len, horizons)``
    """
    is_var_len = isinstance(X, list)

    if is_var_len:
        # -- Variable-length trajectory mode (list of arrays) --
        if not isinstance(U, list):
            raise SidError(
                "bad_input",
                "When X is a list, U must also be a list.",
            )
        L = len(X)
        if len(U) != L:
            raise SidError(
                "dim_mismatch",
                f"X has {L} trajectories but U has {len(U)}.",
            )
        if L == 0:
            raise SidError("bad_input", "Trajectory lists must not be empty.")

        # Ensure each element is an ndarray
        X = [np.asarray(x, dtype=np.float64) for x in X]
        U = [np.asarray(u, dtype=np.float64) for u in U]

        # Handle 1-D state arrays: reshape (N+1,) -> (N+1, 1)
        for i in range(L):
            if X[i].ndim == 1:
                X[i] = X[i][:, np.newaxis]
            if U[i].ndim == 1:
                U[i] = U[i][:, np.newaxis]

        p = X[0].shape[1]
        q = U[0].shape[1]

        horizons = np.empty(L, dtype=np.intp)
        for i in range(L):
            if X[i].shape[1] != p:
                raise SidError(
                    "dim_mismatch",
                    f"Trajectory {i} has {X[i].shape[1]} state dims, expected {p}.",
                )
            if U[i].shape[1] != q:
                raise SidError(
                    "dim_mismatch",
                    f"Trajectory {i} has {U[i].shape[1]} input dims, expected {q}.",
                )
            Nl = X[i].shape[0] - 1
            if U[i].shape[0] != Nl:
                raise SidError(
                    "dim_mismatch",
                    f"Trajectory {i}: U has {U[i].shape[0]} rows but X "
                    f"has {X[i].shape[0]} (need N and N+1).",
                )
            if Nl < 1:
                raise SidError(
                    "too_short",
                    f"Trajectory {i} has fewer than 2 state measurements.",
                )
            horizons[i] = Nl

            # Check for NaN/Inf
            if not np.all(np.isfinite(X[i])):
                raise SidError(
                    "non_finite",
                    f"State data X[{i}] contains NaN or Inf.",
                )
            if not np.all(np.isfinite(U[i])):
                raise SidError(
                    "non_finite",
                    f"Input data U[{i}] contains NaN or Inf.",
                )

        N = int(np.max(horizons))
        if N < 2:
            raise SidError(
                "too_short",
                "Need at least 3 state measurements (N >= 2).",
            )

    else:
        # -- Uniform-horizon mode (ndarray) --
        horizons = None

        X = np.asarray(X, dtype=np.float64)
        U = np.asarray(U, dtype=np.float64)

        # Ensure 3D: (N+1, p, L)
        if X.ndim == 1:
            X = X[:, np.newaxis, np.newaxis]
        elif X.ndim == 2:
            X = X[:, :, np.newaxis]

        if U.ndim == 1:
            U = U[:, np.newaxis, np.newaxis]
        elif U.ndim == 2:
            U = U[:, :, np.newaxis]

        N = X.shape[0] - 1
        p = X.shape[1]
        q = U.shape[1]
        L = X.shape[2]

        # Validate dimensions
        if U.shape[0] != N:
            raise SidError(
                "dim_mismatch",
                f"U must have {N} rows (N), but has {U.shape[0]}. X has N+1 = {N + 1} rows.",
            )
        if U.shape[2] != L:
            raise SidError(
                "dim_mismatch",
                f"X has {L} trajectories but U has {U.shape[2]}.",
            )
        if N < 2:
            raise SidError(
                "too_short",
                "Need at least 3 state measurements (N >= 2).",
            )

        # Check for NaN/Inf
        if not np.all(np.isfinite(X)):
            raise SidError("non_finite", "State data X contains NaN or Inf.")
        if not np.all(np.isfinite(U)):
            raise SidError("non_finite", "Input data U contains NaN or Inf.")

    # -- Validate algorithm --
    algo = algorithm.lower()
    if algo != "cosmic":
        raise SidError(
            "bad_algorithm",
            f"Only 'cosmic' is supported in v1.0. Got '{algorithm}'.",
        )

    # -- If noise_cov is a matrix, enable uncertainty automatically --
    if not isinstance(noise_cov, str):
        noise_cov = np.asarray(noise_cov, dtype=np.float64)
        uncertainty = True

    # -- Validate lambda --
    if not isinstance(lambda_, str):
        lam = np.asarray(lambda_, dtype=np.float64)
        if lam.ndim == 0:
            # Scalar
            if lam <= 0:
                raise SidError("bad_lambda", "Lambda must be positive.")
        else:
            lam = lam.ravel()
            if lam.size != N - 1:
                raise SidError(
                    "bad_lambda",
                    f"Lambda vector must have N-1 = {N - 1} elements, got {lam.size}.",
                )
            if np.any(lam <= 0):
                raise SidError(
                    "bad_lambda",
                    "All lambda values must be positive.",
                )
        lambda_ = lam

    # -- Validate noise_cov --
    if not isinstance(noise_cov, str):
        if noise_cov.ndim != 2:
            raise SidError(
                "bad_noise_cov",
                "NoiseCov must be a p x p matrix or 'estimate'.",
            )
        if noise_cov.shape[0] != p or noise_cov.shape[1] != p:
            raise SidError(
                "bad_noise_cov",
                f"NoiseCov must be {p} x {p} (matching state dimension p), "
                f"got {noise_cov.shape[0]} x {noise_cov.shape[1]}.",
            )
        if not np.all(np.isfinite(noise_cov)):
            raise SidError(
                "bad_noise_cov",
                "NoiseCov contains NaN or Inf.",
            )

    # -- Validate covariance_mode --
    cov_mode = covariance_mode.lower()
    if cov_mode not in ("full", "diagonal", "isotropic"):
        raise SidError(
            "bad_cov_mode",
            f"covariance_mode must be 'full', 'diagonal', or 'isotropic'. Got '{covariance_mode}'.",
        )

    return (
        X,
        U,
        lambda_,
        precondition,
        algo,
        uncertainty,
        noise_cov,
        cov_mode,
        lambda_grid,
        N,
        p,
        q,
        L,
        is_var_len,
        horizons,
    )


def _lcurve_lambda(
    D: np.ndarray,
    Xl: np.ndarray,
    N: int,
    p: int,
    q: int,
    lambda_grid: np.ndarray | None,
) -> float:
    """Select regularization parameter via L-curve maximum curvature.

    Parameters
    ----------
    D : ndarray
        Stacked data matrices from :func:`build_data_matrices`.
    Xl : ndarray
        Stacked leading-state matrices.
    N : int
        Number of time steps.
    p : int
        State dimension.
    q : int
        Input dimension.
    lambda_grid : ndarray or None
        Candidate grid; defaults to ``logspace(-3, 15, 50)``.

    Returns
    -------
    float
        Selected scalar lambda value.
    """
    if lambda_grid is None:
        grid = np.logspace(-3, 15, 50)
    else:
        grid = np.asarray(lambda_grid, dtype=np.float64).ravel()

    n_grid = len(grid)
    F = np.zeros(n_grid)  # data fidelity
    R = np.zeros(n_grid)  # unweighted variation (for L-curve)

    for j in range(n_grid):
        lam = np.full(N - 1, grid[j])
        S, T = build_block_terms(D, Xl, lam, N, p, q)
        C, _ = cosmic_solve(S, T, lam, N, p, q)

        A = C[:p, :, :].transpose(1, 0, 2)
        B = C[p:, :, :].transpose(1, 0, 2)
        _, F[j], _, R[j] = evaluate_cost(A, B, D, Xl, lam, N, p, q)

    # L-curve: find corner of maximum curvature in log-log space
    _eps = np.finfo(np.float64).tiny
    lf = np.log10(np.maximum(F, _eps))
    lr = np.log10(np.maximum(R, _eps))

    # Curvature via finite differences
    kappa = np.zeros(n_grid)
    for j in range(1, n_grid - 1):
        df1 = lf[j] - lf[j - 1]
        df2 = lf[j + 1] - lf[j]
        dr1 = lr[j] - lr[j - 1]
        dr2 = lr[j + 1] - lr[j]

        ddf = df2 - df1
        ddr = dr2 - dr1

        num = abs(ddf * (dr1 + dr2) / 2 - ddr * (df1 + df2) / 2)
        den = ((df1 + df2) ** 2 / 4 + (dr1 + dr2) ** 2 / 4) ** 1.5

        if den > 0:
            kappa[j] = num / den

    idx = int(np.argmax(kappa))
    return float(grid[idx])
