# Copyright (c) 2026 Pedro Lourenco. All rights reserved.
# This code is released under the MIT License. See LICENSE file in the
# project root for full license information.
#
# This module is part of the Open Source System Identification Toolbox (SID).
# https://github.com/pdlourenco/sid

"""Lambda tuning for LTV state-space identification (validation or frequency-based)."""

from __future__ import annotations

import warnings

import numpy as np

from sid._exceptions import SidError
from sid.ltv_disc import ltv_disc
from sid.ltv_disc_frozen import ltv_disc_frozen
from sid.freq_map import freq_map


def ltv_disc_tune(
    *args,
    method: str = "validation",
    lambda_grid: np.ndarray | None = None,
    precondition: bool = False,
    algorithm: str = "cosmic",
    # Frequency method options
    segment_length: int | None = None,
    consistency_threshold: float = 0.90,
    coherence_threshold: float = 0.3,
) -> tuple:
    """Tune the regularization parameter lambda for :func:`sid.ltv_disc`.

    Two tuning strategies are available:

    * **Validation** (default) -- grid search over lambda evaluated by
      trajectory prediction RMSE on held-out state data.
    * **Frequency** -- compare the COSMIC frozen transfer function against
      a non-parametric :func:`sid.freq_map` estimate using a
      Mahalanobis-like consistency score.  Selects the largest lambda
      whose posterior bands are consistent with the non-parametric bands.

    Parameters
    ----------
    *args
        Positional data arguments.  Interpretation depends on *method*:

        **Validation method** (4 positional args):

        - ``X_train`` -- Training state data, shape ``(N+1, p, L_train)``.
        - ``U_train`` -- Training input data, shape ``(N, q, L_train)``.
        - ``X_val`` -- Validation state data, shape ``(N+1, p, L_val)``.
        - ``U_val`` -- Validation input data, shape ``(N, q, L_val)``.

        **Frequency method** (2 positional args):

        - ``X`` -- State data, shape ``(N+1, p, L)`` or list of arrays.
        - ``U`` -- Input data, shape ``(N, q, L)`` or list of arrays.
    method : str, optional
        ``'validation'`` (default) or ``'frequency'``.
    lambda_grid : ndarray or None, optional
        Vector of candidate lambda values.  Defaults:

        - Validation: ``logspace(-3, 15, 50)``
        - Frequency: ``logspace(0, 10, 25)``
    precondition : bool, optional
        Passed through to :func:`sid.ltv_disc`.  Default is ``False``.
    algorithm : str, optional
        Passed through to :func:`sid.ltv_disc`.  Default is ``'cosmic'``.
    segment_length : int or None, optional
        *(Frequency only)* Outer segment length for :func:`sid.freq_map`.
        Default: ``min(N // 4, 256)``.
    consistency_threshold : float, optional
        *(Frequency only)* Fraction of ``(omega, t)`` grid points
        required to be consistent.  Default is ``0.90``.
    coherence_threshold : float, optional
        *(Frequency only)* Minimum coherence for a grid point to be
        included in the consistency test.  Default is ``0.3``.

    Returns
    -------
    tuple
        Return type depends on *method*:

        **Validation method** -- ``(best_result, best_lambda, all_losses)``

        - *best_result* (:class:`~sid._results.LTVResult`) -- Result at
          optimal lambda.
        - *best_lambda* (*float*) -- Optimal scalar lambda.
        - *all_losses* (*ndarray, shape (n_grid,)*) -- Trajectory RMSE
          at each grid lambda.

        **Frequency method** -- ``(best_result, best_lambda, info_dict)``

        - *best_result* (:class:`~sid._results.LTVResult`) -- Result at
          optimal lambda (with uncertainty).
        - *best_lambda* (*float*) -- Optimal scalar lambda.
        - *info_dict* (*dict*) -- Dictionary with keys:

          - ``'lambda_grid'`` -- sorted lambda grid, shape ``(n_grid,)``.
          - ``'fractions'`` -- consistency fraction per lambda, shape
            ``(n_grid,)``.
          - ``'best_fraction'`` -- consistency fraction at selected lambda.
          - ``'freq_map_results'`` -- list of
            :class:`~sid._results.FreqMapResult`, one per state channel.
          - ``'chi2_threshold'`` -- chi-square threshold used (5.991).

    Raises
    ------
    SidError
        If *method* is not ``'validation'`` or ``'frequency'``
        (code: ``'bad_method'``).
    SidError
        If the number of positional arguments is wrong for the chosen
        method (code: ``'bad_args'``).

    Examples
    --------
    Validation-based tuning:

    >>> import numpy as np
    >>> import sid  # doctest: +SKIP
    >>> best, lam, losses = sid.ltv_disc_tune(  # doctest: +SKIP
    ...     X_train, U_train, X_val, U_val)

    Frequency-based tuning (no validation data needed):

    >>> best, lam, info = sid.ltv_disc_tune(  # doctest: +SKIP
    ...     X, U, method='frequency')

    Custom lambda grid:

    >>> grid = np.logspace(-1, 12, 30)  # doctest: +SKIP
    >>> best, lam, losses = sid.ltv_disc_tune(  # doctest: +SKIP
    ...     X_train, U_train, X_val, U_val, lambda_grid=grid)

    Notes
    -----
    **Validation method algorithm:**

    1. For each lambda in the grid, run :func:`sid.ltv_disc` on training
       data and compute the trajectory prediction RMSE on validation data.
    2. Select the lambda with minimum RMSE.
    3. Re-run :func:`sid.ltv_disc` at the optimal lambda.

    **Frequency method algorithm (SPEC.md S8.11):**

    1. Run :func:`sid.freq_map` per state component (SISO) to obtain
       non-parametric frequency response estimates with uncertainty.
    2. For each lambda, run :func:`sid.ltv_disc` with uncertainty and
       compute the frozen transfer function via
       :func:`sid.ltv_disc_frozen`.
    3. Compute a Mahalanobis-like distance at each ``(omega, t)`` grid
       point: ``d^2 = |G_frozen - G_data|^2 / (std_frozen^2 + std_data^2)``.
    4. Count the fraction of grid points (above the coherence threshold)
       where ``d^2 < chi2_threshold`` (5.991 for 95% confidence, 2 DOF).
    5. Select the largest lambda achieving >= *consistency_threshold*
       fraction.  Falls back to the lambda with the best fraction if
       none meet the threshold.

    **Specification:** SPEC.md S8.4, S8.11

    References
    ----------
    .. [1] Carvalho, Soares, Lourenco, Ventura. "COSMIC: fast closed-form
       identification from large-scale data for LTV systems."
       arXiv:2112.04355, 2022.
    .. [2] Ljung, L. "System Identification", 2nd ed., Prentice Hall, 1999.

    See Also
    --------
    sid.ltv_disc : LTV state-space identification.
    sid.ltv_disc_frozen : Frozen transfer function from LTV model.
    sid.freq_map : Time-varying frequency response estimation.

    Changelog
    ---------
    2026-04-08 : First version (Python port) by Pedro Lourenco.
    """
    m = method.lower()
    if m == "validation":
        return _validation_tune(
            args,
            lambda_grid=lambda_grid,
            precondition=precondition,
            algorithm=algorithm,
        )
    elif m == "frequency":
        return _frequency_tune(
            args,
            lambda_grid=lambda_grid,
            precondition=precondition,
            algorithm=algorithm,
            segment_length=segment_length,
            consistency_threshold=consistency_threshold,
            coherence_threshold=coherence_threshold,
        )
    else:
        raise SidError(
            "bad_method",
            f"Method must be 'validation' or 'frequency'. Got '{method}'.",
        )


# ======================================================================
#  Validation-based tuning
# ======================================================================


def _validation_tune(
    positional_args: tuple,
    *,
    lambda_grid: np.ndarray | None,
    precondition: bool,
    algorithm: str,
) -> tuple:
    """Grid search over lambda, evaluated by trajectory RMSE on validation data."""
    if len(positional_args) != 4:
        raise SidError(
            "bad_args",
            "Validation method requires 4 positional args: "
            "X_train, U_train, X_val, U_val. "
            f"Got {len(positional_args)}.",
        )

    X_train, U_train, X_val, U_val = positional_args

    # ---- Default grid ----
    if lambda_grid is None:
        grid = np.logspace(-3, 15, 50)
    else:
        grid = np.asarray(lambda_grid, dtype=np.float64).ravel()

    n_grid = len(grid)

    # ---- Build extra kwargs for ltv_disc ----
    extra_kwargs: dict = {}
    if precondition:
        extra_kwargs["precondition"] = precondition
    if algorithm != "cosmic":
        extra_kwargs["algorithm"] = algorithm

    # ---- Dimensions from validation data ----
    X_val = np.asarray(X_val, dtype=np.float64)
    U_val = np.asarray(U_val, dtype=np.float64)
    if X_val.ndim == 2:
        X_val = X_val[:, :, np.newaxis]
    if U_val.ndim == 2:
        U_val = U_val[:, :, np.newaxis]

    N = X_val.shape[0] - 1
    p = X_val.shape[1]
    L_val = X_val.shape[2]

    # ---- Grid search ----
    all_losses = np.zeros(n_grid)

    for j in range(n_grid):
        res = ltv_disc(X_train, U_train, lambda_=grid[j], **extra_kwargs)
        all_losses[j] = _trajectory_rmse(res.a, res.b, X_val, U_val, N, p, L_val)

    # ---- Select best ----
    best_idx = int(np.argmin(all_losses))
    best_lambda = float(grid[best_idx])
    best_result = ltv_disc(X_train, U_train, lambda_=best_lambda, **extra_kwargs)

    return (best_result, best_lambda, all_losses)


# ======================================================================
#  Frequency-response consistency tuning
# ======================================================================


def _frequency_tune(
    positional_args: tuple,
    *,
    lambda_grid: np.ndarray | None,
    precondition: bool,
    algorithm: str,
    segment_length: int | None,
    consistency_threshold: float,
    coherence_threshold: float,
) -> tuple:
    """Select lambda via frequency-response consistency (SPEC.md S8.11)."""
    if len(positional_args) != 2:
        raise SidError(
            "bad_args",
            f"Frequency method requires 2 positional args: X, U. Got {len(positional_args)}.",
        )

    X, U = positional_args

    # ---- Default grid ----
    if lambda_grid is None:
        grid = np.logspace(0, 10, 25)
    else:
        grid = np.asarray(lambda_grid, dtype=np.float64).ravel()

    grid = np.sort(grid)
    n_grid = len(grid)

    # ---- Build extra kwargs for ltv_disc ----
    extra_kwargs: dict = {}
    if precondition:
        extra_kwargs["precondition"] = precondition
    if algorithm != "cosmic":
        extra_kwargs["algorithm"] = algorithm

    # ---- Dimensions ----
    if isinstance(X, list):
        N = X[0].shape[0] - 1
        p = X[0].shape[1] if np.asarray(X[0]).ndim > 1 else 1
    else:
        X_arr = np.asarray(X)
        N = X_arr.shape[0] - 1
        p = X_arr.shape[1] if X_arr.ndim > 1 else 1

    seg_len = segment_length if segment_length is not None else min(N // 4, 256)

    # ---- Step 1: Run freq_map per state component (SISO) ----
    # MIMO mode produces NaN for ResponseStd in v1.0, so SISO per-channel
    # is needed to get valid uncertainty estimates for the Mahalanobis test.
    fmap_results = []
    for ch in range(p):
        if isinstance(X, list):
            # Build list of single-channel outputs from all trajectories
            y_ch = [np.asarray(xl)[:-1, ch] for xl in X]
        else:
            X_arr = np.asarray(X)
            if X_arr.ndim == 3:
                n_traj = X_arr.shape[2]
                if n_traj > 1:
                    y_ch = X_arr[:N, ch : ch + 1, :]  # (N, 1, L) — single channel, L trajectories
                else:
                    y_ch = X_arr[:N, ch, 0]  # single traj
            elif X_arr.ndim == 2:
                y_ch = X_arr[:N, ch]
            else:
                y_ch = X_arr[:N]

        u_freq = U
        fmap_results.append(freq_map(y_ch, u_freq, segment_length=seg_len))

    fmap_freqs = fmap_results[0].frequency
    Ts = fmap_results[0].sample_time

    # ---- Align time grids ----
    seg_center_samples = fmap_results[0].time / Ts
    k_nearest = np.clip(np.round(seg_center_samples), 0, N - 1).astype(int)
    nk = len(k_nearest)
    nf = len(fmap_freqs)

    # Chi-square threshold for 95% confidence, 2 DOF (complex scalar SISO)
    chi2_threshold = 5.991

    # ---- Step 2: Grid search with Mahalanobis scoring ----
    fractions = np.zeros(n_grid)

    for j in range(n_grid):
        # Run COSMIC with uncertainty
        res = ltv_disc(X, U, lambda_=grid[j], uncertainty=True, **extra_kwargs)

        # Frozen transfer function at aligned time steps and matching frequencies
        frz = ltv_disc_frozen(res, frequencies=fmap_freqs, time_steps=k_nearest)

        # Average per-channel consistency fraction across state components
        chan_fracs = np.zeros(p)
        for ch in range(p):
            # Extract SISO frozen TF for channel ch: response shape (nf, p, q, nk)
            G_frz_ch = frz.response[:, ch, :, :].reshape(nf, nk)
            if frz.response_std is not None:
                GStd_frz_ch = frz.response_std[:, ch, :, :].reshape(nf, nk)
            else:
                GStd_frz_ch = np.zeros((nf, nk))

            # freq_map SISO data for this channel: response shape (nf, K)
            G_dat_ch = fmap_results[ch].response
            GStd_dat_ch = fmap_results[ch].response_std

            if fmap_results[ch].coherence is not None:
                coh_mask = fmap_results[ch].coherence >= coherence_threshold
            else:
                coh_mask = np.ones((nf, nk), dtype=bool)

            chan_fracs[ch] = _compute_consistency_siso(
                G_dat_ch,
                GStd_dat_ch,
                G_frz_ch,
                GStd_frz_ch,
                coh_mask,
                chi2_threshold,
            )
        fractions[j] = np.mean(chan_fracs)

    # ---- Step 3: Select largest lambda with sufficient consistency ----
    consistent = np.where(fractions >= consistency_threshold)[0]
    if len(consistent) > 0:
        best_idx = int(consistent[-1])  # largest lambda (grid sorted ascending)
        best_lambda = float(grid[best_idx])
    else:
        # Fallback: no lambda meets threshold, use best available
        best_idx = int(np.argmax(fractions))
        best_lambda = float(grid[best_idx])
        warnings.warn(
            f"No lambda achieved {consistency_threshold * 100:.0f}% consistency. "
            f"Using best ({fractions[best_idx] * 100:.1f}% at lambda={best_lambda:.2e}).",
            stacklevel=2,
        )

    # ---- Re-run at optimal lambda with uncertainty ----
    best_result = ltv_disc(X, U, lambda_=best_lambda, uncertainty=True, **extra_kwargs)

    # ---- Pack info dict ----
    info_dict = {
        "lambda_grid": grid,
        "fractions": fractions,
        "best_fraction": fractions[best_idx],
        "freq_map_results": fmap_results,
        "chi2_threshold": chi2_threshold,
    }

    return (best_result, best_lambda, info_dict)


# ======================================================================
#  Helper functions
# ======================================================================


def _trajectory_rmse(
    A: np.ndarray,
    B: np.ndarray,
    X_val: np.ndarray,
    U_val: np.ndarray,
    N: int,
    p: int,
    L_val: int,
) -> float:
    """Average trajectory prediction RMSE over the validation set.

    Parameters
    ----------
    A : ndarray, shape (p, p, N)
        Time-varying dynamics matrices.
    B : ndarray, shape (p, q, N)
        Time-varying input matrices.
    X_val : ndarray, shape (N+1, p, L_val)
        Validation state trajectories.
    U_val : ndarray, shape (N, q, L_val)
        Validation input trajectories.
    N : int
        Number of time steps.
    p : int
        State dimension.
    L_val : int
        Number of validation trajectories.

    Returns
    -------
    float
        Mean RMSE across validation trajectories.
    """
    total = 0.0
    for traj in range(L_val):
        x_hat = np.zeros((N + 1, p))
        x_hat[0] = X_val[0, :, traj]
        for k in range(N):
            x_hat[k + 1] = A[:, :, k] @ x_hat[k] + B[:, :, k] @ U_val[k, :, traj]
        err = x_hat - X_val[:, :, traj]
        total += np.sqrt(np.mean(np.sum(err**2, axis=1)))
    return total / L_val


def _compute_consistency_siso(
    G_data: np.ndarray,
    GStd_data: np.ndarray,
    G_frozen: np.ndarray,
    GStd_frozen: np.ndarray,
    coh_mask: np.ndarray,
    chi2_threshold: float,
) -> float:
    """Mahalanobis-like consistency fraction for one SISO channel.

    Parameters
    ----------
    G_data : ndarray, shape (nf, K)
        Complex frequency response from :func:`sid.freq_map`.
    GStd_data : ndarray, shape (nf, K)
        Standard deviation of *G_data*.
    G_frozen : ndarray, shape (nf, nk)
        Complex frozen transfer function from :func:`sid.ltv_disc_frozen`.
    GStd_frozen : ndarray, shape (nf, nk)
        Standard deviation of *G_frozen*.
    coh_mask : ndarray, shape (nf, K)
        Boolean mask where coherence is sufficient.
    chi2_threshold : float
        Chi-square threshold (5.991 for 95% confidence, 2 DOF).

    Returns
    -------
    float
        Fraction of valid grid points that are consistent.
    """
    denominator = GStd_frozen**2 + GStd_data**2
    denominator = np.where(
        denominator < np.finfo(np.float64).eps, np.finfo(np.float64).eps, denominator
    )
    d2 = np.abs(G_frozen - G_data) ** 2 / denominator

    is_consistent = (d2 < chi2_threshold) & coh_mask
    is_valid = coh_mask

    n_valid = np.sum(is_valid)
    if n_valid == 0:
        return 0.0
    return float(np.sum(is_consistent) / n_valid)
