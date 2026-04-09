# Copyright (c) 2026 Pedro Lourenco. All rights reserved.
# This code is released under the MIT License. See LICENSE file in the
# project root for full license information.
#
# This module is part of the Open Source System Identification Toolbox (SID).
# https://github.com/pdlourenco/sid-matlab

"""Model residual analysis and diagnostic tests."""

from __future__ import annotations

import numpy as np

from sid._exceptions import SidError
from sid._internal.cov import sid_cov
from sid._internal.freq_domain_sim import freq_domain_sim
from sid._results import ResidualResult


def residual(
    model: object,
    y: np.ndarray,
    u: np.ndarray | None = None,
    *,
    max_lag: int | None = None,
    plot: bool = False,
) -> ResidualResult:
    """Compute model residuals and perform diagnostic tests.

    This is the Python port of ``sidResidual.m``.

    Computes residuals from an estimated model and performs whiteness
    (autocorrelation) and independence (cross-correlation) tests to
    assess model quality.

    Parameters
    ----------
    model : object
        Result from any ``sid`` estimator.  Must have either ``a`` and
        ``b`` attributes (state-space / COSMIC) or a ``response``
        attribute (frequency-domain).
    y : ndarray, shape ``(N, ny)`` or ``(N+1, p, L)``
        Measured output data.  For state-space models the array contains
        state trajectories ``(N+1, p)`` or ``(N+1, p, L)``.
    u : ndarray or None, optional
        Input data, shape ``(N, nu)`` or ``(N, q, L)``.  ``None`` for
        time-series models (no input).
    max_lag : int or None, optional
        Maximum lag for correlation tests.
        Default: ``min(25, N // 5)``.
    plot : bool, optional
        If ``True``, display a diagnostic plot (requires matplotlib).
        Default: ``False``.

    Returns
    -------
    dict
        Dictionary with the following keys:

        - ``'residual'`` -- ``(N, ny)`` residual time series.
        - ``'auto_corr'`` -- ``(max_lag+1,)`` normalised autocorrelation
          of the first channel (for quick inspection).
        - ``'auto_corr_all'`` -- ``(max_lag+1, ny)`` per-channel
          autocorrelation.
        - ``'cross_corr'`` -- ``(2*max_lag+1, ny*nu)`` normalised
          cross-correlation, or empty array for time-series.
        - ``'confidence_bound'`` -- scalar 99% bound ``2.58/sqrt(N)``.
        - ``'whiteness_pass'`` -- bool, ``True`` if all channels pass.
        - ``'whiteness_pass_all'`` -- ``(ny,)`` bool array, per-channel.
        - ``'independence_pass'`` -- bool, ``True`` if all pairs pass
          (always ``True`` for time-series).
        - ``'independence_pass_all'`` -- ``(ny*nu,)`` bool array
          (omitted for time-series).
        - ``'data_length'`` -- int, effective number of samples *N*.

    Raises
    ------
    SidError
        If the model type cannot be determined (code: ``'bad_model'``).

    Examples
    --------
    >>> import sid  # doctest: +SKIP
    >>> G = sid.freq_bt(y, u)  # doctest: +SKIP
    >>> res = sid.residual(G, y, u)  # doctest: +SKIP

    Notes
    -----
    **Specification:** (Model residual analysis -- not yet in SPEC.md)

    The whiteness test checks that the normalised autocorrelation
    ``|r_ee(tau)| < 2.58 / sqrt(N)`` for ``tau = 1, ..., max_lag``.
    The independence test checks the same bound for the normalised
    cross-correlation between residuals and inputs at all lags in
    ``[-max_lag, max_lag]``.  Both tests use a 99% confidence level.

    See Also
    --------
    sid.compare : Model output comparison.
    sid._internal.cov.sid_cov : Biased cross-covariance estimator.

    Changelog
    ---------
    2026-04-09 : First version (Python port) by Pedro Lourenco.
    """
    is_time_series: bool = u is None

    # ------------------------------------------------------------------
    # Dispatch on model type
    # ------------------------------------------------------------------
    if hasattr(model, "a") and hasattr(model, "b"):
        e, N_eff = _compute_residual_ss(model, y, u)
    elif hasattr(model, "response"):
        e, N_eff = _compute_residual_freq(model, y, u)
    else:
        raise SidError(
            "bad_model",
            "Model must have a 'response' attribute (freq-domain) "
            "or 'a'/'b' attributes (state-space).",
        )

    # ------------------------------------------------------------------
    # Default max_lag
    # ------------------------------------------------------------------
    if max_lag is None:
        max_lag = min(25, N_eff // 5)

    # ------------------------------------------------------------------
    # Per-channel whiteness and independence tests
    # ------------------------------------------------------------------
    ny: int = e.shape[1]
    conf_bound: float = 2.58 / np.sqrt(N_eff)

    # Prepare 2-D input for cross-correlation
    if not is_time_series:
        u_arr = np.asarray(u)
        if u_arr.ndim == 3:
            nu: int = u_arr.shape[1]
            u_2d = u_arr[:, :, 0]  # first trajectory
        else:
            if u_arr.ndim == 1:
                u_2d = u_arr[:, np.newaxis]
            else:
                u_2d = u_arr
            nu = u_2d.shape[1]
        N_use: int = min(e.shape[0], u_2d.shape[0])
    else:
        nu = 0
        u_2d = None
        N_use = e.shape[0]

    # Storage
    auto_corr_all = np.zeros((max_lag + 1, ny))
    whiteness_pass_all = np.ones(ny, dtype=bool)

    if not is_time_series:
        n_pairs: int = ny * nu
        cross_corr_all = np.zeros((2 * max_lag + 1, n_pairs))
        indep_pass_all = np.ones(n_pairs, dtype=bool)

    for ch in range(ny):
        e_ch = e[:N_use, ch : ch + 1]  # (N_use, 1)

        # -- Whiteness: normalised autocorrelation --
        R_ee = sid_cov(e_ch, e_ch, max_lag)  # (max_lag+1,) scalar signal
        R_ee0 = R_ee[0]
        if R_ee0 > 0:
            auto_corr_all[:, ch] = R_ee / R_ee0
        whiteness_pass_all[ch] = bool(np.all(np.abs(auto_corr_all[1:, ch]) < conf_bound))

        # -- Independence: normalised cross-correlation per input --
        if not is_time_series:
            for iu in range(nu):
                pair_idx = ch * nu + iu
                u_ch = u_2d[:N_use, iu : iu + 1]  # (N_use, 1)

                R_eu_pos = sid_cov(e_ch, u_ch, max_lag)  # (max_lag+1,)
                R_ue_pos = sid_cov(u_ch, e_ch, max_lag)  # (max_lag+1,)

                R_uu0 = float(u_ch.ravel() @ u_ch.ravel()) / N_use
                denom = np.sqrt(R_ee0 * R_uu0)

                if denom > 0:
                    cc_pos = R_eu_pos / denom
                    cc_neg = R_ue_pos / denom
                else:
                    cc_pos = np.zeros(max_lag + 1)
                    cc_neg = np.zeros(max_lag + 1)

                # Assemble: [cc_neg(M)..cc_neg(1), cc_pos(0)..cc_pos(M)]
                cross_corr_all[:, pair_idx] = np.concatenate([cc_neg[max_lag:0:-1], cc_pos])
                indep_pass_all[pair_idx] = bool(
                    np.all(np.abs(cross_corr_all[:, pair_idx]) < conf_bound)
                )

    # Aggregate
    whiteness_pass: bool = bool(np.all(whiteness_pass_all))
    if not is_time_series:
        independence_pass: bool = bool(np.all(indep_pass_all))
        cross_corr = cross_corr_all
    else:
        independence_pass = True
        cross_corr = np.array([])

    # First-channel summary for convenience
    auto_corr = auto_corr_all[:, 0]

    # ------------------------------------------------------------------
    # Plot (optional)
    # ------------------------------------------------------------------
    if plot:
        _plot_residual_diagnostics(auto_corr, cross_corr, conf_bound, max_lag, is_time_series)

    # ------------------------------------------------------------------
    # Pack result
    # ------------------------------------------------------------------
    indep_pass_all_out: np.ndarray | None = indep_pass_all if not is_time_series else None

    return ResidualResult(
        residual=e,
        auto_corr=auto_corr,
        auto_corr_all=auto_corr_all,
        cross_corr=cross_corr,
        confidence_bound=conf_bound,
        whiteness_pass=whiteness_pass,
        whiteness_pass_all=whiteness_pass_all,
        independence_pass=independence_pass,
        independence_pass_all=indep_pass_all_out,
        data_length=N_eff,
    )


# ======================================================================
# Private helpers
# ======================================================================


def _compute_residual_ss(
    model: object,
    X: np.ndarray,
    U: np.ndarray | None,
) -> tuple[np.ndarray, int]:
    """Residuals from a state-space model (COSMIC).

    Computes ``e[k] = x[k+1] - A[k] @ x[k] - B[k] @ u[k]``, averaged
    over trajectories.

    Parameters
    ----------
    model : object
        Must expose ``a``, ``b``, ``data_length``, ``state_dim``.
    X : ndarray
        State trajectories, ``(N+1, p)`` or ``(N+1, p, L)``.
    U : ndarray or None
        Input data, ``(N, q)`` or ``(N, q, L)``.

    Returns
    -------
    e : ndarray, shape ``(N, p)``
        Residual time series.
    N : int
        Number of time steps.
    """
    X = np.asarray(X, dtype=np.float64)
    Nm: int = model.data_length
    p: int = model.state_dim

    if X.ndim == 3:
        L: int = X.shape[2]
    else:
        L = 1

    e_all = np.zeros((Nm, p), dtype=np.float64)

    for traj in range(L):
        if L > 1:
            Xl = X[:, :, traj]
            Ul = U[:, :, traj] if U is not None else None
        else:
            Xl = X
            Ul = U

        for k in range(Nm):
            x_pred = model.a[:, :, k] @ Xl[k, :].T + model.b[:, :, k] @ Ul[k, :].T
            e_all[k, :] += Xl[k + 1, :] - x_pred.T

    e = e_all / L
    return e, Nm


def _compute_residual_freq(
    model: object,
    y: np.ndarray,
    u: np.ndarray | None,
) -> tuple[np.ndarray, int]:
    """Residuals from a frequency-domain model via IFFT.

    Parameters
    ----------
    model : object
        Must expose ``response`` and ``frequency``.
    y : ndarray, shape ``(N, ny)``
        Measured output.
    u : ndarray or None
        Input signal, ``(N, nu)`` or ``None``.

    Returns
    -------
    e : ndarray, shape ``(N, ny)``
        Residual time series.
    N : int
        Number of samples.
    """
    y = np.asarray(y, dtype=np.float64)
    if y.ndim == 1:
        y = y[:, np.newaxis]

    N: int = y.shape[0]
    ny: int = y.shape[1]

    if u is None:
        # Time-series: residual = y itself
        return y.copy(), N

    u_arr = np.asarray(u, dtype=np.float64)
    if u_arr.ndim == 1:
        u_arr = u_arr[:, np.newaxis]

    nu: int = u_arr.shape[1]
    G_model = model.response

    # Ensure G_model is 3-D
    G_model = np.asarray(G_model)
    if G_model.ndim == 1:
        G_model = G_model[:, np.newaxis, np.newaxis]
    elif G_model.ndim == 2:
        G_model = G_model.reshape(G_model.shape[0], ny, nu)

    y_pred = freq_domain_sim(G_model, model.frequency, u_arr, N)
    e = y - y_pred
    return e, N


def _plot_residual_diagnostics(
    auto_corr: np.ndarray,
    cross_corr: np.ndarray,
    conf_bound: float,
    max_lag: int,
    is_time_series: bool,
) -> None:
    """Two-panel diagnostic plot for residual analysis.

    Parameters
    ----------
    auto_corr : ndarray, shape ``(max_lag+1,)``
        Normalised autocorrelation.
    cross_corr : ndarray
        Normalised cross-correlation (first pair), or empty.
    conf_bound : float
        99% confidence bound.
    max_lag : int
        Maximum lag.
    is_time_series : bool
        Whether the model is time-series (no input).
    """
    import matplotlib.pyplot as plt

    n_panels = 1 if is_time_series else 2
    _, axes = plt.subplots(n_panels, 1, squeeze=False)

    # Top panel: autocorrelation
    ax = axes[0, 0]
    lags_auto = np.arange(max_lag + 1)
    ax.bar(lags_auto, auto_corr, width=0.5, color=(0.3, 0.5, 0.8))
    ax.axhline(conf_bound, color="r", linestyle="--", linewidth=1)
    ax.axhline(-conf_bound, color="r", linestyle="--", linewidth=1)
    violations = np.abs(auto_corr[1:]) >= conf_bound
    if np.any(violations):
        v_idx = np.where(violations)[0] + 1
        ax.bar(
            v_idx,
            auto_corr[v_idx],
            width=0.5,
            color=(0.9, 0.2, 0.2),
        )
    ax.set_xlabel("Lag")
    ax.set_ylabel(r"$r_{ee}(\tau)$")
    ax.set_title("Residual Autocorrelation (Whiteness Test)")

    # Bottom panel: cross-correlation
    if not is_time_series and cross_corr.size > 0:
        ax2 = axes[1, 0]
        lags_cross = np.arange(-max_lag, max_lag + 1)
        cc_first = cross_corr[:, 0] if cross_corr.ndim == 2 else cross_corr
        ax2.bar(lags_cross, cc_first, width=0.5, color=(0.3, 0.5, 0.8))
        ax2.axhline(conf_bound, color="r", linestyle="--", linewidth=1)
        ax2.axhline(-conf_bound, color="r", linestyle="--", linewidth=1)
        violations_c = np.abs(cc_first) >= conf_bound
        if np.any(violations_c):
            v_idx_c = np.where(violations_c)[0]
            ax2.bar(
                lags_cross[v_idx_c],
                cc_first[v_idx_c],
                width=0.5,
                color=(0.9, 0.2, 0.2),
            )
        ax2.set_xlabel("Lag")
        ax2.set_ylabel(r"$r_{eu}(\tau)$")
        ax2.set_title("Residual-Input Cross-Correlation (Independence Test)")

    plt.tight_layout()
    plt.show()
