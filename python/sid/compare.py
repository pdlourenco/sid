# Copyright (c) 2026 Pedro Lourenco. All rights reserved.
# This code is released under the MIT License. See LICENSE file in the
# project root for full license information.
#
# This module is part of the Open Source System Identification Toolbox (SID).
# https://github.com/pdlourenco/sid-matlab

"""Compare model predicted output to measured data."""

from __future__ import annotations

import numpy as np

from sid._exceptions import SidError
from sid._internal.freq_domain_sim import freq_domain_sim


def compare(
    model: object,
    y: np.ndarray,
    u: np.ndarray | None = None,
    *,
    initial_state: np.ndarray | None = None,
    plot: bool = False,
) -> dict:
    """Compare model predicted output to measured data.

    This is the Python port of ``sidCompare.m``.

    Simulates the model's predicted output given the input signal and
    compares it to the measured output using the NRMSE fit metric:

    .. math::

        \\text{fit} = 100 \\left(1 -
            \\frac{\\|y - \\hat y\\|}{\\|y - \\bar y\\|}\\right)

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
    initial_state : ndarray or None, optional
        Initial state vector for state-space simulation, shape ``(p,)``.
        Default: first row of *y*.
    plot : bool, optional
        If ``True``, display a comparison plot (requires matplotlib).
        Default: ``False``.

    Returns
    -------
    dict
        Dictionary with the following keys:

        - ``'predicted'`` -- ``(N, ny)`` model-predicted output.
        - ``'measured'`` -- ``(N, ny)`` measured output (copy).
        - ``'fit'`` -- ``(ny,)`` NRMSE fit percentage per channel.
          100% is perfect, 0% is no better than the mean, negative
          values indicate worse than the mean.
        - ``'residual'`` -- ``(N, ny)`` residual ``y - y_pred``.
        - ``'method'`` -- str, method of the source model.

    Raises
    ------
    SidError
        If the model type cannot be determined (code: ``'bad_model'``).

    Examples
    --------
    >>> import sid  # doctest: +SKIP
    >>> G = sid.freq_bt(y, u)  # doctest: +SKIP
    >>> result = sid.compare(G, y, u)  # doctest: +SKIP

    Notes
    -----
    **Specification:** (Model output comparison -- not yet in SPEC.md)

    For state-space models, the simulation propagates:

    .. math::

        x(k+1) = A(k)\\,x(k) + B(k)\\,u(k)

    from the initial state *x0* (default: first row of *y*).  For
    multi-trajectory data the fit is averaged across trajectories.

    See Also
    --------
    sid.residual : Residual analysis and diagnostic tests.
    sid._internal.freq_domain_sim.freq_domain_sim : Frequency-domain
        simulation helper.

    Changelog
    ---------
    2026-04-09 : First version (Python port) by Pedro Lourenco.
    """
    # ------------------------------------------------------------------
    # Dispatch on model type
    # ------------------------------------------------------------------
    if hasattr(model, "a") and hasattr(model, "b"):
        y_pred, y_meas = _simulate_ss(model, y, u, initial_state)
    elif hasattr(model, "response"):
        y_pred, y_meas = _simulate_freq(model, y, u)
    else:
        raise SidError(
            "bad_model",
            "Model must have a 'response' attribute (freq-domain) "
            "or 'a'/'b' attributes (state-space).",
        )

    # ------------------------------------------------------------------
    # Compute NRMSE fit per channel
    # ------------------------------------------------------------------
    ny: int = y_meas.shape[1]
    fit_vec = np.zeros(ny)

    for ch in range(ny):
        ym = y_meas[:, ch]
        yp = y_pred[:, ch]
        denom = np.linalg.norm(ym - np.mean(ym))
        if denom > 0:
            fit_vec[ch] = 100.0 * (1.0 - np.linalg.norm(ym - yp) / denom)
        else:
            fit_vec[ch] = np.nan

    # ------------------------------------------------------------------
    # Determine method name
    # ------------------------------------------------------------------
    if hasattr(model, "method"):
        method_name: str = model.method
    else:
        method_name = "unknown"

    # ------------------------------------------------------------------
    # Plot (optional)
    # ------------------------------------------------------------------
    if plot:
        _plot_comparison(y_meas, y_pred, fit_vec, method_name)

    # ------------------------------------------------------------------
    # Pack result
    # ------------------------------------------------------------------
    return {
        "predicted": y_pred,
        "measured": y_meas,
        "fit": fit_vec,
        "residual": y_meas - y_pred,
        "method": method_name,
    }


# ======================================================================
# Private helpers
# ======================================================================


def _simulate_ss(
    model: object,
    X: np.ndarray,
    U: np.ndarray | None,
    x0: np.ndarray | None,
) -> tuple[np.ndarray, np.ndarray]:
    """Simulate state-space model and return predicted vs measured.

    Parameters
    ----------
    model : object
        Must expose ``a``, ``b``, ``data_length``, ``state_dim``.
    X : ndarray
        State trajectories, ``(N+1, p)`` or ``(N+1, p, L)``.
    U : ndarray or None
        Input data, ``(N, q)`` or ``(N, q, L)``.
    x0 : ndarray or None
        Initial state, ``(p,)`` or ``None``.

    Returns
    -------
    y_pred : ndarray, shape ``(N, p)``
        Predicted state (averaged over trajectories).
    y_meas : ndarray, shape ``(N, p)``
        Measured state (averaged over trajectories).
    """
    X = np.asarray(X, dtype=np.float64)
    Nm: int = model.data_length
    p: int = model.state_dim

    if X.ndim == 3:
        L: int = X.shape[2]
    else:
        L = 1

    y_pred_sum = np.zeros((Nm, p), dtype=np.float64)
    y_meas_sum = np.zeros((Nm, p), dtype=np.float64)

    for traj in range(L):
        if L > 1:
            Xl = X[:, :, traj]
            Ul = U[:, :, traj] if U is not None else None
        else:
            Xl = X
            Ul = U

        # Initial state
        if x0 is not None:
            xk = np.asarray(x0, dtype=np.float64).ravel().copy()
        else:
            xk = Xl[0, :].copy()

        x_hat = np.zeros((Nm, p), dtype=np.float64)
        for k in range(Nm):
            x_next = model.a[:, :, k] @ xk + model.b[:, :, k] @ Ul[k, :]
            x_hat[k, :] = x_next
            xk = x_next

        # Measured: x[1:N+1], Predicted: x_hat[0:N]
        y_meas_sum += Xl[1:, :]
        y_pred_sum += x_hat

    y_meas = y_meas_sum / L
    y_pred = y_pred_sum / L
    return y_pred, y_meas


def _simulate_freq(
    model: object,
    y: np.ndarray,
    u: np.ndarray | None,
) -> tuple[np.ndarray, np.ndarray]:
    """Simulate frequency-domain model output via IFFT.

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
    y_pred : ndarray, shape ``(N, ny)``
        Predicted output.
    y_meas : ndarray, shape ``(N, ny)``
        Measured output (reshaped copy).
    """
    y = np.asarray(y, dtype=np.float64)
    if y.ndim == 1:
        y = y[:, np.newaxis]

    N: int = y.shape[0]
    ny: int = y.shape[1]

    if u is None:
        # Time-series: no input to filter, predicted = 0
        return np.zeros((N, ny), dtype=np.float64), y.copy()

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
    return y_pred, y.copy()


def _plot_comparison(
    y_meas: np.ndarray,
    y_pred: np.ndarray,
    fit_vec: np.ndarray,
    method_name: str,
) -> None:
    """Overlay measured and predicted outputs.

    Parameters
    ----------
    y_meas : ndarray, shape ``(N, ny)``
        Measured output.
    y_pred : ndarray, shape ``(N, ny)``
        Predicted output.
    fit_vec : ndarray, shape ``(ny,)``
        NRMSE fit per channel.
    method_name : str
        Model method identifier.
    """
    import matplotlib.pyplot as plt

    ny: int = y_meas.shape[1]
    N: int = y_meas.shape[0]
    t = np.arange(1, N + 1)

    _, axes = plt.subplots(ny, 1, squeeze=False)
    for ch in range(ny):
        ax = axes[ch, 0]
        ax.plot(t, y_meas[:, ch], "b-", linewidth=1, label="Measured")
        ax.plot(t, y_pred[:, ch], "r--", linewidth=1, label="Predicted")
        if ny > 1:
            ax.set_title(f"Channel {ch + 1} - Fit: {fit_vec[ch]:.1f}%")
        else:
            ax.set_title(f"Model: {method_name} - Fit: {fit_vec[ch]:.1f}%")
        ax.set_xlabel("Sample")
        ax.set_ylabel("Output")
        ax.legend(loc="best")

    plt.tight_layout()
    plt.show()
