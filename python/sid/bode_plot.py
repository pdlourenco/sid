# Copyright (c) 2026 Pedro Lourenco. All rights reserved.
# This code is released under the MIT License. See LICENSE file in the
# project root for full license information.
#
# This module is part of the Open Source System Identification Toolbox (SID).
# https://github.com/pdlourenco/sid-matlab

"""Bode diagram with confidence bands."""

from __future__ import annotations

import numpy as np

from sid._exceptions import SidError
from sid._results import FreqResult

# Machine epsilon for clamping
_EPS = float(np.finfo(np.float64).eps)

# Default MATLAB blue
_DEFAULT_COLOR = "#0072BD"


def bode_plot(
    result: FreqResult,
    *,
    confidence: float = 3.0,
    frequency_unit: str = "rad/s",
    color: str | tuple = None,
    line_width: float = 1.5,
    ax: tuple | None = None,
) -> dict:
    """Bode diagram (magnitude and phase) with shaded confidence bands.

    This is the Python port of ``sidBodePlot.m``.

    Plots the magnitude (dB) and phase (degrees) of the estimated
    frequency response.  When uncertainty information is available, a
    shaded region shows the +/- *confidence*-sigma band.

    Parameters
    ----------
    result : FreqResult
        Result struct returned by :func:`sid.freq_bt`,
        :func:`sid.freq_etfe`, or :func:`sid.freq_btfdr`.  Must contain
        a non-``None`` ``response`` field (input-output mode).
    confidence : float, optional
        Number of standard deviations for the shaded confidence band.
        Set to ``0`` to hide the bands.  Default is ``3.0``.
    frequency_unit : str, optional
        ``'rad/s'`` (default) or ``'Hz'``.
    color : str or tuple, optional
        Line and fill colour.  Default is ``'#0072BD'`` (MATLAB blue).
    line_width : float, optional
        Line width.  Default is ``1.5``.
    ax : tuple of (ax_mag, ax_phase) or None, optional
        Existing matplotlib axes to plot into.  If ``None``, a new
        figure with two vertically-stacked subplots is created.

    Returns
    -------
    dict
        Dictionary with the following keys:

        - ``'fig'`` -- matplotlib Figure handle.
        - ``'ax_mag'`` -- Axes handle for the magnitude subplot.
        - ``'ax_phase'`` -- Axes handle for the phase subplot.
        - ``'line_mag'`` -- Line2D handle for the magnitude trace.
        - ``'line_phase'`` -- Line2D handle for the phase trace.

    Raises
    ------
    SidError
        If ``result.response`` is ``None`` (time-series mode).
        Use :func:`sid.spectrum_plot` instead (code: ``'no_response'``).

    Examples
    --------
    >>> import numpy as np
    >>> import sid  # doctest: +SKIP
    >>> N = 1000; u = np.random.randn(N)
    >>> y = np.convolve(u, [1, -0.9])[:N] + 0.1 * np.random.randn(N)
    >>> result = sid.freq_bt(y, u)  # doctest: +SKIP
    >>> h = sid.bode_plot(result, confidence=3)  # doctest: +SKIP

    Notes
    -----
    **Specification:** (Bode plotting -- not yet in SPEC.md)

    For MIMO systems only the first SISO channel pair ``(0, 0)`` is
    plotted.

    See Also
    --------
    sid.freq_bt : Blackman-Tukey spectral analysis.
    sid.spectrum_plot : Noise / output spectrum plot.

    Changelog
    ---------
    2026-04-09 : First version (Python port) by Pedro Lourenco.
    """
    import matplotlib.pyplot as plt

    # ---- Validate ----
    if result.response is None:
        raise SidError(
            "no_response",
            "Result contains no frequency response (time-series mode). Use spectrum_plot instead.",
        )

    if color is None:
        color = _DEFAULT_COLOR

    # ---- Frequency axis ----
    if frequency_unit.lower() == "hz":
        freq = result.frequency_hz
        freq_label = "Frequency (Hz)"
    else:
        freq = result.frequency / result.sample_time
        freq_label = "Frequency (rad/s)"

    # ---- Extract first SISO channel (ny=0, nu=0 for MIMO) ----
    G = np.asarray(result.response)
    if G.ndim == 3:
        G = G[:, 0, 0]
    elif G.ndim == 2:
        G = G[:, 0]

    G_std = None
    if result.response_std is not None:
        G_std = np.asarray(result.response_std)
        if G_std.ndim == 3:
            G_std = G_std[:, 0, 0]
        elif G_std.ndim == 2:
            G_std = G_std[:, 0]

    mag = np.abs(G)
    mag_db = 20.0 * np.log10(np.maximum(mag, _EPS))
    phase = np.angle(G) * 180.0 / np.pi

    # ---- Create or reuse axes ----
    if ax is None:
        fig, (ax_mag, ax_phase) = plt.subplots(2, 1)
    else:
        ax_mag, ax_phase = ax
        fig = ax_mag.get_figure()

    # ---- Magnitude plot ----
    (line_mag,) = ax_mag.semilogx(
        freq,
        mag_db,
        color=color,
        linewidth=line_width,
    )

    if confidence > 0 and G_std is not None:
        mag_upper = 20.0 * np.log10(np.maximum(mag + confidence * G_std, _EPS))
        mag_lower = 20.0 * np.log10(np.maximum(mag - confidence * G_std, _EPS))
        ax_mag.fill_between(
            freq,
            mag_lower,
            mag_upper,
            color=color,
            alpha=0.15,
            edgecolor="none",
        )

    ax_mag.set_ylabel("Magnitude (dB)")
    ax_mag.set_title(f"Bode Diagram ({result.method})")
    ax_mag.grid(True)

    # ---- Phase plot ----
    (line_phase,) = ax_phase.semilogx(
        freq,
        phase,
        color=color,
        linewidth=line_width,
    )

    if confidence > 0 and G_std is not None:
        phase_std = confidence * G_std / np.maximum(mag, _EPS) * 180.0 / np.pi
        ax_phase.fill_between(
            freq,
            phase - phase_std,
            phase + phase_std,
            color=color,
            alpha=0.15,
            edgecolor="none",
        )

    ax_phase.set_xlabel(freq_label)
    ax_phase.set_ylabel("Phase (deg)")
    ax_phase.grid(True)

    # ---- Return handles ----
    return {
        "fig": fig,
        "ax_mag": ax_mag,
        "ax_phase": ax_phase,
        "line_mag": line_mag,
        "line_phase": line_phase,
    }
