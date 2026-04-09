# Copyright (c) 2026 Pedro Lourenco. All rights reserved.
# This code is released under the MIT License. See LICENSE file in the
# project root for full license information.
#
# This module is part of the Open Source System Identification Toolbox (SID).
# https://github.com/pdlourenco/sid-matlab

"""Time-frequency colour map for sidFreqMap results."""

from __future__ import annotations

import numpy as np

from sid._exceptions import SidError
from sid._results import FreqMapResult

# Machine epsilon for clamping
_EPS = float(np.finfo(np.float64).eps)


def map_plot(
    result: FreqMapResult,
    *,
    plot_type: str = "magnitude",
    frequency_unit: str = "rad/s",
    clim: tuple | None = None,
    ax=None,
) -> dict:
    """Time-frequency colour map for :func:`sid.freq_map` results.

    This is the Python port of ``sidMapPlot.m``.

    Plots the time-varying frequency response (or spectrum / coherence)
    as a ``pcolormesh`` colour map with time on the x-axis and
    frequency (log scale) on the y-axis.

    Parameters
    ----------
    result : FreqMapResult
        Result struct returned by :func:`sid.freq_map`.  Must have
        ``method == 'freq_map'``.
    plot_type : str, optional
        What to plot:

        - ``'magnitude'`` -- ``20 * log10(|G(w,t)|)`` in dB (default).
        - ``'phase'`` -- ``angle(G(w,t))`` in degrees.
        - ``'noise'`` -- ``10 * log10(noise_spectrum)`` in dB.
        - ``'coherence'`` -- squared coherence on ``[0, 1]``.
        - ``'spectrum'`` -- ``10 * log10(noise_spectrum)`` in dB
          (alias for time-series data).
    frequency_unit : str, optional
        ``'rad/s'`` (default) or ``'Hz'``.
    clim : tuple of (float, float) or None, optional
        Colour axis limits ``(cmin, cmax)``.  ``None`` for automatic
        scaling (default).
    ax : matplotlib Axes or None, optional
        Existing axes to plot into.  If ``None``, a new figure is
        created.

    Returns
    -------
    dict
        Dictionary with the following keys:

        - ``'fig'`` -- matplotlib Figure handle.
        - ``'ax'`` -- Axes handle.
        - ``'mesh'`` -- QuadMesh handle from ``pcolormesh``.

    Raises
    ------
    SidError
        If *result* is not a ``FreqMapResult`` (code:
        ``'invalid_result'``).
    SidError
        If *plot_type* requires a frequency response but the result is
        time-series only (code: ``'no_response'``).
    SidError
        If *plot_type* is ``'coherence'`` but no coherence data is
        available (code: ``'no_coherence'``).
    SidError
        If *plot_type* is unrecognised (code: ``'invalid_plot_type'``).

    Examples
    --------
    >>> import numpy as np
    >>> import sid  # doctest: +SKIP
    >>> N = 4000; u = np.random.randn(N)
    >>> y = np.convolve(u, [1, -0.9])[:N] + 0.1 * np.random.randn(N)
    >>> result = sid.freq_map(y, u, segment_length=512)  # doctest: +SKIP
    >>> h = sid.map_plot(result, plot_type='magnitude')  # doctest: +SKIP

    Notes
    -----
    **Specification:** SPEC.md S6.9 -- Visualization: sidMapPlot

    See Also
    --------
    sid.freq_map : Time-varying frequency response estimation.
    sid.spectrogram_plot : Spectrogram colour map.
    sid.bode_plot : Static Bode diagram.

    Changelog
    ---------
    2026-04-09 : First version (Python port) by Pedro Lourenco.
    """
    import matplotlib.pyplot as plt

    # ---- Validate result type ----
    if not hasattr(result, "method") or result.method != "freq_map":
        raise SidError(
            "invalid_result",
            "Input must be a FreqMapResult from freq_map.",
        )

    # ---- Frequency axis ----
    if frequency_unit.lower() == "hz":
        freq = result.frequency_hz
        freq_label = "Frequency (Hz)"
    else:
        freq = result.frequency / result.sample_time
        freq_label = "Frequency (rad/s)"

    T = np.asarray(result.time).ravel()
    F = np.asarray(freq).ravel()

    is_time_series = result.response is None

    # ---- Select data based on plot_type ----
    pt = plot_type.lower()

    if pt == "magnitude":
        if is_time_series:
            raise SidError(
                "no_response",
                "PlotType 'magnitude' requires input-output data (not time series).",
            )
        resp = np.asarray(result.response)
        if resp.ndim > 2:
            resp = resp[:, :, 0]
        Z = 20.0 * np.log10(np.maximum(np.abs(resp), _EPS))
        color_label = "Magnitude (dB)"
        title_str = "Time-Varying Magnitude"

    elif pt == "phase":
        if is_time_series:
            raise SidError(
                "no_response",
                "PlotType 'phase' requires input-output data (not time series).",
            )
        resp = np.asarray(result.response)
        if resp.ndim > 2:
            resp = resp[:, :, 0]
        Z = np.angle(resp) * 180.0 / np.pi
        color_label = "Phase (deg)"
        title_str = "Time-Varying Phase"

    elif pt == "noise":
        ns = np.asarray(result.noise_spectrum)
        if ns.ndim > 2:
            ns = ns[:, :, 0]
        Z = 10.0 * np.log10(np.maximum(ns, _EPS))
        color_label = "Noise PSD (dB)"
        title_str = "Time-Varying Noise Spectrum"

    elif pt == "coherence":
        if result.coherence is None:
            raise SidError(
                "no_coherence",
                "Coherence is only available for SISO input-output data.",
            )
        Z = np.asarray(result.coherence)
        color_label = "Coherence"
        title_str = "Time-Varying Coherence"

    elif pt == "spectrum":
        ns = np.asarray(result.noise_spectrum)
        if ns.ndim > 2:
            ns = ns[:, :, 0]
        Z = 10.0 * np.log10(np.maximum(ns, _EPS))
        color_label = "PSD (dB)"
        title_str = "Time-Varying Power Spectrum"

    else:
        raise SidError(
            "invalid_plot_type",
            f"Unknown plot_type '{plot_type}'. "
            "Use magnitude, phase, noise, coherence, or spectrum.",
        )

    # ---- Create or reuse axes ----
    if ax is None:
        fig, ax = plt.subplots()
    else:
        fig = ax.get_figure()

    # ---- Plot ----
    mesh = ax.pcolormesh(T, F, Z, shading="auto")
    ax.set_yscale("log")
    cb = fig.colorbar(mesh, ax=ax)
    cb.set_label(color_label)

    ax.set_xlabel("Time (s)")
    ax.set_ylabel(freq_label)
    ax.set_title(f"{title_str} (L={result.segment_length}, M={result.window_size})")

    if clim is not None:
        mesh.set_clim(clim)

    ax.set_axisbelow(False)

    # ---- Return handles ----
    return {
        "fig": fig,
        "ax": ax,
        "mesh": mesh,
    }
