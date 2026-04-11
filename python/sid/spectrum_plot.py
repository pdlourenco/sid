# Copyright (c) 2026 Pedro Lourenco. All rights reserved.
# This code is released under the MIT License. See LICENSE file in the
# project root for full license information.
#
# This module is part of the Open Source System Identification Toolbox (SID).
# https://github.com/pdlourenco/sid

"""Power spectrum plot with confidence bands."""

from __future__ import annotations

import numpy as np

from sid._results import FreqResult

# Machine epsilon for clamping
_EPS = float(np.finfo(np.float64).eps)

# Default MATLAB orange
_DEFAULT_COLOR = "#D95319"


def spectrum_plot(
    result: FreqResult,
    *,
    confidence: float = 3.0,
    frequency_unit: str = "rad/s",
    color: str | tuple = None,
    line_width: float = 1.5,
    ax=None,
) -> dict:
    """Power spectrum plot with shaded confidence bands.

    This is the Python port of ``sidSpectrumPlot.m``.

    Plots the noise spectrum (or output spectrum in time-series mode)
    in dB with an optional shaded +/- *confidence*-sigma band.

    Parameters
    ----------
    result : FreqResult
        Result struct returned by :func:`sid.freq_bt`,
        :func:`sid.freq_etfe`, or :func:`sid.freq_btfdr`.
    confidence : float, optional
        Number of standard deviations for the shaded confidence band.
        Set to ``0`` to hide the bands.  Default is ``3.0``.
    frequency_unit : str, optional
        ``'rad/s'`` (default) or ``'Hz'``.
    color : str or tuple, optional
        Line and fill colour.  Default is ``'#D95319'`` (MATLAB orange).
    line_width : float, optional
        Line width.  Default is ``1.5``.
    ax : matplotlib Axes or None, optional
        Existing axes to plot into.  If ``None``, a new figure is
        created.

    Returns
    -------
    dict
        Dictionary with the following keys:

        - ``'fig'`` -- matplotlib Figure handle.
        - ``'ax'`` -- Axes handle.
        - ``'line'`` -- Line2D handle for the spectrum trace.

    Examples
    --------
    >>> import numpy as np
    >>> import sid  # doctest: +SKIP
    >>> y = np.random.randn(500)
    >>> result = sid.freq_bt(y, None)  # doctest: +SKIP
    >>> h = sid.spectrum_plot(result, confidence=3)  # doctest: +SKIP

    Notes
    -----
    **Specification:** (Spectrum plotting -- not yet in SPEC.md)

    For MIMO systems only the first output channel is plotted.

    See Also
    --------
    sid.freq_bt : Blackman-Tukey spectral analysis.
    sid.bode_plot : Bode diagram for frequency-response data.

    Changelog
    ---------
    2026-04-09 : First version (Python port) by Pedro Lourenco.
    """
    import matplotlib.pyplot as plt

    if color is None:
        color = _DEFAULT_COLOR

    # ---- Frequency axis ----
    if frequency_unit.lower() == "hz":
        freq = result.frequency_hz
        freq_label = "Frequency (Hz)"
    else:
        freq = result.frequency / result.sample_time
        freq_label = "Frequency (rad/s)"

    # ---- Extract first output channel ----
    PhiV = np.asarray(result.noise_spectrum)
    if PhiV.ndim > 1:
        PhiV = PhiV[:, 0]

    PhiV_std = None
    if result.noise_spectrum_std is not None:
        PhiV_std = np.asarray(result.noise_spectrum_std)
        if PhiV_std.ndim > 1:
            PhiV_std = PhiV_std[:, 0]

    spec_db = 10.0 * np.log10(np.maximum(PhiV, _EPS))

    # ---- Create or reuse axes ----
    if ax is None:
        fig, ax = plt.subplots()
    else:
        fig = ax.get_figure()

    # ---- Plot ----
    (line,) = ax.semilogx(freq, spec_db, color=color, linewidth=line_width)

    if confidence > 0 and PhiV_std is not None:
        upper = 10.0 * np.log10(np.maximum(PhiV + confidence * PhiV_std, _EPS))
        lower = 10.0 * np.log10(np.maximum(PhiV - confidence * PhiV_std, _EPS))
        ax.fill_between(
            freq,
            lower,
            upper,
            color=color,
            alpha=0.15,
            edgecolor="none",
        )

    # ---- Labels ----
    ax.set_xlabel(freq_label)
    if result.response is None:
        ax.set_ylabel("Output Spectrum (dB)")
        title_str = "Output Power Spectrum"
    else:
        ax.set_ylabel("Noise Spectrum (dB)")
        title_str = "Noise Spectrum"

    ax.set_title(f"{title_str} ({result.method})")
    ax.grid(True)

    # ---- Return handles ----
    return {
        "fig": fig,
        "ax": ax,
        "line": line,
    }
