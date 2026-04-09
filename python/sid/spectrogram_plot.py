# Copyright (c) 2026 Pedro Lourenco. All rights reserved.
# This code is released under the MIT License. See LICENSE file in the
# project root for full license information.
#
# This module is part of the Open Source System Identification Toolbox (SID).
# https://github.com/pdlourenco/sid-matlab

"""Spectrogram colour map plot."""

from __future__ import annotations

import numpy as np

from sid._exceptions import SidError
from sid._results import SpectrogramResult


def spectrogram_plot(
    result: SpectrogramResult,
    *,
    frequency_scale: str = "linear",
    channel: int = 0,
    clim: tuple | None = None,
    ax=None,
) -> dict:
    """Spectrogram colour map plot.

    This is the Python port of ``sidSpectrogramPlot.m``.

    Plots the spectrogram as a time-frequency colour map with power
    in dB.

    Parameters
    ----------
    result : SpectrogramResult
        Result struct returned by :func:`sid.spectrogram`.
    frequency_scale : str, optional
        ``'linear'`` (default) or ``'log'``.
    channel : int, optional
        Zero-based channel index for multi-channel data.  Default is
        ``0``.
    clim : tuple of (float, float) or None, optional
        Colour axis limits ``(cmin, cmax)`` in dB.  ``None`` for
        automatic scaling (default).
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
        If *result* is not a ``SpectrogramResult`` (code:
        ``'invalid_result'``).
    SidError
        If *channel* is out of range (code: ``'invalid_channel'``).

    Examples
    --------
    >>> import numpy as np
    >>> import sid  # doctest: +SKIP
    >>> Fs = 1000; Ts = 1 / Fs; N = 5000
    >>> t = np.arange(N) * Ts
    >>> x = np.cos(2 * np.pi * (50 + 100 * t / t[-1]) * t)
    >>> result = sid.spectrogram(x, window_length=256, sample_time=Ts)  # doctest: +SKIP
    >>> h = sid.spectrogram_plot(result)  # doctest: +SKIP

    Notes
    -----
    **Specification:** SPEC.md S7.5 -- Visualization

    The *channel* parameter is 0-indexed (Python convention), unlike
    MATLAB's 1-indexed ``Channel`` option.

    See Also
    --------
    sid.spectrogram : Compute a spectrogram.
    sid.map_plot : Time-frequency map for freq_map results.

    Changelog
    ---------
    2026-04-09 : First version (Python port) by Pedro Lourenco.
    """
    import matplotlib.pyplot as plt

    # ---- Validate result type ----
    if not hasattr(result, "method") or result.method != "spectrogram":
        raise SidError(
            "invalid_result",
            "Input must be a SpectrogramResult from spectrogram.",
        )

    # ---- Extract data ----
    pdb = np.asarray(result.power_db)

    if pdb.ndim == 3:
        n_ch = pdb.shape[2]
        if channel < 0 or channel >= n_ch:
            raise SidError(
                "invalid_channel",
                f"Channel {channel} out of range (data has {n_ch} channels, 0-indexed).",
            )
        Z = pdb[:, :, channel]
    else:
        if channel != 0:
            raise SidError(
                "invalid_channel",
                f"Channel {channel} out of range (data has 1 channel, 0-indexed).",
            )
        Z = pdb

    T = np.asarray(result.time).ravel()
    F = np.asarray(result.frequency).ravel()

    # ---- Create or reuse axes ----
    if ax is None:
        fig, ax = plt.subplots()
    else:
        fig = ax.get_figure()

    # ---- Plot ----
    mesh = ax.pcolormesh(T, F, Z, shading="auto")
    cb = fig.colorbar(mesh, ax=ax)
    cb.set_label("Power (dB)")

    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Frequency (Hz)")
    ax.set_title(f"Spectrogram (L={result.window_length}, P={result.overlap}, NFFT={result.nfft})")

    if frequency_scale.lower() == "log":
        ax.set_yscale("log")

    if clim is not None:
        mesh.set_clim(clim)

    ax.set_axisbelow(False)

    # ---- Return handles ----
    return {
        "fig": fig,
        "ax": ax,
        "mesh": mesh,
    }
