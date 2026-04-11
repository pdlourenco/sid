# Copyright (c) 2026 Pedro Lourenco. All rights reserved.
# This code is released under the MIT License. See LICENSE file in the
# project root for full license information.
#
# This module is part of the Open Source System Identification Toolbox (SID).
# https://github.com/pdlourenco/sid

"""Hann (Hanning) lag window for spectral analysis."""

from __future__ import annotations

import numpy as np


def hann_win(M: int) -> np.ndarray:
    """Compute the Hann lag window for lags 0, 1, ..., *M*.

    This is the Python port of ``sidHannWin.m``.

    .. math::

        W(\\tau) = 0.5\\,(1 + \\cos(\\pi\\,\\tau / M))

    so that ``W[0] == 1`` and ``W[M] == 0``.

    Parameters
    ----------
    M : int
        Window size (positive integer, ``M >= 2``).

    Returns
    -------
    W : ndarray, shape ``(M + 1,)``
        Window values for lags ``0 .. M``.

    Examples
    --------
    >>> W = hann_win(30)  # doctest: +SKIP
    >>> W[0]  # doctest: +SKIP
    1.0

    Notes
    -----
    **Specification:** SPEC.md §2.4 -- Hann Lag Window

    See Also
    --------
    sid._internal.cov.sid_cov : Uses this window for lag windowing.
    sid._internal.windowed_dft.windowed_dft : Applies window in spectral estimation.

    Changelog
    ---------
    2026-04-08 : First version by Pedro Lourenco.
    """

    tau = np.arange(M + 1, dtype=np.float64)
    return 0.5 * (1.0 + np.cos(np.pi * tau / M))
