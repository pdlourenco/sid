# Copyright (c) 2026 Pedro Lourenco. All rights reserved.
# This code is released under the MIT License. See LICENSE file in the
# project root for full license information.
#
# This module is part of the Open Source System Identification Toolbox (SID).
# https://github.com/pdlourenco/sid-matlab

"""Check whether a frequency vector matches the default 128-point grid."""

from __future__ import annotations

import numpy as np


def is_default_freqs(freqs: np.ndarray) -> bool:
    """Check if *freqs* matches the default 128-point linear grid.

    This is the Python port of ``sidIsDefaultFreqs.m``.

    The default grid is ``k * pi / 128`` for ``k = 1, 2, ..., 128``.
    This function is used to decide whether the FFT fast path can be
    applied.

    Parameters
    ----------
    freqs : ndarray, shape ``(nf,)``
        Frequency vector in rad/sample.

    Returns
    -------
    bool
        ``True`` if *freqs* matches the default grid to within ``1e-12``.

    Examples
    --------
    >>> is_default_freqs(np.arange(1, 129) * np.pi / 128)  # doctest: +SKIP
    True

    See Also
    --------
    sid.freq_bt : Uses this to choose FFT vs direct DFT path.

    Changelog
    ---------
    2026-04-08 : First version by Pedro Lourenco.
    """

    if len(freqs) != 128:
        return False
    default = np.arange(1, 129) * np.pi / 128
    return float(np.max(np.abs(freqs.ravel() - default))) < 1e-12
