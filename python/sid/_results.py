# Copyright (c) 2026 Pedro Lourenco. All rights reserved.
# This code is released under the MIT License. See LICENSE file in the
# project root for full license information.
#
# This module is part of the Open Source System Identification Toolbox (SID).
# https://github.com/pdlourenco/sid-matlab

"""Frozen dataclasses for sid result types."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class FreqResult:
    """Result from frequency-domain estimation (freq_bt, freq_etfe, freq_btfdr).

    All array shapes shown for the SISO case.  For MIMO, ``response`` has
    shape ``(nf, ny, nu)`` and ``noise_spectrum`` has shape ``(nf, ny, ny)``.
    """

    frequency: np.ndarray
    """Frequency vector in rad/sample, shape ``(nf,)``."""

    frequency_hz: np.ndarray
    """Frequency vector in Hz, shape ``(nf,)``."""

    response: np.ndarray | None
    """Complex frequency response, shape ``(nf,)`` or ``(nf, ny, nu)``.
    ``None`` in time-series mode."""

    response_std: np.ndarray | None
    """Standard deviation of ``response``, same shape.  ``None`` in time-series mode."""

    noise_spectrum: np.ndarray
    """Noise (or output) power spectrum, shape ``(nf,)`` or ``(nf, ny, ny)``."""

    noise_spectrum_std: np.ndarray
    """Standard deviation of ``noise_spectrum``, same shape."""

    coherence: np.ndarray | None
    """Squared coherence, shape ``(nf,)``.  SISO only; ``None`` for MIMO or time-series."""

    sample_time: float
    """Sample time in seconds."""

    window_size: int | np.ndarray
    """Lag window size M (scalar for BT/ETFE, array for BTFDR)."""

    data_length: int
    """Number of samples N per trajectory."""

    num_trajectories: int
    """Number of trajectories used."""

    method: str
    """Estimation method identifier (``'freq_bt'``, ``'freq_etfe'``, ``'freq_btfdr'``)."""
