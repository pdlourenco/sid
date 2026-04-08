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


@dataclass(frozen=True)
class SpectrogramResult:
    """Result from short-time FFT spectrogram (spectrogram)."""

    time: np.ndarray
    """Center time of each segment in seconds, shape ``(K,)``."""

    frequency: np.ndarray
    """Frequency vector in Hz, shape ``(n_bins,)``."""

    frequency_rad: np.ndarray
    """Frequency vector in rad/s, shape ``(n_bins,)``."""

    power: np.ndarray
    """Power spectral density, shape ``(n_bins, K)`` or ``(n_bins, K, n_ch)``."""

    power_db: np.ndarray
    """Power in dB (``10 * log10(power)``), same shape."""

    complex_stft: np.ndarray
    """Complex STFT coefficients, same shape."""

    sample_time: float
    """Sample time in seconds."""

    window_length: int
    """Segment length L."""

    overlap: int
    """Overlap P between segments."""

    nfft: int
    """FFT length."""

    num_trajectories: int
    """Number of trajectories used."""

    method: str
    """Always ``'spectrogram'``."""


@dataclass(frozen=True)
class FreqMapResult:
    """Result from time-varying frequency response map (freq_map)."""

    time: np.ndarray
    """Center time of each segment in seconds, shape ``(K,)``."""

    frequency: np.ndarray
    """Frequency vector in rad/sample, shape ``(nf,)``."""

    frequency_hz: np.ndarray
    """Frequency vector in Hz, shape ``(nf,)``."""

    response: np.ndarray | None
    """Time-varying frequency response, shape ``(nf, K)`` or ``(nf, K, ny, nu)``.
    ``None`` in time-series mode."""

    response_std: np.ndarray | None
    """Standard deviation of response, same shape."""

    noise_spectrum: np.ndarray
    """Time-varying noise spectrum, shape ``(nf, K)`` or ``(nf, K, ny, ny)``."""

    noise_spectrum_std: np.ndarray
    """Standard deviation of noise spectrum."""

    coherence: np.ndarray | None
    """Squared coherence, shape ``(nf, K)``.  SISO only; ``None`` otherwise."""

    sample_time: float
    """Sample time in seconds."""

    segment_length: int
    """Segment length L."""

    overlap: int
    """Overlap P between segments."""

    window_size: int | None
    """BT lag window size M, or ``None`` for Welch."""

    algorithm: str
    """``'bt'`` or ``'welch'``."""

    num_trajectories: int
    """Number of trajectories used."""

    method: str
    """Always ``'freq_map'``."""
