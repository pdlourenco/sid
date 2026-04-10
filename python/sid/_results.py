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


@dataclass(frozen=True)
class LTVResult:
    """Result from LTV state-space identification (ltv_disc).

    Contains the identified time-varying system matrices A(k), B(k) and
    optional Bayesian uncertainty estimates.  All array shapes use the
    convention ``(rows, cols, time)`` consistent with MATLAB's
    ``(:, :, k)`` indexing.
    """

    a: np.ndarray
    """Time-varying dynamics matrices, shape ``(p, p, N)``."""

    b: np.ndarray
    """Time-varying input matrices, shape ``(p, q, N)``."""

    a_std: np.ndarray | None
    """Standard deviation of ``a`` entries, shape ``(p, p, N)``.
    ``None`` when uncertainty was not requested."""

    b_std: np.ndarray | None
    """Standard deviation of ``b`` entries, shape ``(p, q, N)``.
    ``None`` when uncertainty was not requested."""

    p_cov: np.ndarray | None
    """Row-wise posterior covariance blocks, shape ``(d, d, N)`` where
    ``d = p + q``.  ``None`` when uncertainty was not requested."""

    noise_cov: np.ndarray | None
    """Noise covariance matrix, shape ``(p, p)``.
    ``None`` when uncertainty was not requested."""

    noise_cov_estimated: bool | None
    """``True`` if ``noise_cov`` was estimated from residuals,
    ``False`` if user-provided.  ``None`` when uncertainty was not
    requested."""

    noise_variance: float | None
    """Scalar noise variance ``trace(noise_cov) / p``.
    ``None`` when uncertainty was not requested."""

    degrees_of_freedom: float | None
    """Effective degrees of freedom used in noise covariance estimation.
    ``NaN`` when ``noise_cov`` was user-provided.  ``None`` when
    uncertainty was not requested."""

    lambda_: np.ndarray
    """Regularization values used, shape ``(N-1,)``."""

    cost: np.ndarray
    """Cost vector ``[total, fidelity, regularization]``, shape ``(3,)``."""

    data_length: int
    """Number of time steps *N*."""

    state_dim: int
    """State dimension *p*."""

    input_dim: int
    """Input dimension *q*."""

    num_trajectories: int
    """Number of trajectories *L*."""

    algorithm: str
    """Identification algorithm (``'cosmic'``)."""

    preconditioned: bool | str
    """Whether block-diagonal preconditioning was applied.

    ``False`` when not requested, ``'not_implemented'`` when requested
    but disabled (v1.0), or ``True`` when implemented and applied.
    """

    method: str
    """Always ``'ltv_disc'``."""


@dataclass(frozen=True)
class LTVIOResult:
    """Result from LTV input-output identification (ltv_disc_io).

    Contains the identified time-varying system matrices A(k), B(k),
    estimated state trajectories, and optional Bayesian uncertainty
    estimates.  All array shapes use the convention ``(rows, cols, time)``
    consistent with MATLAB's ``(:, :, k)`` indexing.
    """

    a: np.ndarray
    """Time-varying dynamics matrices, shape ``(n, n, N)``."""

    b: np.ndarray
    """Time-varying input matrices, shape ``(n, q, N)``."""

    x: np.ndarray | list
    """Estimated state trajectories, shape ``(N+1, n, L)`` or list of
    ``(N_l+1, n)`` arrays for variable-length trajectories."""

    h: np.ndarray
    """Observation matrix, shape ``(py, n)`` (copy)."""

    r: np.ndarray
    """Measurement noise covariance, shape ``(py, py)`` (copy)."""

    cost: np.ndarray
    """Cost history at each iteration, shape ``(n_iter,)``."""

    iterations: int
    """Number of alternating iterations."""

    lambda_: np.ndarray
    """Regularization values used, shape ``(N-1,)``."""

    data_length: int
    """Number of time steps *N*."""

    state_dim: int
    """State dimension *n*."""

    output_dim: int
    """Output dimension *py*."""

    input_dim: int
    """Input dimension *q*."""

    num_trajectories: int
    """Number of trajectories *L*."""

    a_std: np.ndarray | None
    """Standard deviation of ``a`` entries, shape ``(n, n, N)``.
    ``None`` when uncertainty was not computed."""

    b_std: np.ndarray | None
    """Standard deviation of ``b`` entries, shape ``(n, q, N)``.
    ``None`` when uncertainty was not computed."""

    p_cov: np.ndarray | None
    """Row-wise posterior covariance blocks, shape ``(d, d, N)`` where
    ``d = n + q``.  ``None`` when uncertainty was not computed."""

    noise_cov: np.ndarray | None
    """Noise covariance matrix, shape ``(n, n)``.
    ``None`` when uncertainty was not computed."""

    noise_cov_estimated: bool | None
    """``True`` if ``noise_cov`` was estimated from residuals.
    ``None`` when uncertainty was not computed."""

    noise_variance: float | None
    """Scalar noise variance ``trace(noise_cov) / n``.
    ``None`` when uncertainty was not computed."""

    degrees_of_freedom: float | None
    """Effective degrees of freedom used in noise covariance estimation.
    ``None`` when uncertainty was not computed."""

    algorithm: str
    """Identification algorithm (``'cosmic'``)."""

    method: str
    """Always ``'ltv_disc_io'``."""


@dataclass(frozen=True)
class FrozenResult:
    """Result from frozen transfer function computation (ltv_disc_frozen).

    Contains the instantaneous (frozen) frequency response G(w, k) computed
    from time-varying state-space matrices A(k), B(k), along with optional
    uncertainty propagation.
    """

    frequency: np.ndarray
    """Frequency vector in rad/sample, shape ``(nf,)``."""

    frequency_hz: np.ndarray
    """Frequency vector in Hz, shape ``(nf,)``."""

    time_steps: np.ndarray
    """Selected time step indices (0-based), shape ``(nk,)``."""

    response: np.ndarray
    """Complex frozen transfer function, shape ``(nf, p, q, nk)``."""

    response_std: np.ndarray | None
    """Standard deviation of ``response``, shape ``(nf, p, q, nk)``.
    ``None`` when the input ``LTVResult`` has no uncertainty."""

    sample_time: float
    """Sample time in seconds."""

    method: str
    """Always ``'ltv_disc_frozen'``."""


@dataclass(frozen=True)
class ResidualResult:
    """Result from residual analysis (residual).

    Contains the model residuals, normalised auto- and cross-correlation
    functions, and whiteness/independence diagnostic test outcomes.
    """

    residual: np.ndarray
    """Residual time series, shape ``(N, ny)``."""

    auto_corr: np.ndarray
    """Normalised autocorrelation of the first output channel,
    shape ``(max_lag + 1,)``."""

    auto_corr_all: np.ndarray
    """Per-channel normalised autocorrelation, shape ``(max_lag + 1, ny)``."""

    cross_corr: np.ndarray
    """Normalised cross-correlation between residuals and inputs,
    shape ``(2 * max_lag + 1, ny * nu)``.  Empty array for time-series."""

    confidence_bound: float
    """99% confidence bound ``2.58 / sqrt(N)``."""

    whiteness_pass: bool
    """``True`` if all channels pass the whiteness test."""

    whiteness_pass_all: np.ndarray
    """Per-channel whiteness test result, shape ``(ny,)``."""

    independence_pass: bool
    """``True`` if all (output, input) pairs pass the independence test."""

    independence_pass_all: np.ndarray | None
    """Per-pair independence test result, shape ``(ny * nu,)``.
    ``None`` for time-series mode."""

    data_length: int
    """Number of samples used."""


@dataclass(frozen=True)
class CompareResult:
    """Result from model output comparison (compare).

    Contains the predicted and measured outputs, the NRMSE fit metric
    per channel, and the residual.
    """

    predicted: np.ndarray
    """Model-predicted output, shape ``(N, ny)``."""

    measured: np.ndarray
    """Measured output (copy), shape ``(N, ny)``."""

    fit: np.ndarray
    """NRMSE fit percentage per channel, shape ``(ny,)``.
    100% = perfect, 0% = no better than mean predictor."""

    residual: np.ndarray
    """Residual ``measured - predicted``, shape ``(N, ny)``."""

    method: str
    """Method identifier of the source model."""
