# Copyright (c) 2026 Pedro Lourenco. All rights reserved.
# This code is released under the MIT License. See LICENSE file in the
# project root for full license information.
#
# This module is part of the Open Source System Identification Toolbox (SID).
# https://github.com/pdlourenco/sid-matlab

"""Time-varying frequency response map via segmented spectral analysis."""

from __future__ import annotations

import math

import numpy as np

from sid._exceptions import SidError
from sid._internal.validate_data import validate_data
from sid._results import FreqMapResult, FreqResult
from sid.freq_bt import freq_bt


# ---- Welch inner estimator (private) ----------------------------------------


def _welch_estimate(
    y: np.ndarray,
    u: np.ndarray | None,
    is_time_series: bool,
    ny: int,
    nu: int,
    Lsub: int,
    Psub: int,
    nfft: int,
    win_type: str | np.ndarray,
    Ts: float,
) -> FreqResult:
    """Welch averaged periodogram estimate within one outer segment.

    Returns a :class:`FreqResult` with the same field semantics as
    :func:`freq_bt` output.  Supports multi-trajectory input where *y*
    has shape ``(Lseg, ny, L)`` and *u* has shape ``(Lseg, nu, L)``.
    """
    Lseg = y.shape[0]
    n_traj_w = y.shape[2] if y.ndim == 3 else 1

    # ---- Build window ----
    if isinstance(win_type, np.ndarray) or (
        not isinstance(win_type, str) and hasattr(win_type, "__len__")
    ):
        w = np.asarray(win_type, dtype=np.float64).ravel()
        if len(w) != Lsub:
            raise SidError(
                "invalid_window",
                f"Window vector length ({len(w)}) must match SubSegmentLength ({Lsub}).",
            )
    else:
        name = win_type.lower()
        if name == "hann":
            n = np.arange(Lsub, dtype=np.float64)
            w = 0.5 * (1.0 - np.cos(2.0 * np.pi * n / (Lsub - 1)))
        elif name == "hamming":
            n = np.arange(Lsub, dtype=np.float64)
            w = 0.54 - 0.46 * np.cos(2.0 * np.pi * n / (Lsub - 1))
        elif name == "rect":
            w = np.ones(Lsub, dtype=np.float64)
        else:
            raise SidError(
                "invalid_window",
                f"Window must be 'hann', 'hamming', 'rect', or a numeric vector. Got '{win_type}'.",
            )

    S1 = float(np.sum(w**2))  # window power normalisation

    # ---- Sub-segmentation ----
    sub_step = Lsub - Psub
    J = (Lseg - Lsub) // sub_step + 1
    if J < 1:
        raise SidError(
            "too_few_sub_segments",
            f"Segment too short for sub-segmentation with Lsub={Lsub}, Psub={Psub}.",
        )

    # ---- One-sided frequency grid (skip DC, up to Nyquist) ----
    n_bins = nfft // 2  # number of one-sided bins (excluding DC)
    freqs = np.arange(1, n_bins + 1) * (2.0 * np.pi / nfft)  # rad/sample

    # ---- Accumulate averaged periodograms ----
    PhiY = np.zeros((n_bins, ny, ny), dtype=complex)
    if not is_time_series:
        PhiU = np.zeros((n_bins, nu, nu), dtype=complex)
        PhiYU = np.zeros((n_bins, ny, nu), dtype=complex)

    for lt in range(n_traj_w):
        for j in range(J):
            s = j * sub_step
            e = s + Lsub

            # Windowed FFT of output
            if n_traj_w > 1:
                y_seg = y[s:e, :, lt]
            else:
                y_seg = y[s:e, :]
            Yj = np.fft.fft(y_seg * w[:, np.newaxis], n=nfft, axis=0)
            Yj = Yj[1 : n_bins + 1, :]  # one-sided, skip DC

            # Auto-spectrum of y
            for a in range(ny):
                for b in range(ny):
                    PhiY[:, a, b] += Yj[:, a] * np.conj(Yj[:, b])

            if not is_time_series:
                # Windowed FFT of input
                if n_traj_w > 1:
                    u_seg = u[s:e, :, lt]
                else:
                    u_seg = u[s:e, :]
                Uj = np.fft.fft(u_seg * w[:, np.newaxis], n=nfft, axis=0)
                Uj = Uj[1 : n_bins + 1, :]

                # Auto-spectrum of u
                for a in range(nu):
                    for b in range(nu):
                        PhiU[:, a, b] += Uj[:, a] * np.conj(Uj[:, b])

                # Cross-spectrum y*u'
                for a in range(ny):
                    for b in range(nu):
                        PhiYU[:, a, b] += Yj[:, a] * np.conj(Uj[:, b])

    # ---- Average and normalise ----
    # Factor of 2 for one-sided spectrum (positive frequencies only,
    # excluding DC).  Cancels in G = Pyu/Puu and coherence but is needed
    # for correct PSD magnitude.
    J_total = J * n_traj_w
    PhiY = 2.0 * PhiY / (J_total * S1)
    if not is_time_series:
        PhiU = 2.0 * PhiU / (J_total * S1)
        PhiYU = 2.0 * PhiYU / (J_total * S1)

    # ---- Degrees of freedom for uncertainty ----
    # For Hann window at 50% overlap: nu_dof ~ 1.8 * J per trajectory.
    # Multi-trajectory multiplies by n_traj_w (independent realisations).
    overlap_ratio = Psub / Lsub
    if overlap_ratio <= 0:
        nu_dof = 2.0 * J * n_traj_w
    else:
        nu_dof = max(2.0, 1.8 * J * n_traj_w)

    # ---- Form transfer function, noise spectrum, coherence ----
    nf = n_bins
    if is_time_series:
        G: np.ndarray | None = None
        GStd: np.ndarray | None = None
        PhiV = np.real(PhiY).squeeze()
        if ny == 1:
            PhiV = PhiV.ravel()
        PhiVStd = PhiV * math.sqrt(2.0 / nu_dof)
        Coh: np.ndarray | None = None

    elif ny == 1 and nu == 1:
        # SISO
        PhiY_s = PhiY.ravel()
        PhiU_s = PhiU.ravel()
        PhiYU_s = PhiYU.ravel()
        G = PhiYU_s / PhiU_s
        PhiV = np.real(PhiY_s) - np.abs(PhiYU_s) ** 2 / np.real(PhiU_s)
        PhiV = np.maximum(PhiV, 0.0)
        Coh = np.abs(PhiYU_s) ** 2 / (np.real(PhiY_s) * np.real(PhiU_s))
        Coh = np.clip(Coh, 0.0, 1.0)
        GStd = np.abs(G) * np.sqrt((1.0 - Coh) / (Coh * nu_dof))
        GStd = np.where(Coh < 1e-10, np.nan, GStd)
        PhiVStd = PhiV * math.sqrt(2.0 / nu_dof)

    else:
        # MIMO
        G = np.zeros((nf, ny, nu), dtype=complex)
        PhiV = np.zeros((nf, ny, ny))
        for k in range(nf):
            PhiU_k = PhiU[k, :, :].reshape(nu, nu)
            PhiYU_k = PhiYU[k, :, :].reshape(ny, nu)
            PhiY_k = PhiY[k, :, :].reshape(ny, ny)
            # G(k) = Phi_yu(k) * Phi_u(k)^{-1}
            G[k, :, :] = np.linalg.solve(PhiU_k.T, PhiYU_k.T).T
            # Phi_v = Phi_y - Phi_yu * Phi_u^{-1} * Phi_yu'
            Gk = G[k, :, :]
            PhiV[k, :, :] = np.real(PhiY_k - Gk @ PhiYU_k.conj().T)
        PhiV = np.real(PhiV)
        GStd = np.full_like(G, np.nan, dtype=float)
        PhiVStd = np.abs(PhiV) * math.sqrt(2.0 / nu_dof)
        Coh = None

    # ---- Pack result (matching freq_bt output structure) ----
    return FreqResult(
        frequency=freqs,
        frequency_hz=freqs / (2.0 * np.pi * Ts),
        response=G,
        response_std=GStd,
        noise_spectrum=PhiV,
        noise_spectrum_std=PhiVStd,
        coherence=Coh,
        sample_time=Ts,
        window_size=0,
        data_length=Lseg,
        num_trajectories=n_traj_w,
        method="welch",
    )


# ---- Public API --------------------------------------------------------------


def freq_map(
    y: np.ndarray | list,
    u: np.ndarray | list | None = None,
    *,
    segment_length: int | None = None,
    overlap: int | None = None,
    algorithm: str = "bt",
    sample_time: float = 1.0,
    # BT-specific
    window_size: int | None = None,
    frequencies: np.ndarray | None = None,
    # Welch-specific
    sub_segment_length: int | None = None,
    sub_overlap: int | None = None,
    window: str | np.ndarray = "hann",
    nfft: int | None = None,
) -> FreqMapResult:
    """Estimate a time-varying frequency response map.

    Estimates G(omega, t) by applying spectral analysis to overlapping
    segments of input-output data.  For a linear time-invariant (LTI)
    system the map is constant along time; for a linear time-varying
    (LTV) system it reveals how the transfer function, noise spectrum,
    and coherence evolve.

    Two inner estimation algorithms are supported:

    * ``'bt'`` (default) -- Blackman-Tukey correlogram via
      :func:`sid.freq_bt`.
    * ``'welch'`` -- Welch averaged periodogram.

    Parameters
    ----------
    y : ndarray or list of ndarray, shape (N, ny) or (N, ny, L)
        Output data.  A 1-D array is treated as a single-channel signal.
        For multiple trajectories pass a 3-D array ``(N, ny, L)`` or a
        list of 1-D / 2-D arrays.  Spectral estimates within each
        segment are ensemble-averaged.
    u : ndarray, list of ndarray, or None, optional
        Input data.  Same shape conventions as *y*.  Pass ``None`` for
        time-series mode (output spectrum only).  Default is ``None``.
    segment_length : int, optional
        Number of samples per outer segment *L*.  Must be >= 4.
        Default: ``min(N // 4, 256)``.
    overlap : int, optional
        Overlap *P* between segments, ``0 <= P < L``.
        Default: ``L // 2``.
    algorithm : str, optional
        ``'bt'`` (default) or ``'welch'``.
    sample_time : float, optional
        Sample time in seconds.  Must be positive.  Default is ``1.0``.
    window_size : int, optional
        *(BT only)* Hann lag window size *M*.
        Default: ``min(L // 10, 30)``.
    frequencies : ndarray, shape (nf,), optional
        *(BT only)* Frequency vector in rad/sample, in ``(0, pi]``.
        Default: 128-point linear grid.
    sub_segment_length : int, optional
        *(Welch only)* Sub-segment length within each outer segment.
        Default: ``floor(L / 4.5)``.
    sub_overlap : int, optional
        *(Welch only)* Sub-segment overlap.
        Default: ``sub_segment_length // 2``.
    window : str or ndarray, optional
        *(Welch only)* ``'hann'`` (default), ``'hamming'``, ``'rect'``,
        or a numeric window vector of length *sub_segment_length*.
    nfft : int, optional
        *(Welch only)* FFT length.
        Default: ``max(256, 2**ceil(log2(sub_segment_length)))``.

    Returns
    -------
    FreqMapResult
        Frozen dataclass with fields:

        - **time** (*ndarray, shape (K,)*) -- Centre time of each
          segment in seconds.
        - **frequency** (*ndarray, shape (nf,)*) -- Frequency vector,
          rad/sample.
        - **frequency_hz** (*ndarray, shape (nf,)*) -- Frequency
          vector, Hz.
        - **response** (*ndarray or None*) -- Complex time-varying
          frequency response.  Shape ``(nf, K)`` for SISO,
          ``(nf, K, ny, nu)`` for MIMO, or ``None`` in time-series
          mode.
        - **response_std** (*ndarray or None*) -- Standard deviation of
          ``response``, same shape.
        - **noise_spectrum** (*ndarray*) -- Time-varying noise spectrum.
          Shape ``(nf, K)`` for SISO / time-series,
          ``(nf, K, ny, ny)`` for MIMO.
        - **noise_spectrum_std** (*ndarray*) -- Standard deviation of
          ``noise_spectrum``, same shape.
        - **coherence** (*ndarray or None*) -- Squared coherence, shape
          ``(nf, K)``.  SISO only; ``None`` for MIMO or time-series.
        - **sample_time** (*float*) -- Sample time in seconds.
        - **segment_length** (*int*) -- Segment length *L*.
        - **overlap** (*int*) -- Overlap *P* between segments.
        - **window_size** (*int or None*) -- BT lag window size *M*, or
          ``None`` for Welch.
        - **algorithm** (*str*) -- ``'bt'`` or ``'welch'``.
        - **num_trajectories** (*int*) -- Number of trajectories.
        - **method** (*str*) -- Always ``'freq_map'``.

    Raises
    ------
    SidError
        If ``segment_length < 4`` (code: ``'invalid_segment_length'``).
    SidError
        If ``segment_length > N`` (code: ``'segment_too_long'``).
    SidError
        If ``sample_time <= 0`` (code: ``'invalid_sample_time'``).
    SidError
        If ``algorithm`` is not ``'bt'`` or ``'welch'``
        (code: ``'invalid_algorithm'``).
    SidError
        If ``overlap`` is not in ``[0, L-1]``
        (code: ``'invalid_overlap'``).
    SidError
        If data is too short for even one segment
        (code: ``'too_few_segments'``).

    Examples
    --------
    Time-varying frequency map (Blackman-Tukey):

    >>> import numpy as np
    >>> import sid  # doctest: +SKIP
    >>> N = 4000; rng = np.random.default_rng(0)
    >>> u = rng.standard_normal(N)
    >>> from scipy.signal import lfilter  # doctest: +SKIP
    >>> y = lfilter([1], [1, -0.9], u) + 0.1 * rng.standard_normal(N)  # doctest: +SKIP
    >>> result = sid.freq_map(y, u, segment_length=512)  # doctest: +SKIP

    Time-varying spectrum (time-series, Welch):

    >>> result = sid.freq_map(y, algorithm='welch', segment_length=256)  # doctest: +SKIP

    Notes
    -----
    **Algorithm (outer segmentation + inner estimator):**

    1. Validate and orient data; detect time-series / SISO / MIMO.
    2. Divide the data into *K* overlapping segments of length *L*
       with stride ``L - P``.
    3. For each segment, apply either the Blackman-Tukey or Welch
       inner estimator.
    4. Stack the per-segment results into 2-D / 4-D arrays indexed by
       ``(frequency, segment)``.

    **Specification:** SPEC.md S6 -- Time-Varying Frequency Response Map

    See Also
    --------
    sid.freq_bt : Blackman-Tukey spectral analysis (inner estimator).
    sid.map_plot : Surface / image plot of a FreqMapResult.
    sid.spectrogram : Short-time FFT spectrogram.

    Changelog
    ---------
    2026-04-08 : First version (Python port) by Pedro Lourenco.
    """
    # ---- Validate data ----
    y, u, N, ny, nu, is_time_series, n_traj = validate_data(y, u)

    # ---- Defaults ----
    L = segment_length if segment_length is not None else min(N // 4, 256)
    Ts = sample_time
    alg = algorithm.lower()

    # ---- Validate common parameters ----
    if not isinstance(L, (int, np.integer)) or L != int(L) or L < 4:
        raise SidError(
            "invalid_segment_length",
            "SegmentLength must be an integer >= 4.",
        )
    L = int(L)
    if L > N:
        raise SidError(
            "segment_too_long",
            f"SegmentLength ({L}) exceeds data length ({N}).",
        )
    if Ts <= 0:
        raise SidError("invalid_sample_time", "SampleTime must be a positive scalar.")
    if alg not in ("bt", "welch"):
        raise SidError(
            "invalid_algorithm",
            f"Algorithm must be 'bt' or 'welch'. Got '{algorithm}'.",
        )

    # ---- Overlap default and validation ----
    P = overlap if overlap is not None else L // 2
    if not isinstance(P, (int, np.integer)) or P != int(P) or P < 0 or P >= L:
        raise SidError("invalid_overlap", "Overlap must be an integer in [0, L-1].")
    P = int(P)

    # ---- Algorithm-specific parameter setup ----
    M: int | None
    if alg == "bt":
        # BT-specific defaults and validation
        if window_size is None:
            M = min(L // 10, 30)
        else:
            M = int(window_size)
        if M < 2:
            raise SidError(
                "invalid_window_size",
                "WindowSize must be an integer >= 2.",
            )
        if L <= 2 * M:
            raise SidError(
                "segment_too_short",
                f"SegmentLength ({L}) must be greater than 2*WindowSize ({2 * M}).",
            )

        bt_kwargs: dict = {"window_size": M, "sample_time": Ts}
        if frequencies is not None:
            bt_kwargs["frequencies"] = frequencies
    else:
        # Welch-specific defaults and validation
        M = None
        Lsub = (
            int(sub_segment_length) if sub_segment_length is not None else int(math.floor(L / 4.5))
        )
        Psub = int(sub_overlap) if sub_overlap is not None else Lsub // 2
        nfft_val = int(nfft) if nfft is not None else max(256, 1 << math.ceil(math.log2(Lsub)))
        win_type = window

        if Lsub < 2:
            raise SidError(
                "invalid_sub_segment_length",
                "SubSegmentLength must be an integer >= 2.",
            )
        if Lsub > L:
            raise SidError(
                "sub_segment_too_long",
                f"SubSegmentLength ({Lsub}) exceeds SegmentLength ({L}).",
            )
        if Psub < 0 or Psub >= Lsub:
            raise SidError(
                "invalid_sub_overlap",
                "SubOverlap must be in [0, SubSegmentLength-1].",
            )

    # ---- Outer segmentation (SPEC.md S6.2) ----
    step = L - P
    K = (N - L) // step + 1
    if K < 1:
        raise SidError(
            "too_few_segments",
            f"Data too short for even one segment with L={L}, P={P}.",
        )

    # ---- Helper to extract a segment ----
    def _extract_segment(start: int) -> tuple[np.ndarray, np.ndarray | None]:
        end = start + L
        if n_traj > 1:
            yk = y[start:end, :, :]
            uk = None if is_time_series else u[start:end, :, :]
        else:
            yk = y[start:end, :]
            uk = None if is_time_series else u[start:end, :]
        return yk, uk

    # ---- Helper to run inner estimator on one segment ----
    def _estimate_segment(yk: np.ndarray, uk: np.ndarray | None) -> FreqResult:
        if alg == "bt":
            return freq_bt(yk, uk, **bt_kwargs)
        else:
            return _welch_estimate(
                yk,
                uk,
                is_time_series,
                ny,
                nu,
                Lsub,
                Psub,
                nfft_val,
                win_type,
                Ts,
            )

    # ---- First segment (to learn output dimensions) ----
    yk0, uk0 = _extract_segment(0)
    r1 = _estimate_segment(yk0, uk0)
    nf = len(r1.frequency)

    # ---- Pre-allocate arrays ----
    if is_time_series:
        G_all: np.ndarray | None = None
        GStd_all: np.ndarray | None = None
        Coh_all: np.ndarray | None = None
        if ny == 1:
            NS_all = np.zeros((nf, K))
            NSStd_all = np.zeros((nf, K))
        else:
            NS_all = np.zeros((nf, K, ny, ny))
            NSStd_all = np.zeros((nf, K, ny, ny))
    else:
        if ny == 1 and nu == 1:
            # SISO
            G_all = np.zeros((nf, K), dtype=complex)
            GStd_all = np.zeros((nf, K))
            NS_all = np.zeros((nf, K))
            NSStd_all = np.zeros((nf, K))
            if r1.coherence is not None:
                Coh_all = np.zeros((nf, K))
            else:
                Coh_all = None
        else:
            # MIMO
            G_all = np.zeros((nf, K, ny, nu), dtype=complex)
            GStd_all = np.zeros((nf, K, ny, nu))
            NS_all = np.zeros((nf, K, ny, ny))
            NSStd_all = np.zeros((nf, K, ny, ny))
            Coh_all = None

    # ---- Helper to store one segment result ----
    def _store_segment(idx: int, rk: FreqResult) -> None:
        if not is_time_series:
            if ny == 1 and nu == 1:
                G_all[:, idx] = rk.response.ravel()
                GStd_all[:, idx] = rk.response_std.ravel()
            else:
                G_all[:, idx, :, :] = rk.response.reshape(nf, ny, nu)
                GStd_all[:, idx, :, :] = rk.response_std.reshape(nf, ny, nu)
            if Coh_all is not None:
                Coh_all[:, idx] = rk.coherence

        if ny == 1 or is_time_series:
            NS_all[:, idx] = rk.noise_spectrum.ravel()
            NSStd_all[:, idx] = rk.noise_spectrum_std.ravel()
        else:
            NS_all[:, idx, :, :] = rk.noise_spectrum.reshape(nf, ny, ny)
            NSStd_all[:, idx, :, :] = rk.noise_spectrum_std.reshape(nf, ny, ny)

    # ---- Store first segment ----
    _store_segment(0, r1)

    # ---- Loop over remaining segments (SPEC.md S6.3) ----
    for k in range(1, K):
        start_idx = k * step
        yk, uk = _extract_segment(start_idx)
        rk = _estimate_segment(yk, uk)
        _store_segment(k, rk)

    # ---- Time vector (SPEC.md S6.2) ----
    # t_k = (k * step + L/2) * Ts -- centre of each segment
    time_vec = (np.arange(K, dtype=np.float64) * step + L / 2.0) * Ts

    # ---- Pack result ----
    return FreqMapResult(
        time=time_vec,
        frequency=r1.frequency,
        frequency_hz=r1.frequency_hz,
        response=G_all,
        response_std=GStd_all,
        noise_spectrum=NS_all,
        noise_spectrum_std=NSStd_all,
        coherence=Coh_all,
        sample_time=Ts,
        segment_length=L,
        overlap=P,
        window_size=M,
        algorithm=alg,
        num_trajectories=n_traj,
        method="freq_map",
    )
