# Copyright (c) 2026 Pedro Lourenco. All rights reserved.
# This code is released under the MIT License. See LICENSE file in the
# project root for full license information.
#
# This module is part of the Open Source System Identification Toolbox (SID).
# https://github.com/pdlourenco/sid-matlab

"""Blackman-Tukey spectral analysis for frequency response estimation."""

from __future__ import annotations

import warnings

import numpy as np

from sid._exceptions import SidError
from sid._internal.cov import sid_cov
from sid._internal.hann_win import hann_win
from sid._internal.is_default_freqs import is_default_freqs
from sid._internal.uncertainty import sid_uncertainty
from sid._internal.validate_data import validate_data
from sid._internal.windowed_dft import windowed_dft
from sid._results import FreqResult


def freq_bt(
    y: np.ndarray | list,
    u: np.ndarray | list | None = None,
    *,
    window_size: int | None = None,
    frequencies: np.ndarray | None = None,
    sample_time: float = 1.0,
) -> FreqResult:
    """Estimate frequency response via Blackman-Tukey spectral analysis.

    Estimates the frequency response G(e^{jw}) and noise power spectrum
    from time-domain input/output data using the Blackman-Tukey method.
    The algorithm computes biased sample covariances, applies a Hann lag
    window, and Fourier-transforms to obtain spectral estimates.

    This is an open-source replacement for the System Identification
    Toolbox function ``spa``.

    Parameters
    ----------
    y : ndarray or list of ndarray, shape (N, ny) or (N, ny, L)
        Output data.  A 1-D array is treated as a single-channel signal.
        For multiple trajectories pass a 3-D array ``(N, ny, L)`` or a
        list of 1-D / 2-D arrays (variable-length trajectories are trimmed
        to the shortest).  Spectral estimates are ensemble-averaged across
        trajectories.
    u : ndarray, list of ndarray, or None, optional
        Input data.  Same shape conventions as *y*.  Pass ``None`` for
        time-series mode (output spectrum only).  Default is ``None``.
    window_size : int, optional
        Hann lag window size *M*.  Must satisfy ``M >= 2``.  If *M*
        exceeds ``N // 2`` it is silently reduced.
        Default is ``min(N // 10, 30)``.
    frequencies : ndarray, shape (nf,), optional
        Frequency vector in rad/sample.  All values must lie in the
        interval ``(0, pi]``.  Default is 128 linearly spaced values
        ``k * pi / 128`` for ``k = 1, ..., 128``.
    sample_time : float, optional
        Sample time in seconds.  Must be positive.  Default is ``1.0``.

    Returns
    -------
    FreqResult
        Frozen dataclass with fields:

        - **frequency** (*ndarray, shape (nf,)*) -- Frequency vector,
          rad/sample.
        - **frequency_hz** (*ndarray, shape (nf,)*) -- Frequency vector,
          Hz.
        - **response** (*ndarray or None*) -- Complex frequency response.
          Shape ``(nf,)`` for SISO, ``(nf, ny, nu)`` for MIMO, or
          ``None`` in time-series mode.
        - **response_std** (*ndarray or None*) -- Standard deviation of
          ``response``, same shape.  ``None`` in time-series mode.
        - **noise_spectrum** (*ndarray*) -- Noise power spectrum (or
          output spectrum in time-series mode).  Shape ``(nf,)`` for
          SISO / time-series, ``(nf, ny, ny)`` for MIMO.
        - **noise_spectrum_std** (*ndarray*) -- Standard deviation of
          ``noise_spectrum``, same shape.
        - **coherence** (*ndarray or None*) -- Squared coherence, shape
          ``(nf,)``.  SISO only; ``None`` for MIMO or time-series.
        - **sample_time** (*float*) -- Sample time in seconds.
        - **window_size** (*int*) -- Lag window size *M* used.
        - **data_length** (*int*) -- Number of samples *N* per trajectory.
        - **num_trajectories** (*int*) -- Number of trajectories.
        - **method** (*str*) -- ``'freq_bt'``.

    Raises
    ------
    SidError
        If *sample_time* is not positive (code: ``'bad_ts'``).
    SidError
        If *window_size* is less than 2 (code: ``'bad_window_size'``).
    SidError
        If any frequency is outside ``(0, pi]`` (code: ``'bad_freqs'``).
    SidError
        If data contains NaN/Inf (code: ``'non_finite'``), is complex
        (code: ``'complex_data'``), or is too short
        (code: ``'too_short'``).  These are raised by the data
        validation layer.

    Examples
    --------
    Basic SISO system identification:

    >>> import numpy as np
    >>> import sid  # doctest: +SKIP
    >>> N = 1000; rng = np.random.default_rng(0)
    >>> u = rng.standard_normal(N)
    >>> from scipy.signal import lfilter  # doctest: +SKIP
    >>> y = lfilter([1], [1, -0.9], u) + 0.1 * rng.standard_normal(N)  # doctest: +SKIP
    >>> result = sid.freq_bt(y, u)  # doctest: +SKIP
    >>> result.response.shape  # doctest: +SKIP
    (128,)

    Time-series spectrum estimation:

    >>> y = rng.standard_normal(500)  # doctest: +SKIP
    >>> result = sid.freq_bt(y)  # doctest: +SKIP
    >>> result.noise_spectrum.shape  # doctest: +SKIP
    (128,)

    Custom window size and frequencies:

    >>> w = np.linspace(0.01, np.pi, 256)  # doctest: +SKIP
    >>> result = sid.freq_bt(y, u, window_size=50, frequencies=w)  # doctest: +SKIP

    Multi-trajectory (ensemble-averaged):

    >>> L = 5; N = 1000  # doctest: +SKIP
    >>> u3d = rng.standard_normal((N, 1, L))  # doctest: +SKIP
    >>> y3d = np.zeros_like(u3d)  # doctest: +SKIP
    >>> for l in range(L):  # doctest: +SKIP
    ...     y3d[:, 0, l] = lfilter([1], [1, -0.9], u3d[:, 0, l]) + 0.1 * rng.standard_normal(N)
    >>> result = sid.freq_bt(y3d, u3d)  # doctest: +SKIP

    Notes
    -----
    **Algorithm:**

    1. Validate and orient input data; detect time-series vs. SISO vs. MIMO.
    2. Compute biased sample covariances R_yy, R_uu, R_yu for lags
       0 .. M (ensemble-averaged for multi-trajectory data).
    3. Apply Hann lag window and Fourier-transform to obtain spectral
       estimates Phi_y, Phi_u, Phi_yu.  Uses the FFT fast path when the
       frequency grid matches the default 128-point linear grid, and
       direct DFT otherwise.
    4. Form the transfer function estimate G = Phi_yu / Phi_u (SISO) or
       G = Phi_yu * Phi_u^{-1} (MIMO), and the noise spectrum
       Phi_v = Phi_y - |Phi_yu|^2 / Phi_u (SISO) or
       Phi_v = Phi_y - Phi_yu * Phi_u^{-1} * Phi_yu' (MIMO).
    5. Compute asymptotic standard deviations using Ljung (1999) formulae.

    **Specification:** SPEC.md S2 -- Blackman-Tukey Spectral Analysis

    References
    ----------
    .. [1] Ljung, L., "System Identification: Theory for the User",
       2nd ed., Prentice Hall, 1999. Sections 2.3, 6.3--6.4.

    See Also
    --------
    sid.freq_etfe : Empirical Transfer Function Estimate (no lag window).
    sid.freq_btfdr : Blackman-Tukey with frequency-dependent resolution.
    sid.bode_plot : Bode magnitude/phase plot of a FreqResult.
    sid.spectrum_plot : Plot the noise/output spectrum of a FreqResult.

    Changelog
    ---------
    2026-04-08 : First version (Python port) by Pedro Lourenco.
    """
    # ---- Validate data ----
    y, u, N, ny, nu, is_time_series, n_traj = validate_data(y, u)

    # ---- Apply defaults ----
    M = window_size
    freqs = frequencies
    Ts = sample_time

    if M is None:
        M = min(N // 10, 30)
    if freqs is None:
        freqs = np.arange(1, 129) * (np.pi / 128)
    else:
        freqs = np.asarray(freqs, dtype=np.float64).ravel()

    # ---- Validate parameters ----
    if Ts <= 0:
        raise SidError("bad_ts", "Sample time must be positive.")
    if M < 2:
        raise SidError("bad_window_size", "Window size M must be at least 2.")
    if M > N // 2:
        warnings.warn(
            f"Window size {M} exceeds N/2 = {N // 2}. Reduced to {N // 2}.",
            stacklevel=2,
        )
        M = N // 2
    if np.any(freqs <= 0) or np.any(freqs > np.pi):
        raise SidError(
            "bad_freqs",
            "Frequencies must be in the range (0, pi] rad/sample.",
        )

    nf = len(freqs)
    # Use FFT fast path when frequencies match the default 128-point grid
    use_fft = (frequencies is None) or is_default_freqs(freqs)

    # ---- Compute biased sample covariances (SPEC.md S2.3) ----
    # Ryy corresponds to R-hat_yy(tau) for lags tau = 0..M
    Ryy = sid_cov(y, y, M)  # (M+1, ny, ny) or (M+1,) for scalar

    if not is_time_series:
        Ruu = sid_cov(u, u, M)  # (M+1, nu, nu) or (M+1,)
        Ryu = sid_cov(y, u, M)  # (M+1, ny, nu) or (M+1,)
        Ruy = sid_cov(u, y, M)  # (M+1, nu, ny) or (M+1,) -- negative lags

    # ---- Hann lag window (SPEC.md S2.4) ----
    W = hann_win(M)  # (M+1,) for lags 0..M

    # ---- Windowed DFT: covariances -> spectral estimates (SPEC.md S2.5) ----
    # Phi_y(w) = sum_{tau=-M}^{M} W(|tau|) * R_yy(tau) * e^{-jw*tau}
    PhiY = windowed_dft(Ryy, W, freqs, use_fft, Ryy)  # (nf, ny, ny) or (nf,)

    if not is_time_series:
        PhiU = windowed_dft(Ruu, W, freqs, use_fft, Ruu)  # (nf, nu, nu) or (nf,)
        PhiYU = windowed_dft(Ryu, W, freqs, use_fft, Ruy)  # (nf, ny, nu) or (nf,)

    # ---- Form transfer function and noise spectrum (SPEC.md S2.6-2.7) ----
    if is_time_series:
        G: np.ndarray | None = None
        PhiV = np.real(PhiY)
        Coh: np.ndarray | None = None

    elif ny == 1 and nu == 1:
        # SISO: G(w) = Phi_yu(w) / Phi_u(w)  (SPEC.md S2.6)
        # Regularization: if |Phi_u(w_k)| < eps * max(|Phi_u|), set G = NaN
        eps_reg = 1e-10
        PhiU_abs = np.abs(np.real(PhiU))
        PhiU_max = float(np.max(PhiU_abs)) if PhiU_abs.size > 0 else 0.0
        singular_mask = PhiU_abs < eps_reg * PhiU_max

        G = np.empty_like(PhiYU)
        with np.errstate(divide="ignore", invalid="ignore"):
            G[:] = PhiYU / PhiU
        if np.any(singular_mask):
            G[singular_mask] = np.nan + 1j * np.nan
            warnings.warn(
                "Input spectrum Phi_u is near-singular at some "
                "frequencies. G set to NaN at those points.",
                stacklevel=2,
            )

        # Phi_v(w) = Phi_y(w) - |Phi_yu(w)|^2 / Phi_u(w)  -- noise spectrum
        with np.errstate(divide="ignore", invalid="ignore"):
            PhiV = np.real(PhiY) - np.abs(PhiYU) ** 2 / np.real(PhiU)
        if np.any(singular_mask):
            PhiV[singular_mask] = np.real(PhiY[singular_mask])
        PhiV = np.maximum(PhiV, 0.0)

        # gamma^2(w) = |Phi_yu|^2 / (Phi_y * Phi_u)  -- squared coherence
        with np.errstate(divide="ignore", invalid="ignore"):
            Coh = np.abs(PhiYU) ** 2 / (np.real(PhiY) * np.real(PhiU))
        if np.any(singular_mask):
            Coh[singular_mask] = 0.0
        Coh = np.clip(Coh, 0.0, 1.0)

    else:
        # MIMO: G(w) = Phi_yu(w) * Phi_u(w)^{-1} per frequency (SPEC.md S2.6)
        # Ensure spectral arrays are 3D even when one dimension is 1
        # (sid_cov squeezes scalar signals to 1D)
        if PhiY.ndim == 1:
            PhiY = PhiY[:, np.newaxis, np.newaxis]
        if PhiU.ndim == 1:
            PhiU = PhiU[:, np.newaxis, np.newaxis]
        if PhiYU.ndim == 1:
            PhiYU = PhiYU[:, np.newaxis, np.newaxis]

        G = np.zeros((nf, ny, nu), dtype=complex)
        PhiV = np.zeros((nf, ny, ny))
        eps_reg = 1e-10
        warned_singular = False

        for k in range(nf):
            PhiU_k = PhiU[k, :, :].reshape(nu, nu)  # (nu, nu)
            PhiYU_k = PhiYU[k, :, :].reshape(ny, nu)  # (ny, nu)
            PhiY_k = PhiY[k, :, :].reshape(ny, ny)  # (ny, ny)

            # Regularization: check condition of Phi_u (SPEC.md S2.6)
            rc = 1.0 / np.linalg.cond(PhiU_k)
            if rc < eps_reg:
                G[k, :, :] = np.nan
                PhiV[k, :, :] = np.real(PhiY_k)
                if not warned_singular:
                    warnings.warn(
                        "Input spectrum Phi_u is near-singular at some "
                        "frequencies. G set to NaN at those points.",
                        stacklevel=2,
                    )
                    warned_singular = True
            else:
                # MATLAB: PhiYU_k / PhiU_k  == PhiYU_k * inv(PhiU_k)
                # Python: solve(PhiU_k.T, PhiYU_k.T).T
                G[k, :, :] = np.linalg.solve(PhiU_k.T, PhiYU_k.T).T

                # Phi_v = Phi_y - Phi_yu * Phi_u^{-1} * Phi_yu'
                Gk = G[k, :, :]
                PhiV[k, :, :] = np.real(PhiY_k - Gk @ PhiYU_k.conj().T)

        PhiV = np.real(PhiV)

        # Clamp MIMO noise spectrum to PSD (SPEC.md S2.7):
        # zero any negative eigenvalues at each frequency.
        for k in range(nf):
            Vk = PhiV[k, :, :].reshape(ny, ny)
            Vk = (Vk + Vk.T) / 2  # enforce symmetry
            eigvals, eigvecs = np.linalg.eigh(Vk)
            if np.any(eigvals < 0):
                eigvals = np.maximum(eigvals, 0.0)
                Vk = eigvecs @ np.diag(eigvals) @ eigvecs.T
                PhiV[k, :, :] = np.real(Vk)

        Coh = None

    # ---- Asymptotic uncertainty (SPEC.md S3) ----
    if is_time_series:
        _, PhiVStd = sid_uncertainty(None, PhiV, None, N, W, n_traj)
        GStd: np.ndarray | None = None
    elif ny == 1 and nu == 1:
        GStd, PhiVStd = sid_uncertainty(G, PhiV, Coh, N, W, n_traj)
    else:
        # MIMO: pass PhiU for diagonal uncertainty approximation
        GStd, PhiVStd = sid_uncertainty(G, PhiV, Coh, N, W, n_traj, PhiU)

    # ---- Pack result ----
    return FreqResult(
        frequency=freqs,
        frequency_hz=freqs / (2.0 * np.pi * Ts),
        response=G,
        response_std=GStd,
        noise_spectrum=PhiV,
        noise_spectrum_std=PhiVStd,
        coherence=Coh,
        sample_time=Ts,
        window_size=M,
        data_length=N,
        num_trajectories=n_traj,
        method="freq_bt",
    )
