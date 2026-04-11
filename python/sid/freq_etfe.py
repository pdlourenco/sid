# Copyright (c) 2026 Pedro Lourenco. All rights reserved.
# This code is released under the MIT License. See LICENSE file in the
# project root for full license information.
#
# This module is part of the Open Source System Identification Toolbox (SID).
# https://github.com/pdlourenco/sid

"""Empirical Transfer Function Estimate (ETFE) for frequency response estimation."""

from __future__ import annotations

import warnings

import numpy as np

from sid._exceptions import SidError
from sid._internal.dft import sid_dft
from sid._internal.is_default_freqs import is_default_freqs
from sid._internal.validate_data import validate_data
from sid._results import FreqResult


# ---------------------------------------------------------------------------
# Private helper
# ---------------------------------------------------------------------------


def _boxcar_smooth(x: np.ndarray, S: int) -> np.ndarray:
    """Apply length-*S* boxcar (moving average) smoothing.

    At boundaries the window shrinks to the available neighbors so
    that no zero-padding is needed and the output length equals the
    input length.

    Parameters
    ----------
    x : ndarray, shape ``(nf,)``
        1-D complex or real array to smooth.
    S : int
        Window length (positive odd integer).

    Returns
    -------
    ndarray, shape ``(nf,)``
        Smoothed array.
    """
    nf = len(x)
    out = np.empty_like(x)
    half = (S - 1) // 2
    for k in range(nf):
        lo = max(0, k - half)
        hi = min(nf, k + half + 1)
        out[k] = np.mean(x[lo:hi])
    return out


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def freq_etfe(
    y: np.ndarray | list,
    u: np.ndarray | list | None = None,
    *,
    smoothing: int = 1,
    frequencies: np.ndarray | None = None,
    sample_time: float = 1.0,
) -> FreqResult:
    """Estimate the frequency response via the Empirical Transfer Function Estimate.

    Computes the frequency response as the ratio of the output and input
    discrete Fourier transforms.  Provides maximum frequency resolution
    but high variance.  Optional smoothing reduces variance at the cost
    of resolution.

    This is an open-source replacement for the System Identification
    Toolbox function ``etfe``.

    Parameters
    ----------
    y : ndarray or list of ndarray, shape (N, ny) or (N, ny, L)
        Output data.  A 1-D array is treated as a single-channel signal.
        For multiple trajectories pass a 3-D array ``(N, ny, L)`` or a
        list of 1-D / 2-D arrays (variable-length trajectories are trimmed
        to the shortest).  Cross-periodograms are ensemble-averaged across
        trajectories.
    u : ndarray, list of ndarray, or None, optional
        Input data.  Same shape conventions as *y*.  Pass ``None`` for
        time-series mode (output periodogram only).  Default is ``None``.
    smoothing : int, optional
        Smoothing window length *S*.  Must be a positive integer; when
        *S* > 1 it must be odd.  A length-*S* boxcar (moving average)
        filter is applied to the raw ETFE.  Default is ``1`` (no
        smoothing).
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
          ``response``, same shape.  Filled with ``NaN`` (ETFE has no
          closed-form asymptotic variance).  ``None`` in time-series mode.
        - **noise_spectrum** (*ndarray*) -- Noise power spectrum (or
          output periodogram in time-series mode).  Shape ``(nf,)`` for
          SISO / single-channel time-series, ``(nf, ny, ny)`` for MIMO
          or multi-channel time-series.
        - **noise_spectrum_std** (*ndarray*) -- Standard deviation of
          ``noise_spectrum``, same shape.  Filled with ``NaN``.
        - **coherence** (*None*) -- Always ``None`` for ETFE.
        - **sample_time** (*float*) -- Sample time in seconds.
        - **window_size** (*int*) -- Data length *N*.
        - **data_length** (*int*) -- Number of samples *N* per trajectory.
        - **num_trajectories** (*int*) -- Number of trajectories.
        - **method** (*str*) -- ``'freq_etfe'``.

    Raises
    ------
    SidError
        If *sample_time* is not positive (code: ``'bad_ts'``).
    SidError
        If *smoothing* is not a positive integer, or is even and > 1
        (code: ``'bad_smoothing'``).
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
    >>> result = sid.freq_etfe(y, u, smoothing=5)  # doctest: +SKIP
    >>> result.response.shape  # doctest: +SKIP
    (128,)

    Time-series periodogram:

    >>> y = rng.standard_normal(500)  # doctest: +SKIP
    >>> result = sid.freq_etfe(y)  # doctest: +SKIP
    >>> result.noise_spectrum.shape  # doctest: +SKIP
    (128,)

    Custom frequencies:

    >>> w = np.linspace(0.01, np.pi, 256)  # doctest: +SKIP
    >>> result = sid.freq_etfe(y, u, frequencies=w)  # doctest: +SKIP

    Multi-trajectory (ensemble-averaged):

    >>> L = 5; N = 1000  # doctest: +SKIP
    >>> u3d = rng.standard_normal((N, 1, L))  # doctest: +SKIP
    >>> y3d = np.zeros_like(u3d)  # doctest: +SKIP
    >>> for l in range(L):  # doctest: +SKIP
    ...     y3d[:, 0, l] = lfilter([1], [1, -0.9], u3d[:, 0, l]) + 0.1 * rng.standard_normal(N)
    >>> result = sid.freq_etfe(y3d, u3d)  # doctest: +SKIP

    Notes
    -----
    **Algorithm:**

    1. Validate and orient input data; detect time-series vs. SISO vs. MIMO.
    2. Compute discrete Fourier transforms Y(w) and U(w) of the data.
    3. Form the raw ETFE: G(w) = Y(w) / U(w) (SISO) or
       G(w) = Phi_yu(w) * Phi_u(w)^{-1} (MIMO), using the H1 estimator
       for multi-trajectory data.
    4. Optionally smooth G with a length-*S* boxcar (moving average) window.
    5. Noise spectrum: Phi_v(w) = (1/N) * |Y(w) - G(w) * U(w)|^2.
       Time-series mode: periodogram Phi_y(w) = (1/N) * |Y(w)|^2.
    6. ETFE has no closed-form asymptotic variance; standard deviations
       are returned as NaN arrays.

    **Specification:** SPEC.md S4 -- Empirical Transfer Function Estimate

    References
    ----------
    .. [1] Ljung, L., "System Identification: Theory for the User",
       2nd ed., Prentice Hall, 1999. Sections 2.3, 6.3.

    See Also
    --------
    sid.freq_bt : Blackman-Tukey spectral analysis (lower variance).
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
    S = smoothing
    freqs = frequencies
    Ts = sample_time

    # ---- Validate parameters ----
    if Ts <= 0:
        raise SidError("bad_ts", "Sample time must be positive.")

    if S < 1 or S != round(S):
        raise SidError(
            "bad_smoothing",
            "Smoothing parameter S must be a positive integer.",
        )
    if S > 1 and S % 2 == 0:
        raise SidError(
            "bad_smoothing",
            "Smoothing parameter S must be odd.",
        )

    # ---- Frequency grid ----
    if freqs is None:
        freqs = np.arange(1, 129) * (np.pi / 128)
        use_fft = True
    else:
        freqs = np.asarray(freqs, dtype=np.float64).ravel()
        use_fft = is_default_freqs(freqs)

    nf = len(freqs)

    if np.any(freqs <= 0) or np.any(freqs > np.pi):
        raise SidError(
            "bad_freqs",
            "Frequencies must be in the range (0, pi] rad/sample.",
        )

    # ---- Compute DFTs (SPEC.md S4.1) ----
    # Y(w) = sum_{t=1}^{N} y(t) * exp(-jw*t), U(w) analogously.
    # sid_dft returns shape (nf, p).
    if n_traj == 1:
        Ydft = sid_dft(y, freqs, use_fft)  # (nf, ny)
        if not is_time_series:
            assert u is not None
            Udft = sid_dft(u, freqs, use_fft)  # (nf, nu)
    else:
        # Multi-trajectory: DFTs are computed per trajectory below.
        Ydft = None
        Udft = None

    # ---- Form transfer function and noise spectrum (SPEC.md S4.2-4.3) ----
    if is_time_series:
        # Periodogram: Phi_y(w) = (1/L) * sum_l (1/N) |Y_l(w)|^2
        G: np.ndarray | None = None

        if n_traj == 1:
            assert Ydft is not None
            if ny == 1:
                # (nf, 1) -> (nf,)
                PhiV = (1.0 / N) * np.abs(Ydft[:, 0]) ** 2
            else:
                PhiV = np.zeros((nf, ny, ny))
                for kk in range(nf):
                    Yk = Ydft[kk, :][:, np.newaxis]  # (ny, 1)
                    PhiV[kk, :, :] = (1.0 / N) * np.real(Yk @ Yk.conj().T)
        else:
            assert u is None
            if ny == 1:
                PhiV = np.zeros(nf)
                for ll in range(n_traj):
                    Yl = sid_dft(y[:, :, ll], freqs, use_fft)  # (nf, 1)
                    PhiV += (1.0 / N) * np.abs(Yl[:, 0]) ** 2
                PhiV /= n_traj
            else:
                PhiV = np.zeros((nf, ny, ny))
                for ll in range(n_traj):
                    Yl = sid_dft(y[:, :, ll], freqs, use_fft)  # (nf, ny)
                    for kk in range(nf):
                        Yk = Yl[kk, :][:, np.newaxis]  # (ny, 1)
                        PhiV[kk, :, :] += (1.0 / N) * np.real(Yk @ Yk.conj().T)
                PhiV /= n_traj

        Coh: np.ndarray | None = None

    elif ny == 1 and nu == 1:
        # ---- SISO: G(w) = Y(w) / U(w) ----
        eps_reg = 1e-10

        if n_traj == 1:
            assert Ydft is not None and Udft is not None
            # Single trajectory: direct ratio with regularization
            Y1 = Ydft[:, 0]  # (nf,)
            U1 = Udft[:, 0]  # (nf,)
            Uabs = np.abs(U1)
            Umax = np.max(Uabs)

            G_arr = np.empty(nf, dtype=np.complex128)
            for kk in range(nf):
                if Uabs[kk] < eps_reg * Umax:
                    G_arr[kk] = np.nan + 1j * np.nan
                else:
                    G_arr[kk] = Y1[kk] / U1[kk]
            G = G_arr
        else:
            # Multi-trajectory H1 estimator (SPEC.md S4.1)
            assert u is not None
            PhiYU = np.zeros(nf, dtype=np.complex128)
            PhiU = np.zeros(nf)
            for ll in range(n_traj):
                Yl = sid_dft(y[:, :, ll], freqs, use_fft)[:, 0]  # (nf,)
                Ul = sid_dft(u[:, :, ll], freqs, use_fft)[:, 0]  # (nf,)
                PhiYU += Yl * Ul.conj()
                PhiU += np.abs(Ul) ** 2
            PhiYU /= n_traj
            PhiU /= n_traj

            Umax = np.max(PhiU)
            G_arr = np.empty(nf, dtype=np.complex128)
            for kk in range(nf):
                if PhiU[kk] < eps_reg * Umax:
                    G_arr[kk] = np.nan + 1j * np.nan
                else:
                    G_arr[kk] = PhiYU[kk] / PhiU[kk]
            G = G_arr

        # Optional smoothing
        if S > 1:
            G = _boxcar_smooth(G, S)

        # Noise spectrum: Phi_v(w) = (1/N) * |Y(w) - G(w) * U(w)|^2
        if n_traj == 1:
            assert Ydft is not None and Udft is not None
            residual = Ydft[:, 0] - G * Udft[:, 0]
            PhiV = (1.0 / N) * np.abs(residual) ** 2
        else:
            assert u is not None
            PhiV = np.zeros(nf)
            for ll in range(n_traj):
                Yl = sid_dft(y[:, :, ll], freqs, use_fft)[:, 0]
                Ul = sid_dft(u[:, :, ll], freqs, use_fft)[:, 0]
                res_l = Yl - G * Ul
                PhiV += (1.0 / N) * np.abs(res_l) ** 2
            PhiV /= n_traj

        PhiV = np.maximum(PhiV, 0.0)
        Coh = None

    else:
        # ---- MIMO: G(w) = Phi_yu(w) * Phi_u(w)^{-1} ----
        assert u is not None
        eps_reg = 1e-10

        if n_traj == 1:
            assert Ydft is not None and Udft is not None
            # Single trajectory: rank-1 cross/auto-periodograms
            PhiYU = np.zeros((nf, ny, nu), dtype=np.complex128)
            PhiU = np.zeros((nf, nu, nu), dtype=np.complex128)
            for kk in range(nf):
                Yk = Ydft[kk, :][:, np.newaxis]  # (ny, 1)
                Uk = Udft[kk, :][:, np.newaxis]  # (nu, 1)
                PhiYU[kk, :, :] = Yk @ Uk.conj().T
                PhiU[kk, :, :] = Uk @ Uk.conj().T
        else:
            # Multi-trajectory: average cross-periodograms (SPEC.md S4.1)
            PhiYU = np.zeros((nf, ny, nu), dtype=np.complex128)
            PhiU = np.zeros((nf, nu, nu), dtype=np.complex128)
            for ll in range(n_traj):
                Yl = sid_dft(y[:, :, ll], freqs, use_fft)  # (nf, ny)
                Ul = sid_dft(u[:, :, ll], freqs, use_fft)  # (nf, nu)
                for kk in range(nf):
                    Yk = Yl[kk, :][:, np.newaxis]  # (ny, 1)
                    Uk = Ul[kk, :][:, np.newaxis]  # (nu, 1)
                    PhiYU[kk, :, :] += Yk @ Uk.conj().T
                    PhiU[kk, :, :] += Uk @ Uk.conj().T
            PhiYU /= n_traj
            PhiU /= n_traj

        G_mimo = np.zeros((nf, ny, nu), dtype=np.complex128)
        warned_singular = False
        for kk in range(nf):
            PhiU_k = PhiU[kk, :, :].reshape(nu, nu)
            PhiYU_k = PhiYU[kk, :, :].reshape(ny, nu)

            if nu == 1:
                # Scalar input spectrum
                if np.abs(PhiU_k[0, 0]) < eps_reg * np.max(np.abs(PhiU)):
                    G_mimo[kk, :, :] = np.nan
                else:
                    G_mimo[kk, :, :] = PhiYU_k / PhiU_k[0, 0]
            else:
                # MIMO: check condition number
                try:
                    rc = 1.0 / np.linalg.cond(PhiU_k)
                except (np.linalg.LinAlgError, FloatingPointError):
                    rc = 0.0
                if rc < eps_reg:
                    G_mimo[kk, :, :] = np.nan
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
                    G_mimo[kk, :, :] = np.linalg.solve(PhiU_k.T, PhiYU_k.T).T

        G = G_mimo

        # Optional smoothing (element-wise)
        if S > 1:
            for ii in range(ny):
                for jj in range(nu):
                    G[:, ii, jj] = _boxcar_smooth(G[:, ii, jj], S)

        # Noise spectrum: Phi_v(w) = (1/N) * res @ res.conj().T
        PhiV = np.zeros((nf, ny, ny), dtype=np.complex128)
        if n_traj == 1:
            assert Ydft is not None and Udft is not None
            for kk in range(nf):
                Yk = Ydft[kk, :][:, np.newaxis]  # (ny, 1)
                Uk = Udft[kk, :][:, np.newaxis]  # (nu, 1)
                Gk = G[kk, :, :].reshape(ny, nu)
                res = Yk - Gk @ Uk
                PhiV[kk, :, :] = (1.0 / N) * (res @ res.conj().T)
        else:
            for ll in range(n_traj):
                Yl = sid_dft(y[:, :, ll], freqs, use_fft)  # (nf, ny)
                Ul = sid_dft(u[:, :, ll], freqs, use_fft)  # (nf, nu)
                for kk in range(nf):
                    Yk = Yl[kk, :][:, np.newaxis]
                    Uk = Ul[kk, :][:, np.newaxis]
                    Gk = G[kk, :, :].reshape(ny, nu)
                    res = Yk - Gk @ Uk
                    PhiV[kk, :, :] += (1.0 / N) * (res @ res.conj().T)
            PhiV /= n_traj

        PhiV = np.real(PhiV)
        Coh = None

    # ---- Uncertainty (SPEC.md S4.5) ----
    # ETFE has no closed-form asymptotic variance -- return NaN.
    if is_time_series:
        GStd: np.ndarray | None = None
        PhiVStd = np.full_like(PhiV, np.nan)
    else:
        assert G is not None
        GStd = np.full_like(G, np.nan, dtype=np.float64)
        PhiVStd = np.full_like(PhiV, np.nan, dtype=np.float64)

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
        window_size=N,
        data_length=N,
        num_trajectories=n_traj,
        method="freq_etfe",
    )
