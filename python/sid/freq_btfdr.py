# Copyright (c) 2026 Pedro Lourenco. All rights reserved.
# This code is released under the MIT License. See LICENSE file in the
# project root for full license information.
#
# This module is part of the Open Source System Identification Toolbox (SID).
# https://github.com/pdlourenco/sid

"""Blackman-Tukey spectral analysis with frequency-dependent resolution."""

from __future__ import annotations


import numpy as np

from sid._exceptions import SidError
from sid._internal.cov import sid_cov
from sid._internal.hann_win import hann_win
from sid._internal.validate_data import validate_data
from sid._results import FreqResult


def freq_btfdr(
    y: np.ndarray | list,
    u: np.ndarray | list | None = None,
    *,
    resolution: float | np.ndarray | None = None,
    frequencies: np.ndarray | None = None,
    sample_time: float = 1.0,
) -> FreqResult:
    """Estimate frequency response via Blackman-Tukey with frequency-dependent resolution.

    Like :func:`sid.freq_bt`, but the Hann lag-window size varies across
    frequencies.  The user specifies a *resolution* parameter (in
    rad/sample) instead of a fixed window size.  Finer resolution
    (smaller value) uses a larger window and gives lower variance but
    coarser frequency detail, while coarser resolution (larger value)
    uses a smaller window.

    This is an open-source replacement for the System Identification
    Toolbox function ``spafdr``.

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
    resolution : float, ndarray, or None, optional
        Frequency resolution in rad/sample.  A scalar value applies
        uniformly to all frequencies; a vector of the same length as
        *frequencies* sets per-frequency resolution.  Must be positive.
        Default is ``2 * pi / min(N // 10, 30)`` (clamped so the
        implied window size is at least 2).
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
        - **window_size** (*ndarray, shape (nf,)*) -- Per-frequency lag
          window sizes ``M_k``.
        - **data_length** (*int*) -- Number of samples *N* per trajectory.
        - **num_trajectories** (*int*) -- Number of trajectories.
        - **method** (*str*) -- ``'freq_btfdr'``.

    Raises
    ------
    SidError
        If *sample_time* is not positive (code: ``'bad_ts'``).
    SidError
        If any *resolution* value is not positive
        (code: ``'bad_resolution'``).
    SidError
        If *resolution* is a vector whose length does not match the
        frequency grid (code: ``'bad_resolution'``).
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
    >>> result = sid.freq_btfdr(y, u, resolution=0.3)  # doctest: +SKIP
    >>> result.response.shape  # doctest: +SKIP
    (128,)

    Time-series spectrum estimation:

    >>> y = rng.standard_normal(500)  # doctest: +SKIP
    >>> result = sid.freq_btfdr(y)  # doctest: +SKIP
    >>> result.noise_spectrum.shape  # doctest: +SKIP
    (128,)

    Per-frequency resolution vector:

    >>> w = np.linspace(0.01, np.pi, 64)  # doctest: +SKIP
    >>> R = np.linspace(0.1, 1.0, 64)  # doctest: +SKIP
    >>> result = sid.freq_btfdr(y, u, resolution=R, frequencies=w)  # doctest: +SKIP

    Multi-trajectory (ensemble-averaged):

    >>> L = 5; N = 1000  # doctest: +SKIP
    >>> u3d = rng.standard_normal((N, 1, L))  # doctest: +SKIP
    >>> y3d = np.zeros_like(u3d)  # doctest: +SKIP
    >>> for l in range(L):  # doctest: +SKIP
    ...     y3d[:, 0, l] = lfilter([1], [1, -0.9], u3d[:, 0, l]) + 0.1 * rng.standard_normal(N)
    >>> result = sid.freq_btfdr(y3d, u3d, resolution=0.3)  # doctest: +SKIP

    Notes
    -----
    **Algorithm:**

    1. Validate and orient input data; detect time-series vs. SISO vs. MIMO.
    2. Convert the resolution parameter to per-frequency window sizes
       ``M_k = ceil(2 * pi / R_k)``, clamped to ``[2, N // 2]``.
    3. Pre-compute biased sample covariances up to ``max(M_k)`` using
       ensemble averaging for multi-trajectory data.
    4. For each frequency ``w_k``:

       a. Truncate covariances to lag ``M_k`` and compute a local Hann
          window of size ``M_k``.
       b. Compute windowed spectral estimates via a direct single-frequency
          DFT (not FFT, since the window size varies per frequency).
       c. Form G, noise spectrum, coherence, and asymptotic uncertainty
          using the local window norm ``C_W``.

    The per-frequency approach gives smooth bias-variance trade-offs
    across the frequency axis, at the cost of a slower per-frequency
    loop compared to the fixed-window :func:`sid.freq_bt`.

    **Specification:** SPEC.md S5 -- Frequency-Dependent Resolution

    References
    ----------
    .. [1] Ljung, L., "System Identification: Theory for the User",
       2nd ed., Prentice Hall, 1999. Sections 6.3--6.4.

    See Also
    --------
    sid.freq_bt : Blackman-Tukey with fixed window size.
    sid.freq_etfe : Empirical Transfer Function Estimate (no lag window).
    sid.bode_plot : Bode magnitude/phase plot of a FreqResult.
    sid.spectrum_plot : Plot the noise/output spectrum of a FreqResult.

    Changelog
    ---------
    2026-04-08 : First version (Python port) by Pedro Lourenco.
    """
    # ---- Validate data ----
    y, u, N, ny, nu, is_time_series, n_traj = validate_data(y, u)
    Neff = N * n_traj  # effective sample size for variance scaling

    # ---- Apply defaults ----
    R = resolution
    freqs = frequencies
    Ts = sample_time

    if freqs is None:
        freqs = np.arange(1, 129) * (np.pi / 128)
    else:
        freqs = np.asarray(freqs, dtype=np.float64).ravel()

    nf = len(freqs)

    if R is None:
        M_default = max(min(N // 10, 30), 2)
        R = 2.0 * np.pi / M_default

    # ---- Validate parameters ----
    if Ts <= 0:
        raise SidError("bad_ts", "Sample time must be positive.")

    if np.any(freqs <= 0) or np.any(freqs > np.pi):
        raise SidError(
            "bad_freqs",
            "Frequencies must be in the range (0, pi] rad/sample.",
        )

    R = np.atleast_1d(np.asarray(R, dtype=np.float64))
    if np.any(R <= 0):
        raise SidError("bad_resolution", "Resolution must be positive.")

    # Expand scalar R to vector
    if R.size == 1:
        R = np.full(nf, R[0])
    else:
        R = R.ravel()
        if len(R) != nf:
            raise SidError(
                "bad_resolution",
                f"Resolution vector length ({len(R)}) must match frequency vector length ({nf}).",
            )

    # ---- Resolution to window size (SPEC.md S5.2) ----
    Mk = np.clip(np.ceil(2.0 * np.pi / R).astype(int), 2, N // 2)

    # ---- Pre-compute biased covariances up to max(Mk) (SPEC.md S2.3) ----
    Mmax = int(np.max(Mk))
    Ryy = sid_cov(y, y, Mmax)  # (Mmax+1, ny, ny) or (Mmax+1,) for scalar

    if not is_time_series:
        Ruu = sid_cov(u, u, Mmax)  # (Mmax+1, nu, nu) or (Mmax+1,)
        Ryu = sid_cov(y, u, Mmax)  # (Mmax+1, ny, nu) or (Mmax+1,)
        Ruy = sid_cov(u, y, Mmax)  # (Mmax+1, nu, ny) or (Mmax+1,) neg lags

    # ---- Per-frequency spectral estimation (SPEC.md S5.2) ----
    is_siso = ny == 1 and nu == 1 and not is_time_series

    if is_time_series:
        # ---- Time-series path ----
        G: np.ndarray | None = None
        GStd: np.ndarray | None = None
        Coh: np.ndarray | None = None

        PhiV = np.zeros((nf, ny, ny))
        PhiVStd = np.zeros((nf, ny, ny))

        for k in range(nf):
            Mk_k = Mk[k]
            W = hann_win(Mk_k)
            Ryy_k = _truncate_cov(Ryy, Mk_k)

            # Direct DFT at single frequency
            PhiY_k = _matrix_single_freq_dft(Ryy_k, W, freqs[k], ny, ny)
            PhiV[k, :, :] = np.real(PhiY_k)

            # Uncertainty: Var{Phi_y} = (2*C_W/Neff) * Phi_y^2
            CW = float(W[0] ** 2 + 2.0 * np.sum(W[1:] ** 2))
            PhiVStd[k, :, :] = np.sqrt(2.0 * CW / Neff) * np.abs(PhiY_k)

        # Squeeze if scalar output
        if ny == 1:
            PhiV = PhiV.ravel()
            PhiVStd = PhiVStd.ravel()

    elif is_siso:
        # ---- SISO path ----
        # Covariances are 1-D (scalar signals)
        G_arr = np.zeros(nf, dtype=complex)
        PhiV_arr = np.zeros(nf)
        GStd_arr = np.zeros(nf)
        PhiVStd_arr = np.zeros(nf)
        Coh_arr = np.zeros(nf)

        eps_reg = 1e-10

        # Pass 1: compute all spectral estimates
        PhiY_all = np.zeros(nf, dtype=complex)
        PhiU_all = np.zeros(nf, dtype=complex)
        PhiYU_all = np.zeros(nf, dtype=complex)
        W_store: list[np.ndarray] = []

        for k in range(nf):
            Mk_k = Mk[k]
            W = hann_win(Mk_k)
            W_store.append(W)

            PhiY_all[k] = _scalar_single_freq_dft(
                Ryy[: Mk_k + 1],
                W,
                freqs[k],
            )
            PhiU_all[k] = _scalar_single_freq_dft(
                Ruu[: Mk_k + 1],
                W,
                freqs[k],
            )
            PhiYU_all[k] = _scalar_single_freq_dft(
                Ryu[: Mk_k + 1],
                W,
                freqs[k],
                R_neg=Ruy[: Mk_k + 1],
            )

        PhiU_max = np.max(np.abs(PhiU_all))

        # Pass 2: form G, PhiV, uncertainty
        for k in range(nf):
            W = W_store[k]
            PhiY_k = PhiY_all[k]
            PhiU_k = PhiU_all[k]
            PhiYU_k = PhiYU_all[k]

            if np.abs(PhiU_k) < eps_reg * PhiU_max:
                G_arr[k] = np.nan + 1j * np.nan
                PhiV_arr[k] = np.real(PhiY_k)
                Coh_arr[k] = 0.0
                GStd_arr[k] = np.inf
            else:
                G_arr[k] = PhiYU_k / PhiU_k
                PhiV_arr[k] = max(
                    np.real(PhiY_k) - np.abs(PhiYU_k) ** 2 / np.real(PhiU_k),
                    0.0,
                )
                Coh_arr[k] = float(
                    np.clip(
                        np.abs(PhiYU_k) ** 2 / (np.real(PhiY_k) * np.real(PhiU_k)),
                        0.0,
                        1.0,
                    )
                )

                # Uncertainty with local window norm (SPEC.md S5.3)
                CW = float(W[0] ** 2 + 2.0 * np.sum(W[1:] ** 2))
                coh_safe = max(Coh_arr[k], eps_reg)
                GStd_arr[k] = np.sqrt(
                    (CW / Neff) * np.abs(G_arr[k]) ** 2 * (1.0 - coh_safe) / coh_safe
                )

            # Noise uncertainty
            CW = float(W[0] ** 2 + 2.0 * np.sum(W[1:] ** 2))
            PhiVStd_arr[k] = np.sqrt(2.0 * CW / Neff) * np.abs(PhiV_arr[k])

        G = G_arr
        PhiV = PhiV_arr
        GStd = GStd_arr
        PhiVStd = PhiVStd_arr
        Coh = Coh_arr

    else:
        # ---- MIMO path ----
        G_mat = np.zeros((nf, ny, nu), dtype=complex)
        PhiV_mat = np.zeros((nf, ny, ny))
        GStd_mat = np.zeros((nf, ny, nu))
        PhiVStd_mat = np.zeros((nf, ny, ny))
        Coh = None

        eps_floor = 1e-10

        for k in range(nf):
            Mk_k = Mk[k]
            W = hann_win(Mk_k)

            Ryy_k = _truncate_cov(Ryy, Mk_k)
            Ruu_k = _truncate_cov(Ruu, Mk_k)
            Ryu_k = _truncate_cov(Ryu, Mk_k)
            Ruy_k = _truncate_cov(Ruy, Mk_k)

            PhiY_k = _matrix_single_freq_dft(Ryy_k, W, freqs[k], ny, ny)
            PhiU_k = _matrix_single_freq_dft(Ruu_k, W, freqs[k], nu, nu)
            PhiYU_k = _matrix_single_freq_dft(
                Ryu_k,
                W,
                freqs[k],
                ny,
                nu,
                R_neg=Ruy_k,
            )

            # G(w_k) = Phi_yu * Phi_u^{-1}  (MATLAB: PhiYU_k / PhiU_k)
            G_mat[k, :, :] = np.linalg.solve(PhiU_k.T, PhiYU_k.T).T

            # Phi_v = Phi_y - Phi_yu * Phi_u^{-1} * Phi_yu'
            Gk = G_mat[k, :, :]
            PhiV_k = PhiY_k - Gk @ PhiYU_k.conj().T
            PhiV_mat[k, :, :] = np.real(PhiV_k)

            # Noise uncertainty (SPEC.md S5.3)
            CW = float(W[0] ** 2 + 2.0 * np.sum(W[1:] ** 2))
            PhiVStd_mat[k, :, :] = np.sqrt(2.0 * CW / Neff) * np.abs(PhiV_k)

            # Diagonal MIMO G uncertainty:
            # Var{G_{ij}} ~ C_W / Neff * Phi_v_{ii} / Phi_u_{jj}
            for ii in range(ny):
                phiV_ii = max(np.real(PhiV_k[ii, ii]), 0.0)
                for jj in range(nu):
                    phiU_jj = np.real(PhiU_k[jj, jj])
                    if phiU_jj > eps_floor:
                        GStd_mat[k, ii, jj] = np.sqrt(CW / Neff * phiV_ii / phiU_jj)
                    else:
                        GStd_mat[k, ii, jj] = np.inf

        G = G_mat
        PhiV = PhiV_mat
        GStd = GStd_mat
        PhiVStd = PhiVStd_mat

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
        window_size=Mk,
        data_length=N,
        num_trajectories=n_traj,
        method="freq_btfdr",
    )


# ---------------------------------------------------------------------------
# Private helpers
# ---------------------------------------------------------------------------


def _truncate_cov(R: np.ndarray, Mk: int) -> np.ndarray:
    """Extract covariances for lags ``0 .. Mk`` from a larger array.

    Parameters
    ----------
    R : ndarray
        Full covariance array, shape ``(Mmax+1,)`` or ``(Mmax+1, p, q)``.
    Mk : int
        Local window size (number of non-zero lags).

    Returns
    -------
    ndarray
        Truncated covariance, shape ``(Mk+1,)`` or ``(Mk+1, p, q)``.
    """
    if R.ndim == 1:
        return R[: Mk + 1]
    return R[: Mk + 1, :, :]


def _scalar_single_freq_dft(
    R: np.ndarray,
    W: np.ndarray,
    w: float,
    R_neg: np.ndarray | None = None,
) -> complex:
    """Windowed DFT at a single frequency for a scalar covariance sequence.

    Computes:

    .. math::

        \\Phi(\\omega) = W_0 R_0 + \\sum_{\\tau=1}^{M}
            W_\\tau \\bigl(R_\\tau e^{-j\\omega\\tau}
            + R^{-}_\\tau e^{+j\\omega\\tau}\\bigr)

    where :math:`R^{-}_\\tau = \\overline{R_\\tau}` (auto-covariance) or
    :math:`R^{-}_\\tau = R_{\\text{neg},\\tau}` (cross-covariance).

    Parameters
    ----------
    R : ndarray, shape ``(M+1,)``
        Covariance for lags ``0 .. M``.
    W : ndarray, shape ``(M+1,)``
        Hann window values.
    w : float
        Frequency in rad/sample.
    R_neg : ndarray or None, optional
        Covariance for negative lags (cross-covariance case).

    Returns
    -------
    complex
        Scalar spectral estimate.
    """
    M = len(W) - 1

    # Lag 0 contribution
    val = W[0] * R[0]

    # Lags 1 .. M
    for tau in range(1, M + 1):
        e = np.exp(-1j * w * tau)
        if R_neg is None:
            val += W[tau] * (R[tau] * e + np.conj(R[tau]) * np.conj(e))
        else:
            val += W[tau] * (R[tau] * e + R_neg[tau] * np.conj(e))

    return complex(val)


def _matrix_single_freq_dft(
    R: np.ndarray,
    W: np.ndarray,
    w: float,
    p: int,
    q: int,
    R_neg: np.ndarray | None = None,
) -> np.ndarray:
    """Windowed DFT at a single frequency for matrix-valued covariances.

    Loops over the ``(p, q)`` elements and delegates to
    :func:`_scalar_single_freq_dft`.

    Parameters
    ----------
    R : ndarray
        Covariance array, shape ``(M+1,)`` or ``(M+1, p, q)``.
    W : ndarray, shape ``(M+1,)``
        Hann window values.
    w : float
        Frequency in rad/sample.
    p : int
        Number of rows of the spectral matrix.
    q : int
        Number of columns of the spectral matrix.
    R_neg : ndarray or None, optional
        Reverse covariance for negative lags, shape ``(M+1,)`` or
        ``(M+1, q, p)`` (note transposed indices).

    Returns
    -------
    Phi : ndarray, shape ``(p, q)``, complex
        Spectral matrix at the given frequency.
    """
    Phi = np.zeros((p, q), dtype=complex)

    for ii in range(p):
        for jj in range(q):
            if p == 1 and q == 1:
                R_vec = R.ravel()
                R_neg_vec = None if R_neg is None else R_neg.ravel()
            else:
                R_vec = R[:, ii, jj]
                R_neg_vec = None if R_neg is None else R_neg[:, jj, ii]

            Phi[ii, jj] = _scalar_single_freq_dft(R_vec, W, w, R_neg_vec)

    return Phi
