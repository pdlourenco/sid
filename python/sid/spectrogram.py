# Copyright (c) 2026 Pedro Lourenco. All rights reserved.
# This code is released under the MIT License. See LICENSE file in the
# project root for full license information.
#
# This module is part of the Open Source System Identification Toolbox (SID).
# https://github.com/pdlourenco/sid

"""Short-time FFT spectrogram for time-frequency analysis."""

from __future__ import annotations

import math

import numpy as np

from sid._exceptions import SidError
from sid._results import SpectrogramResult


# ---------------------------------------------------------------------------
# Private helpers
# ---------------------------------------------------------------------------


def _build_window(window: str | np.ndarray, L: int) -> np.ndarray:
    """Create a window vector of length *L*.

    Parameters
    ----------
    window : str or ndarray
        Window specification.  One of ``'hann'``, ``'hamming'``,
        ``'rect'`` / ``'rectangular'``, or a 1-D NumPy array of length *L*.
    L : int
        Window (segment) length.

    Returns
    -------
    ndarray, shape ``(L,)``
        The window vector.

    Raises
    ------
    SidError
        If the window name is unrecognised (code: ``'invalid_window'``) or
        a numeric array has the wrong length (code: ``'window_size_mismatch'``).
    """
    if isinstance(window, np.ndarray):
        w = np.asarray(window, dtype=np.float64).ravel()
        if len(w) != L:
            raise SidError(
                "window_size_mismatch",
                f"Custom window vector length ({len(w)}) must equal window_length ({L}).",
            )
        return w

    if not isinstance(window, str):
        raise SidError(
            "invalid_window",
            "Window must be 'hann', 'hamming', 'rect', or a numeric array.",
        )

    n = np.arange(L, dtype=np.float64)
    name = window.lower()

    if name == "hann":
        return 0.5 * (1.0 - np.cos(2.0 * np.pi * n / (L - 1)))
    if name == "hamming":
        return 0.54 - 0.46 * np.cos(2.0 * np.pi * n / (L - 1))
    if name in ("rect", "rectangular"):
        return np.ones(L, dtype=np.float64)

    raise SidError(
        "invalid_window",
        f"Unknown window type '{window}'. Use 'hann', 'hamming', or 'rect'.",
    )


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def spectrogram(
    x: np.ndarray,
    *,
    window_length: int = 256,
    overlap: int | None = None,
    nfft: int | None = None,
    window: str | np.ndarray = "hann",
    sample_time: float = 1.0,
) -> SpectrogramResult:
    """Compute the short-time FFT spectrogram of one or more signals.

    Divides the signal into overlapping segments, applies a window, and
    computes the one-sided power spectral density via FFT.  This is a
    dependency-free replacement for the Signal Processing Toolbox
    ``spectrogram`` function.

    Parameters
    ----------
    x : ndarray, shape ``(N,)``, ``(N, n_ch)``, or ``(N, n_ch, n_traj)``
        Real-valued signal data.  A 1-D array is treated as a
        single-channel signal.  Each column of a 2-D array is a separate
        channel.  A 3-D array ``(N, n_ch, n_traj)`` provides multiple
        trajectories; the power spectral density is ensemble-averaged
        across trajectories within each segment.
    window_length : int, optional
        Segment length *L*.  Must be a positive integer.  Default is
        ``256``.
    overlap : int or None, optional
        Overlap *P* between consecutive segments, ``0 <= P < L``.
        Default is ``floor(L / 2)``.
    nfft : int or None, optional
        FFT length.  Must be an integer >= *L*.  Default is
        ``max(256, 2**ceil(log2(L)))``.
    window : str or ndarray, optional
        Window function applied to each segment before FFT.  Accepted
        strings are ``'hann'`` (default), ``'hamming'``, ``'rect'``
        (or ``'rectangular'``).  Alternatively, pass a 1-D NumPy array
        of length *L*.
    sample_time : float, optional
        Sample time in seconds.  Must be positive.  Default is ``1.0``.

    Returns
    -------
    SpectrogramResult
        Frozen dataclass with fields:

        - **time** (*ndarray, shape (K,)*) -- Center time of each
          segment in seconds.
        - **frequency** (*ndarray, shape (n_bins,)*) -- Frequency
          vector in Hz.
        - **frequency_rad** (*ndarray, shape (n_bins,)*) -- Frequency
          vector in rad/s.
        - **power** (*ndarray, shape (n_bins, K)* or
          *(n_bins, K, n_ch)*) -- One-sided power spectral density.
        - **power_db** (*ndarray, shape (n_bins, K)* or
          *(n_bins, K, n_ch)*) -- Power in dB,
          ``10 * log10(max(power, eps))``.
        - **complex_stft** (*ndarray, shape (n_bins, K)* or
          *(n_bins, K, n_ch)*) -- Complex STFT coefficients
          (ensemble-averaged across trajectories).
        - **sample_time** (*float*) -- Sample time in seconds.
        - **window_length** (*int*) -- Segment length *L*.
        - **overlap** (*int*) -- Overlap *P*.
        - **nfft** (*int*) -- FFT length.
        - **num_trajectories** (*int*) -- Number of trajectories used.
        - **method** (*str*) -- Always ``'spectrogram'``.

    Raises
    ------
    SidError
        If *x* is complex (code: ``'complex_data'``).
    SidError
        If *x* contains NaN or Inf (code: ``'non_finite'``).
    SidError
        If *window_length* is not a positive integer
        (code: ``'invalid_window_length'``).
    SidError
        If the signal is shorter than *window_length*
        (code: ``'too_short'``).
    SidError
        If *sample_time* is not positive
        (code: ``'invalid_sample_time'``).
    SidError
        If *overlap* is not an integer in ``[0, L-1]``
        (code: ``'invalid_overlap'``).
    SidError
        If *nfft* is not an integer >= *L*
        (code: ``'invalid_nfft'``).
    SidError
        If the data is too short for even one segment
        (code: ``'too_few_segments'``).
    SidError
        If the window name is unrecognised
        (code: ``'invalid_window'``) or a custom window array has
        the wrong length (code: ``'window_size_mismatch'``).

    Examples
    --------
    Spectrogram of a chirp signal:

    >>> import numpy as np
    >>> import sid  # doctest: +SKIP
    >>> Fs = 1000; Ts = 1 / Fs; N = 5000
    >>> t = np.arange(N) * Ts
    >>> x = np.cos(2 * np.pi * (50 + 100 * t / t[-1]) * t)
    >>> result = sid.spectrogram(x, window_length=256, sample_time=Ts)  # doctest: +SKIP
    >>> result.power.shape  # doctest: +SKIP
    (129, 18)

    Multi-channel input:

    >>> x2 = np.column_stack([x, 0.5 * x])  # doctest: +SKIP
    >>> result = sid.spectrogram(x2, window_length=256, sample_time=Ts)  # doctest: +SKIP
    >>> result.power.shape  # doctest: +SKIP
    (129, 18, 2)

    Multi-trajectory ensemble averaging:

    >>> rng = np.random.default_rng(0)  # doctest: +SKIP
    >>> x3d = rng.standard_normal((1000, 1, 5))  # doctest: +SKIP
    >>> result = sid.spectrogram(x3d, window_length=128)  # doctest: +SKIP
    >>> result.num_trajectories  # doctest: +SKIP
    5

    Notes
    -----
    **Algorithm (SPEC.md S7):**

    1. Divide the signal into ``K`` overlapping segments of length *L*
       with step ``L - P``.
    2. Apply the time-domain window to each segment.
    3. Compute the FFT of each windowed segment (zero-padded to *nfft*).
    4. Compute the one-sided power spectral density per segment:
       ``P_k(m) = (1 / (Fs * S1)) * |X_k(m)|^2``, where
       ``S1 = sum(w**2)`` is the window power.  Positive-frequency bins
       (excluding DC and Nyquist for even *nfft*) are doubled.
    5. For multi-trajectory data, PSD and complex STFT coefficients are
       ensemble-averaged across trajectories.

    The time-domain window reduces spectral leakage; it is distinct from
    the lag-domain Hann window used in :func:`sid.freq_bt`.

    **Specification:** SPEC.md S7 -- Short-Time Spectral Analysis

    References
    ----------
    .. [1] Oppenheim, A.V. and Schafer, R.W., "Discrete-Time Signal
       Processing", 3rd ed., Prentice Hall, 2010.

    See Also
    --------
    sid.freq_bt : Blackman-Tukey spectral analysis.
    sid.freq_etfe : Empirical transfer function estimate.
    sid.spectrogram_plot : Plot a spectrogram result.

    Changelog
    ---------
    2026-04-08 : First version (Python port) by Pedro Lourenco.
    """

    # ---- Coerce to ndarray ----
    x = np.asarray(x)

    # ---- Validate signal (complex check before dtype cast) ----
    if np.iscomplexobj(x):
        raise SidError("complex_data", "Complex data is not supported. Signal x must be real.")

    x = x.astype(np.float64, copy=False)

    # ---- Reshape 1-D to column ----
    if x.ndim == 1:
        x = x[:, np.newaxis]

    N = x.shape[0]
    n_ch = x.shape[1]
    n_traj = x.shape[2] if x.ndim == 3 else 1

    L = window_length
    Ts = sample_time

    if not np.all(np.isfinite(x)):
        raise SidError("non_finite", "Signal x contains NaN or Inf values.")

    if not isinstance(L, (int, np.integer)) or L < 1:
        raise SidError(
            "invalid_window_length",
            "WindowLength must be a positive integer.",
        )

    if N < L:
        raise SidError(
            "too_short",
            f"Signal length ({N}) is shorter than window_length ({L}).",
        )

    if Ts <= 0:
        raise SidError(
            "invalid_sample_time",
            "SampleTime must be a positive scalar.",
        )

    # ---- Defaults that depend on L ----
    P = overlap if overlap is not None else L // 2

    if nfft is not None:
        nfft_val = nfft
    else:
        nfft_val = max(256, 2 ** math.ceil(math.log2(L)))

    # ---- Validate derived parameters ----
    if not isinstance(P, (int, np.integer)) or P < 0 or P >= L:
        raise SidError(
            "invalid_overlap",
            "Overlap must be an integer in [0, L-1].",
        )

    if not isinstance(nfft_val, (int, np.integer)) or nfft_val < L:
        raise SidError(
            "invalid_nfft",
            "NFFT must be an integer >= window_length.",
        )

    # ---- Build window vector ----
    w = _build_window(window, L)

    # ---- Segmentation (SPEC.md S7.2) ----
    step = L - P
    K = (N - L) // step + 1

    if K < 1:
        raise SidError(
            "too_few_segments",
            f"Data too short for even one segment with L={L}, P={P}.",
        )

    Fs = 1.0 / Ts
    n_bins = nfft_val // 2 + 1
    S1 = np.sum(w**2)

    # ---- Pre-allocate ----
    stft_coeffs = np.zeros((n_bins, K, n_ch), dtype=np.complex128)
    Pxx = np.zeros((n_bins, K, n_ch), dtype=np.float64)

    # ---- Compute STFT (SPEC.md S7.2-7.3) ----
    for ch in range(n_ch):
        for k in range(K):
            start_idx = k * step

            Pk_sum = np.zeros(n_bins, dtype=np.float64)
            X_sum = np.zeros(n_bins, dtype=np.complex128)

            for lt in range(n_traj):
                if n_traj > 1:
                    seg = x[start_idx : start_idx + L, ch, lt] * w
                else:
                    seg = x[start_idx : start_idx + L, ch] * w

                X = np.fft.fft(seg, nfft_val)[:n_bins]
                X_sum += X

                # One-sided PSD: P(w) = (1 / (Fs * S1)) * |X(w)|^2
                Pk = (1.0 / (Fs * S1)) * np.abs(X) ** 2

                # Double positive-frequency bins for one-sided spectrum
                if nfft_val % 2 == 0:
                    Pk[1:-1] *= 2.0  # all except DC and Nyquist
                else:
                    Pk[1:] *= 2.0  # all except DC (no Nyquist for odd nfft)

                Pk_sum += Pk

            # Ensemble-average across trajectories
            stft_coeffs[:, k, ch] = X_sum / n_traj
            Pxx[:, k, ch] = Pk_sum / n_traj

    # ---- Time vector (center of each segment) ----
    time_vec = (np.arange(K, dtype=np.float64) * step + L / 2.0) * Ts

    # ---- Frequency vectors ----
    freq_hz = np.arange(n_bins, dtype=np.float64) * Fs / nfft_val
    freq_rad = 2.0 * np.pi * freq_hz

    # ---- Power in dB ----
    Pxx_db = 10.0 * np.log10(np.maximum(Pxx, np.finfo(np.float64).eps))

    # ---- Squeeze single-channel outputs to 2-D ----
    if n_ch == 1:
        Pxx = Pxx[:, :, 0]
        Pxx_db = Pxx_db[:, :, 0]
        stft_coeffs = stft_coeffs[:, :, 0]

    # ---- Pack result ----
    return SpectrogramResult(
        time=time_vec,
        frequency=freq_hz,
        frequency_rad=freq_rad,
        power=Pxx,
        power_db=Pxx_db,
        complex_stft=stft_coeffs,
        sample_time=Ts,
        window_length=L,
        overlap=P,
        nfft=nfft_val,
        num_trajectories=n_traj,
        method="spectrogram",
    )
