# Copyright (c) 2026 Pedro Lourenco. All rights reserved.
# This code is released under the MIT License. See LICENSE file in the
# project root for full license information.
#
# This module is part of the Open Source System Identification Toolbox (SID).
# https://github.com/pdlourenco/sid-matlab

"""Discrete Fourier transform at specified frequencies."""

from __future__ import annotations

import math

import numpy as np


def sid_dft(x: np.ndarray, freqs: np.ndarray, use_fft: bool) -> np.ndarray:
    """Compute the DFT of a time-domain signal at specified frequencies.

    This is the Python port of ``sidDFT.m``.

    Computes:

    .. math::

        X(\\omega_k) = \\sum_{t=1}^{N} x(t)\\,e^{-j\\omega_k t}

    Parameters
    ----------
    x : ndarray, shape ``(N,)`` or ``(N, p)``
        Real data matrix (*p* channels).
    freqs : ndarray, shape ``(nf,)``
        Frequency vector in rad/sample, in ``(0, pi]``.
    use_fft : bool
        If ``True``, use the FFT fast path (requires the default linear
        grid of 128 points).

    Returns
    -------
    X : ndarray, shape ``(nf, p)``, complex
        DFT values at the requested frequencies.

    Examples
    --------
    >>> X = sid_dft(x, freqs, use_fft=True)  # doctest: +SKIP

    See Also
    --------
    sid._internal.windowed_dft.windowed_dft : Windowed DFT of covariances.

    Changelog
    ---------
    2026-04-08 : First version by Pedro Lourenco.
    """

    # Handle 1-D input
    squeeze_output = False
    if x.ndim == 1:
        x = x[:, np.newaxis]
        squeeze_output = True

    N: int = x.shape[0]
    p: int = x.shape[1]
    nf: int = len(freqs)

    if use_fft:
        # FFT fast path
        L = 2 * nf
        Nfft = max(N, L)
        # Round up to next power of 2
        Nfft = 2 ** math.ceil(math.log2(Nfft))

        Xfft = np.fft.fft(x, n=Nfft, axis=0)  # (Nfft, p)

        # Compute bin indices (0-based)
        bin_idx = np.round(freqs * Nfft / (2.0 * np.pi)).astype(int)
        # Clamp to valid range [1, Nfft//2]
        bin_idx = np.clip(bin_idx, 1, Nfft // 2)

        # Extract the desired bins (0-based indexing)
        X = Xfft[bin_idx, :]  # (nf, p)

        # Phase correction: MATLAB fft gives sum x(t)*exp(-jw*(t-1)),
        # we want sum x(t)*exp(-jw*t) = exp(-jw) * fft result
        w_actual = bin_idx * 2.0 * np.pi / Nfft
        correction = np.exp(-1j * w_actual)[:, np.newaxis]  # (nf, 1)
        X = X * correction
    else:
        # Direct DFT at arbitrary frequencies
        X = np.zeros((nf, p), dtype=np.complex128)
        t = np.arange(1, N + 1, dtype=np.float64)  # 1-based time

        for k in range(nf):
            e = np.exp(-1j * freqs[k] * t)  # (N,)
            X[k, :] = e @ x  # (p,)

    # Keep 2-D output even for scalar input — matches MATLAB convention
    if squeeze_output:
        pass  # still return (nf, 1) shape — caller can squeeze if needed

    return X
