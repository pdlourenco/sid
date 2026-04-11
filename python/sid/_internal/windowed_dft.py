# Copyright (c) 2026 Pedro Lourenco. All rights reserved.
# This code is released under the MIT License. See LICENSE file in the
# project root for full license information.
#
# This module is part of the Open Source System Identification Toolbox (SID).
# https://github.com/pdlourenco/sid

"""Windowed Fourier transform of covariance estimates."""

from __future__ import annotations

import numpy as np


def windowed_dft(
    R: np.ndarray,
    W: np.ndarray,
    freqs: np.ndarray,
    use_fft: bool,
    R_neg: np.ndarray | None = None,
) -> np.ndarray:
    """Windowed Fourier transform of covariance estimates.

    This is the Python port of ``sidWindowedDFT.m``.

    Computes the spectral estimate:

    .. math::

        \\Phi(\\omega) = \\sum_{\\tau=-M}^{M}
            R(\\tau)\\,W(|\\tau|)\\,e^{-j\\omega\\tau}

    For auto-covariance, :math:`R(-\\tau) = \\overline{R(\\tau)}` is used
    automatically.  For cross-covariance, pass *R_neg* so that
    :math:`R(-\\tau) = R_{\\text{neg}}(\\tau)`.

    Parameters
    ----------
    R : ndarray
        Covariance estimates for lags ``0 .. M``.  Shape ``(M+1,)`` for
        scalar signals or ``(M+1, p, q)`` for matrix signals.
    W : ndarray, shape ``(M+1,)``
        Hann window values for lags ``0 .. M``.
    freqs : ndarray, shape ``(nf,)``
        Frequency vector in rad/sample.
    use_fft : bool
        If ``True``, use the FFT fast path (requires the default 128-point
        linear grid).
    R_neg : ndarray or None, optional
        Covariance for negative lags.  Same size as *R*.  When provided,
        ``R_neg[tau]`` is used for :math:`R(-\\tau)` instead of
        ``conj(R[tau])``.  Required for cross-covariance estimation.
        For matrix signals, note that negative-lag indices are transposed:
        ``R_neg[:, j, i]`` is used for the ``(i, j)`` element.

    Returns
    -------
    Phi : ndarray
        Spectral estimate at each frequency.  Shape ``(nf,)`` for scalar
        signals or ``(nf, p, q)`` for matrix signals.

    Examples
    --------
    >>> Phi = windowed_dft(R, W, freqs, use_fft=True)  # doctest: +SKIP

    Notes
    -----
    **Specification:** SPEC.md §2.5 -- Windowed Spectral Estimates

    See Also
    --------
    sid._internal.cov.sid_cov : Computes the covariances fed to this function.
    sid._internal.hann_win.hann_win : Computes the window fed to this function.

    Changelog
    ---------
    2026-04-08 : First version by Pedro Lourenco.
    """

    nf: int = len(freqs)

    # Determine signal dimensions
    if R.ndim == 1:
        # Scalar signals
        if use_fft:
            return _fft_path(R, W, nf, R_neg)
        return _direct_path(R, W, freqs, R_neg)

    # Matrix signals: R has shape (M+1, p, q)
    p: int = R.shape[1]
    q: int = R.shape[2]
    Phi = np.zeros((nf, p, q), dtype=np.complex128)

    for ii in range(p):
        for jj in range(q):
            R_ij = R[:, ii, jj]
            if R_neg is None:
                R_neg_ij = None
            else:
                # Transposed indices for negative lags: R_neg[:, j, i]
                R_neg_ij = R_neg[:, jj, ii]

            if use_fft:
                Phi[:, ii, jj] = _fft_path(R_ij, W, nf, R_neg_ij)
            else:
                Phi[:, ii, jj] = _direct_path(R_ij, W, freqs, R_neg_ij)

    return Phi


def _fft_path(
    R: np.ndarray,
    W: np.ndarray,
    nf: int,
    R_neg: np.ndarray | None,
) -> np.ndarray:
    """FFT fast path for the default linear frequency grid.

    Outputs the spectral estimate at the default-grid frequencies
    ``omega_k = k * pi / nf`` for ``k = 1 .. nf``.  The FFT length ``L``
    is chosen per SPEC.md S2.5.1 to satisfy two constraints simultaneously:

    1. ``L >= 2*M+1`` so that positive lags ``s[0..M]`` and the wrapped
       negative lags ``s[L-M..L-1]`` do not collide in the buffer.
    2. ``L`` is an integer multiple of ``2*nf`` so that the FFT bins
       remain aligned with the default frequency grid.

    The smallest ``L`` satisfying both is ``L = 2*nf*K`` where
    ``K = ceil((2*M + 1) / (2*nf))``.  For the typical regime ``M << nf``
    (e.g. the default ``M = min(N // 10, 30)`` with ``nf = 128``),
    ``K = 1`` and ``L = 256``, identical to the historical fast path.
    For larger ``M`` the FFT length grows just enough to fit the lag
    sequence, and the desired bins are recovered by striding the FFT
    output by ``K``.

    Parameters
    ----------
    R : ndarray, shape ``(M+1,)``
        Covariance for lags ``0 .. M``.
    W : ndarray, shape ``(M+1,)``
        Window values.
    nf : int
        Number of frequency points.
    R_neg : ndarray or None
        Covariance for negative lags (scalar sequence).

    Returns
    -------
    Phi : ndarray, shape ``(nf,)``, complex
    """

    M: int = len(W) - 1

    # K = ceil((2M+1) / (2*nf)) using integer arithmetic; L = 2*nf*K.
    # See SPEC.md S2.5.1 step 2.
    two_nf: int = 2 * nf
    K: int = max(1, (2 * M + 1 + two_nf - 1) // two_nf)
    L: int = two_nf * K

    s = np.zeros(L, dtype=np.complex128)

    # Lag 0
    s[0] = W[0] * R[0]

    # Positive lags 1 .. M  (safe: L >= 2*M+1 implies M < L)
    for tau in range(1, M + 1):
        s[tau] = W[tau] * R[tau]

    # Negative lags -1 .. -M  (wrapped into positions L-1 .. L-M;
    # safe and disjoint from positive lags because L - M >= M + 1)
    for tau in range(1, M + 1):
        if R_neg is None:
            s[L - tau] = W[tau] * np.conj(R[tau])
        else:
            s[L - tau] = W[tau] * R_neg[tau]

    S = np.fft.fft(s, n=L)

    # Extract bins K, 2K, ..., nf*K  (SPEC.md S2.5.1 step 4).
    # Bin k*K of a length-L FFT sits at frequency k*K * 2*pi/L
    # = k*K * 2*pi/(2*nf*K) = k * pi/nf = omega_k, as required.
    return S[K : K * (nf + 1) : K]


def _direct_path(
    R: np.ndarray,
    W: np.ndarray,
    freqs: np.ndarray,
    R_neg: np.ndarray | None,
) -> np.ndarray:
    """Direct DFT at arbitrary frequencies.

    Parameters
    ----------
    R : ndarray, shape ``(M+1,)``
        Covariance for lags ``0 .. M``.
    W : ndarray, shape ``(M+1,)``
        Window values.
    freqs : ndarray, shape ``(nf,)``
        Frequency vector in rad/sample.
    R_neg : ndarray or None
        Covariance for negative lags (scalar sequence).

    Returns
    -------
    Phi : ndarray, shape ``(nf,)``, complex
    """

    M: int = len(W) - 1
    nf: int = len(freqs)
    Phi = np.zeros(nf, dtype=np.complex128)

    for k in range(nf):
        w = freqs[k]

        # Lag 0 contribution
        val = W[0] * R[0]

        # Lags 1 .. M: positive and negative combined
        for tau in range(1, M + 1):
            e = np.exp(-1j * w * tau)
            if R_neg is None:
                val += W[tau] * (R[tau] * e + np.conj(R[tau]) * np.conj(e))
            else:
                val += W[tau] * (R[tau] * e + R_neg[tau] * np.conj(e))

        Phi[k] = val

    return Phi
