# Copyright (c) 2026 Pedro Lourenco. All rights reserved.
# This code is released under the MIT License. See LICENSE file in the
# project root for full license information.
#
# This module is part of the Open Source System Identification Toolbox (SID).
# https://github.com/pdlourenco/sid

"""Frequency-domain model simulation via FFT/IFFT."""

from __future__ import annotations

import numpy as np


def freq_domain_sim(
    G_model: np.ndarray,
    freqs_model: np.ndarray,
    u: np.ndarray,
    N: int,
) -> np.ndarray:
    """Simulate frequency-domain model output via IFFT.

    This is the Python port of ``sidFreqDomainSim.m``.

    Filters the input signal *u* through the frequency response *G_model*
    by computing the FFT of *u*, interpolating *G_model* onto the FFT
    frequency grid (using log-magnitude and unwrapped phase), multiplying,
    and taking the IFFT.

    Frequencies outside the model grid are set to zero (no extrapolation).

    Parameters
    ----------
    G_model : ndarray, shape ``(nf, ny, nu)``
        Complex frequency response.  If 1-D or 2-D, it is broadcast to
        3-D internally.
    freqs_model : ndarray, shape ``(nf,)``
        Frequency vector in rad/sample, in ``(0, pi]``.
    u : ndarray, shape ``(N, nu)``
        Real input signal.
    N : int
        Number of samples (must equal ``u.shape[0]``).

    Returns
    -------
    Y_pred : ndarray, shape ``(N, ny)``
        Predicted output signal.

    Examples
    --------
    >>> import numpy as np
    >>> from sid._internal.freq_domain_sim import freq_domain_sim
    >>> G = np.ones((64, 1, 1), dtype=complex)
    >>> freqs = np.linspace(0.01, np.pi, 64)
    >>> u = np.random.randn(128, 1)
    >>> y = freq_domain_sim(G, freqs, u, 128)  # doctest: +SKIP

    Notes
    -----
    **Specification:** (Frequency-domain simulation helper -- not a standalone SPEC.md section)

    See Also
    --------
    sid.compare : Model output comparison using this helper.
    sid.residual : Residual analysis using this helper.

    Changelog
    ---------
    2026-04-09 : First version (Python port) by Pedro Lourenco.
    """
    # ------------------------------------------------------------------
    # 1. Ensure G_model is 3-D: (nf, ny, nu)
    # ------------------------------------------------------------------
    G_model = np.asarray(G_model)
    if G_model.ndim == 1:
        G_model = G_model[:, np.newaxis, np.newaxis]
    elif G_model.ndim == 2:
        G_model = G_model[:, :, np.newaxis]

    ny: int = G_model.shape[1]
    nu: int = G_model.shape[2]

    # ------------------------------------------------------------------
    # 2. Build FFT frequency grid
    # ------------------------------------------------------------------
    nfft: int = N
    npos: int = nfft // 2
    freqs_fft = np.arange(1, npos + 1) * (2.0 * np.pi / nfft)

    # ------------------------------------------------------------------
    # 3. FFT of input
    # ------------------------------------------------------------------
    U_fft = np.fft.fft(u, n=nfft, axis=0)

    # ------------------------------------------------------------------
    # 4. Interpolate G onto FFT grid and multiply
    # ------------------------------------------------------------------
    Y_pred_fft = np.zeros((nfft, ny), dtype=complex)

    for iy in range(ny):
        for iu in range(nu):
            Gij = G_model[:, iy, iu].ravel()
            G_interp = _interp_g_safe(freqs_model.ravel(), Gij, freqs_fft)
            # Bins 1..npos correspond to positive frequencies
            Y_pred_fft[1 : npos + 1, iy] += G_interp * U_fft[1 : npos + 1, iu]

    # ------------------------------------------------------------------
    # 5. Conjugate symmetry for real output
    # ------------------------------------------------------------------
    for iy in range(ny):
        if nfft % 2 == 0:
            Y_pred_fft[npos + 1 :, iy] = np.conj(Y_pred_fft[npos - 1 : 0 : -1, iy])
        else:
            Y_pred_fft[npos + 1 :, iy] = np.conj(Y_pred_fft[npos:0:-1, iy])

    # ------------------------------------------------------------------
    # 6. IFFT to time domain
    # ------------------------------------------------------------------
    Y_pred: np.ndarray = np.real(np.fft.ifft(Y_pred_fft, n=nfft, axis=0))

    return Y_pred


def _interp_g_safe(
    freqs_model: np.ndarray,
    G: np.ndarray,
    freqs_target: np.ndarray,
) -> np.ndarray:
    """Interpolate a complex transfer function with safe boundaries.

    Uses log-magnitude and unwrapped-phase interpolation, which is more
    numerically stable near resonances than real/imaginary interpolation.
    Frequencies outside the model grid are set to zero (no extrapolation).

    Parameters
    ----------
    freqs_model : ndarray, shape ``(nf,)``
        Model frequency vector in rad/sample.
    G : ndarray, shape ``(nf,)``
        Complex transfer function values on *freqs_model*.
    freqs_target : ndarray, shape ``(nt,)``
        Target frequency vector in rad/sample.

    Returns
    -------
    G_interp : ndarray, shape ``(nt,)``, complex
        Interpolated transfer function.  Zero outside the model range.
    """
    mag = np.abs(G)
    ph = np.unwrap(np.angle(G))

    # Avoid log(0): clamp magnitude floor
    mag = np.maximum(mag, np.finfo(float).eps)

    in_range = (freqs_target >= freqs_model[0]) & (freqs_target <= freqs_model[-1])

    G_interp = np.zeros(len(freqs_target), dtype=complex)

    if np.any(in_range):
        logmag_interp = np.interp(freqs_target[in_range], freqs_model, np.log(mag))
        ph_interp = np.interp(freqs_target[in_range], freqs_model, ph)
        G_interp[in_range] = np.exp(logmag_interp) * np.exp(1j * ph_interp)

    return G_interp
