# Copyright (c) 2026 Pedro Lourenco. All rights reserved.
# This code is released under the MIT License. See LICENSE file in the
# project root for full license information.
#
# This module is part of the Open Source System Identification Toolbox (SID).
# https://github.com/pdlourenco/sid-matlab

"""Tests for windowed_dft from sid._internal.windowed_dft."""

from __future__ import annotations

import numpy as np

from sid._internal.cov import sid_cov
from sid._internal.hann_win import hann_win
from sid._internal.windowed_dft import windowed_dft


class TestWindowedDFT:
    """Unit tests for windowed Fourier transform of covariances."""

    def test_fft_vs_direct_scalar(self) -> None:
        """FFT and direct paths agree for scalar auto-covariance."""
        rng = np.random.default_rng(42)
        N = 500
        x = rng.standard_normal((N, 1))
        M = 30
        R = sid_cov(x, x, M)
        W = hann_win(M)
        freqs = np.arange(1, 129) * np.pi / 128

        Phi_fft = windowed_dft(R, W, freqs, use_fft=True)
        Phi_direct = windowed_dft(R, W, freqs, use_fft=False)

        rel_err = np.max(np.abs(Phi_fft - Phi_direct)) / np.max(np.abs(Phi_direct))
        np.testing.assert_allclose(rel_err, 0.0, atol=1e-10)

    def test_auto_spectrum_real(self) -> None:
        """Auto-covariance spectrum has negligible imaginary part."""
        rng = np.random.default_rng(42)
        N = 500
        x = rng.standard_normal((N, 1))
        M = 30
        R = sid_cov(x, x, M)
        W = hann_win(M)
        freqs = np.arange(1, 129) * np.pi / 128

        Phi = windowed_dft(R, W, freqs, use_fft=True)
        assert np.max(np.abs(np.imag(Phi))) < 1e-10

    def test_auto_spectrum_nonneg(self) -> None:
        """Auto-spectrum is non-negative."""
        rng = np.random.default_rng(42)
        N = 500
        x = rng.standard_normal((N, 1))
        M = 30
        R = sid_cov(x, x, M)
        W = hann_win(M)
        freqs = np.arange(1, 129) * np.pi / 128

        Phi = windowed_dft(R, W, freqs, use_fft=True)
        assert np.all(np.real(Phi) > -1e-10)

    def test_output_shape_scalar(self) -> None:
        """Scalar signal output has shape (128,)."""
        rng = np.random.default_rng(42)
        N = 500
        x = rng.standard_normal((N, 1))
        M = 30
        R = sid_cov(x, x, M)
        W = hann_win(M)
        freqs = np.arange(1, 129) * np.pi / 128

        Phi = windowed_dft(R, W, freqs, use_fft=True)
        assert Phi.shape == (128,)

    def test_output_shape_matrix(self) -> None:
        """Matrix signal input (M+1, 2, 3) gives output (nf, 2, 3)."""
        rng = np.random.default_rng(7)
        N = 200
        x = rng.standard_normal((N, 2))
        z = rng.standard_normal((N, 3))
        M = 20
        R = sid_cov(x, z, M)  # (M+1, 2, 3)
        W = hann_win(M)
        freqs = np.arange(1, 65) * np.pi / 64

        Phi = windowed_dft(R, W, freqs, use_fft=False)
        assert Phi.shape == (64, 2, 3)

    def test_fft_vs_direct_matrix(self) -> None:
        """FFT and direct paths agree for matrix auto-covariance."""
        rng = np.random.default_rng(7)
        N = 200
        x = rng.standard_normal((N, 2))
        M = 20
        R = sid_cov(x, x, M)  # (M+1, 2, 2)
        W = hann_win(M)
        freqs = np.arange(1, 129) * np.pi / 128

        Phi_fft = windowed_dft(R, W, freqs, use_fft=True)
        Phi_direct = windowed_dft(R, W, freqs, use_fft=False)

        rel_err = np.max(np.abs(Phi_fft.ravel() - Phi_direct.ravel())) / np.max(
            np.abs(Phi_direct.ravel())
        )
        np.testing.assert_allclose(rel_err, 0.0, atol=1e-10)

    def test_white_noise_spectrum(self) -> None:
        """White noise spectrum is approximately 1 (flat)."""
        rng = np.random.default_rng(99)
        N = 10000
        x = rng.standard_normal((N, 1))
        M = 50
        R = sid_cov(x, x, M)
        W = hann_win(M)
        freqs = np.arange(1, 129) * np.pi / 128

        Phi = np.real(windowed_dft(R, W, freqs, use_fft=True))
        assert np.max(np.abs(Phi - 1)) < 0.25

    def test_custom_frequencies(self) -> None:
        """Non-default frequencies [0.1, 0.5, 1, 2, 3] give shape (5,)."""
        rng = np.random.default_rng(42)
        N = 200
        x = rng.standard_normal((N, 1))
        M = 20
        R = sid_cov(x, x, M)
        W = hann_win(M)
        freqs_custom = np.array([0.1, 0.5, 1.0, 2.0, 3.0])

        Phi = windowed_dft(R, W, freqs_custom, use_fft=False)
        assert Phi.shape == (5,)
        assert np.all(np.real(Phi) > -1e-10)
