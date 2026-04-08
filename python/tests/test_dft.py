# Copyright (c) 2026 Pedro Lourenco. All rights reserved.
# This code is released under the MIT License. See LICENSE file in the
# project root for full license information.
#
# This module is part of the Open Source System Identification Toolbox (SID).
# https://github.com/pdlourenco/sid-matlab

"""Tests for sid_dft from sid._internal.dft."""

from __future__ import annotations

import numpy as np

from sid._internal.dft import sid_dft


class TestSidDFT:
    """Unit tests for the DFT computation."""

    def test_impulse(self) -> None:
        """DFT of impulse x=[1,0,0,0]: X(w) = exp(-jw)."""
        x = np.array([1.0, 0.0, 0.0, 0.0])
        freqs = np.array([np.pi / 4, np.pi / 2, np.pi])
        X = sid_dft(x, freqs, use_fft=False)
        expected = np.exp(-1j * freqs)[:, np.newaxis]
        np.testing.assert_allclose(X, expected, atol=1e-12)

    def test_constant_signal(self) -> None:
        """DFT of constant signal at w=pi/3 matches direct sum."""
        N = 16
        x = np.ones(N)
        w = np.array([np.pi / 3])
        X = sid_dft(x, w, use_fft=False)
        t = np.arange(1, N + 1, dtype=np.float64)
        expected = np.sum(np.exp(-1j * w[0] * t))
        np.testing.assert_allclose(X[0, 0], expected, atol=1e-10)

    def test_fft_vs_direct(self) -> None:
        """FFT and direct DFT agree on default 128-point grid."""
        rng = np.random.default_rng(42)
        N = 200
        x = rng.standard_normal(N)
        freqs = np.arange(1, 129) * np.pi / 128
        X_fft = sid_dft(x, freqs, use_fft=True)
        X_direct = sid_dft(x, freqs, use_fft=False)
        rel_err = np.max(np.abs(X_fft - X_direct)) / np.max(np.abs(X_direct))
        assert rel_err < 0.05

    def test_multi_channel(self) -> None:
        """Multi-channel (50x3) input, 3 freqs gives (3,3) output."""
        rng = np.random.default_rng(7)
        N = 50
        x = rng.standard_normal((N, 3))
        freqs = np.array([0.5, 1.0, 2.0])
        X = sid_dft(x, freqs, use_fft=False)
        assert X.shape == (3, 3)
        # Verify first channel manually
        t = np.arange(1, N + 1, dtype=np.float64)
        X1 = np.zeros(3, dtype=np.complex128)
        for k in range(3):
            X1[k] = np.sum(x[:, 0] * np.exp(-1j * freqs[k] * t))
        np.testing.assert_allclose(X[:, 0], X1, atol=1e-10)

    def test_output_dims(self) -> None:
        """(100x2) input, 64 freqs gives (64, 2) output."""
        rng = np.random.default_rng(1)
        x = rng.standard_normal((100, 2))
        freqs = np.arange(1, 65) * np.pi / 64
        X = sid_dft(x, freqs, use_fft=False)
        assert X.shape == (64, 2)

    def test_single_frequency(self) -> None:
        """N=50, single freq pi/2 gives (1,1) shaped output."""
        rng = np.random.default_rng(42)
        x = rng.standard_normal(50)
        freqs = np.array([np.pi / 2])
        X = sid_dft(x, freqs, use_fft=False)
        assert X.shape == (1, 1)
        t = np.arange(1, 51, dtype=np.float64)
        expected = np.sum(x * np.exp(-1j * np.pi / 2 * t))
        np.testing.assert_allclose(X[0, 0], expected, atol=1e-10)

    def test_parseval_energy(self) -> None:
        """Parseval theorem: (1/N) sum |X(w_k)|^2 = sum |x(t)|^2."""
        rng = np.random.default_rng(88)
        N = 64
        x = rng.standard_normal(N)
        # Fourier frequencies: w_k = 2*pi*k/N for k=1..N-1 (exclude DC)
        freqs_full = np.arange(1, N) * 2 * np.pi / N
        X = sid_dft(x, freqs_full, use_fft=False)
        # DC component manually: X(0) = sum(x)
        X_dc = np.sum(x)
        energy_freq = (np.abs(X_dc) ** 2 + np.sum(np.abs(X) ** 2)) / N
        energy_time = np.sum(x**2)
        np.testing.assert_allclose(energy_freq, energy_time, atol=1e-8)
