# Copyright (c) 2026 Pedro Lourenco. All rights reserved.
# This code is released under the MIT License. See LICENSE file in the
# project root for full license information.
#
# This module is part of the Open Source System Identification Toolbox (SID).
# https://github.com/pdlourenco/sid

"""Unit tests for the frequency-domain simulation helper (SPEC §15).

Port-symmetric with test_sidFreqDomainSim.m (#124). The helper is also
cross-validated against MATLAB via reference_freq_domain_sim.json; these are
independent-oracle unit tests.
"""

from __future__ import annotations

import numpy as np

from sid._internal.freq_domain_sim import freq_domain_sim


class TestFreqDomainSim:
    """Frequency-domain simulation via IFFT."""

    def test_constant_gain_scales_input(self) -> None:
        """A constant gain G(w) = k scales a zero-mean input by k (#124).

        With a dense grid over (0, pi] and zero-mean u, the DC/out-of-grid
        zeroing is negligible, so the round-trip reproduces k*u to machine
        precision — a clean independent oracle for the multiply-in-frequency.
        """
        rng = np.random.default_rng(124)
        N = 256
        u = rng.standard_normal((N, 1))
        u -= u.mean()
        freqs = np.arange(1, N // 2 + 1) * np.pi / (N // 2)
        G = 2.0 * np.ones(len(freqs), dtype=complex)

        Y = np.asarray(freq_domain_sim(G, freqs, u, N)).ravel()

        np.testing.assert_allclose(
            Y,
            2.0 * u.ravel(),
            atol=1e-10,
            err_msg="constant gain should scale the input by 2",
        )

    def test_out_of_grid_frequencies_zeroed(self) -> None:
        """Frequencies outside the model grid contribute nothing (#124)."""
        rng = np.random.default_rng(125)
        N = 128
        u = rng.standard_normal((N, 1))
        # Model grid covers only the low third of (0, pi]; high-freq content of
        # u must be dropped, so the output energy is strictly less than input.
        freqs = np.linspace(0.05, np.pi / 3, 30)
        G = np.ones(len(freqs), dtype=complex)
        Y = np.asarray(freq_domain_sim(G, freqs, u, N)).ravel()
        assert np.sum(Y**2) < np.sum(u.ravel() ** 2), (
            "out-of-grid (high-freq) content should be removed"
        )
