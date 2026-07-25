# Copyright (c) 2026 Pedro Lourenco. All rights reserved.
# This code is released under the MIT License. See LICENSE file in the
# project root for full license information.
#
# This module is part of the Open Source System Identification Toolbox (SID).
# https://github.com/pdlourenco/sid

"""Unit test for the COSMIC block-tridiagonal solver's ill-conditioning
diagnostic (SPEC §8.3.4, #120).

Port-symmetric with test_sidLTVcosmicSolve.m. The §8.3.4 warning fires when a
forward-pass information block ``Lbd(k-1)`` is near-singular; it is a defensive
diagnostic not reachable from normal ``ltv_disc`` inputs, so it is exercised
here by feeding the internal solver a deliberately near-singular first block.
"""

from __future__ import annotations

import warnings

import numpy as np

from sid._internal.ltv_cosmic_solve import cosmic_solve


class TestCosmicSolveDiagnostic:
    """SPEC §8.3.4 ill-conditioning warning."""

    def test_near_singular_block_warns(self) -> None:
        """A near-singular forward-pass block warns while still returning a
        finite result (#120)."""
        p, q, n = 1, 1, 2
        d = p + q
        S = np.zeros((d, d, n))
        S[:, :, 0] = np.diag([1.0, 1e-17])  # rcond ~1e-17 < eps -> warns
        S[:, :, 1] = 5.0 * np.eye(d)
        T = np.ones((d, p, n))
        lam = np.ones(n - 1)

        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            C, _ = cosmic_solve(S, T, lam, n, p, q)

        assert np.all(np.isfinite(C)), "solver must still return a finite result"
        assert any(
            "near-singular" in str(x.message).lower()
            or "rcond" in str(x.message).lower()
            for x in caught
        ), "expected §8.3.4 ill-conditioning warning"
