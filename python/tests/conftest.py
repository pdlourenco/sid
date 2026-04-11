# Copyright (c) 2026 Pedro Lourenco. All rights reserved.
# This code is released under the MIT License. See LICENSE file in the
# project root for full license information.
#
# This module is part of the Open Source System Identification Toolbox (SID).
# https://github.com/pdlourenco/sid

"""Shared fixtures and helpers for sid test suite."""

from __future__ import annotations

import json
import pathlib
import sys

import numpy as np
import pytest

TESTDATA = pathlib.Path(__file__).resolve().parent.parent.parent / "testdata"

# ---- Examples path injection --------------------------------------------
# Several tests (cross-validation, LTV state estimation, the helper's own
# regression suite) call the spring-mass-damper plant simulator defined in
# python/examples/util_msd.py. That module is a sibling of the test tree,
# not part of the installed ``sid`` package, so we prepend its directory
# to sys.path here so test modules can do ``from util_msd import util_msd``
# without manipulating ``sys.path`` themselves. See
# ``spec/EXAMPLES.md`` §2.1 for the binding contract of the helper.
_EXAMPLES_DIR = pathlib.Path(__file__).resolve().parent.parent / "examples"
if str(_EXAMPLES_DIR) not in sys.path:
    sys.path.insert(0, str(_EXAMPLES_DIR))


@pytest.fixture
def rng():
    """Reproducible random number generator seeded at 42."""
    return np.random.default_rng(42)


def load_reference(name: str) -> dict:
    """Load a JSON reference file from the testdata directory.

    Parameters
    ----------
    name : str
        Filename (e.g. ``'reference_siso_bt.json'``).

    Returns
    -------
    dict
        Parsed JSON content.
    """
    path = TESTDATA / name
    with open(path) as f:
        return json.load(f)
