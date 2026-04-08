# Copyright (c) 2026 Pedro Lourenco. All rights reserved.
# This code is released under the MIT License. See LICENSE file in the
# project root for full license information.
#
# This module is part of the Open Source System Identification Toolbox (SID).
# https://github.com/pdlourenco/sid-matlab

"""Shared fixtures and helpers for sid test suite."""

from __future__ import annotations

import json
import pathlib

import numpy as np
import pytest

TESTDATA = pathlib.Path(__file__).resolve().parent.parent.parent / "testdata"


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
