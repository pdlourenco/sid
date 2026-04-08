# Copyright (c) 2026 Pedro Lourenco. All rights reserved.
# This code is released under the MIT License. See LICENSE file in the
# project root for full license information.
#
# This module is part of the Open Source System Identification Toolbox (SID).
# https://github.com/pdlourenco/sid-matlab

"""Custom exception classes for the sid toolbox."""

from __future__ import annotations


class SidError(Exception):
    """Base exception for sid toolbox errors.

    Parameters
    ----------
    code : str
        Machine-readable error code (e.g. ``'too_short'``, ``'non_finite'``).
    message : str
        Human-readable error description.

    Attributes
    ----------
    code : str
        The error code passed at construction time.
    """

    def __init__(self, code: str, message: str) -> None:
        self.code = code
        super().__init__(message)
