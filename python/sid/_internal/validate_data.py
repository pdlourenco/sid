# Copyright (c) 2026 Pedro Lourenco. All rights reserved.
# This code is released under the MIT License. See LICENSE file in the
# project root for full license information.
#
# This module is part of the Open Source System Identification Toolbox (SID).
# https://github.com/pdlourenco/sid-matlab

"""Validate and orient data for frequency-domain estimation functions."""

from __future__ import annotations

import warnings

import numpy as np

from sid._exceptions import SidError


def validate_data(
    y: np.ndarray | list,
    u: np.ndarray | list | None = None,
) -> tuple[np.ndarray, np.ndarray | None, int, int, int, bool, int]:
    """Validate and orient data for ``sid.freq_*`` functions.

    This is the Python port of ``sidValidateData.m``.  It ensures column
    orientation, checks for NaN/Inf and complex data, and verifies size
    consistency.  Multi-trajectory input is supported via 3-D arrays or
    lists of arrays (trimmed to the shortest trajectory length).

    Parameters
    ----------
    y : ndarray or list of ndarray
        Output data.  Accepted shapes: ``(N,)``, ``(N, ny)``,
        ``(N, ny, L)``, or a length-*L* list of 1-D / 2-D arrays.
    u : ndarray, list of ndarray, or None, optional
        Input data.  Same shape conventions as *y*.  Pass ``None`` (the
        default) for time-series (output-only) identification.

    Returns
    -------
    y : ndarray
        Oriented output data, shape ``(N, ny)`` or ``(N, ny, L)``.
    u : ndarray or None
        Oriented input data, or ``None`` for time-series mode.
    N : int
        Number of samples per trajectory.
    ny : int
        Number of output channels.
    nu : int
        Number of input channels (0 for time-series).
    is_time_series : bool
        ``True`` when *u* is ``None``.
    n_traj : int
        Number of trajectories (1 for single-trajectory).

    Raises
    ------
    SidError
        With codes ``'too_short'``, ``'complex_data'``, ``'non_finite'``,
        ``'dim_mismatch'``, ``'size_mismatch'``, ``'traj_mismatch'``, or
        ``'bad_input'``.

    Examples
    --------
    >>> y, u, N, ny, nu, is_ts, n_traj = validate_data(y, u)  # doctest: +SKIP

    Notes
    -----
    **Specification:** (Input validation -- not yet in SPEC.md)

    See Also
    --------
    sid.freq_bt : Main function that uses this validator.

    Changelog
    ---------
    2026-04-08 : First version by Pedro Lourenco.
    """

    # ---- Handle list input (variable-length trajectories) ---------------
    if isinstance(y, list):
        is_time_series = u is None or (isinstance(u, list) and len(u) == 0)
        L = len(y)
        if L == 0:
            raise SidError("bad_input", "Cell arrays must not be empty.")
        if not is_time_series:
            if not isinstance(u, list):
                raise SidError(
                    "bad_input",
                    "When y is a list, u must also be a list or None.",
                )
            if len(u) != L:
                raise SidError(
                    "dim_mismatch",
                    f"y has {L} trajectories but u has {len(u)}.",
                )

        # Ensure each trajectory is at least 2-D and find the shortest
        lengths = np.empty(L, dtype=int)
        for traj in range(L):
            y[traj] = np.asarray(y[traj], dtype=np.float64)
            if y[traj].ndim == 1:
                y[traj] = y[traj][:, np.newaxis]
            lengths[traj] = y[traj].shape[0]

        N_common = int(np.min(lengths))
        ny = y[0].shape[1]

        # Stack into 3-D array
        y_3d = np.zeros((N_common, ny, L), dtype=np.float64)
        for traj in range(L):
            if y[traj].shape[1] != ny:
                raise SidError(
                    "dim_mismatch",
                    f"y[{traj}] has {y[traj].shape[1]} columns, expected {ny}.",
                )
            y_3d[:, :, traj] = y[traj][:N_common, :]
        y = y_3d

        if not is_time_series:
            assert isinstance(u, list)  # for type checker
            nu = None
            for traj in range(L):
                u[traj] = np.asarray(u[traj], dtype=np.float64)
                if u[traj].ndim == 1:
                    u[traj] = u[traj][:, np.newaxis]
                if nu is None:
                    nu = u[traj].shape[1]
                elif u[traj].shape[1] != nu:
                    raise SidError(
                        "dim_mismatch",
                        f"u[{traj}] has {u[traj].shape[1]} columns, expected {nu}.",
                    )
                if u[traj].shape[0] < N_common:
                    raise SidError(
                        "size_mismatch",
                        f"u[{traj}] has {u[traj].shape[0]} samples but y requires "
                        f"at least {N_common}.",
                    )
            assert nu is not None
            u_3d = np.zeros((N_common, nu, L), dtype=np.float64)
            for traj in range(L):
                u_3d[:, :, traj] = u[traj][:N_common, :]
            u = u_3d

        if np.any(lengths != N_common):
            warnings.warn(
                f"Variable-length trajectories trimmed to shortest length N = {N_common}.",
                stacklevel=2,
            )
    else:
        y = np.asarray(y)
        if np.iscomplexobj(y):
            raise SidError(
                "complex_data",
                "Complex data is not supported in v1.0. Input y must be real.",
            )
        y = np.asarray(y, dtype=np.float64)

    # ---- Ensure column orientation --------------------------------------
    if y.ndim == 1:
        y = y[:, np.newaxis]

    is_time_series = u is None
    if not is_time_series:
        u = np.asarray(u, dtype=np.float64)
        if u.ndim == 1:
            u = u[:, np.newaxis]

    # ---- Detect multi-trajectory (3-D arrays) ---------------------------
    if y.ndim == 3:
        n_traj = y.shape[2]
    else:
        n_traj = 1

    N = y.shape[0]
    ny = y.shape[1]

    # ---- Validate data --------------------------------------------------
    if N < 2:
        raise SidError("too_short", "Data must have at least 2 samples.")
    if np.iscomplexobj(y):
        raise SidError(
            "complex_data",
            "Complex data is not supported in v1.0. Input y must be real.",
        )
    if not np.all(np.isfinite(y)):
        raise SidError("non_finite", "Data y contains NaN or Inf values.")

    if not is_time_series:
        assert u is not None
        nu = u.shape[1]
        if u.shape[0] != N:
            raise SidError(
                "size_mismatch",
                f"Input u ({u.shape[0]} samples) and output y ({N} samples) "
                f"must have the same length.",
            )
        if np.iscomplexobj(u):
            raise SidError(
                "complex_data",
                "Complex data is not supported in v1.0. Input u must be real.",
            )
        if not np.all(np.isfinite(u)):
            raise SidError("non_finite", "Data u contains NaN or Inf values.")
        # Multi-trajectory: u must have same number of trajectories
        if n_traj > 1:
            if u.ndim != 3 or u.shape[2] != n_traj:
                raise SidError(
                    "traj_mismatch",
                    f"y has {n_traj} trajectories but u does not match.",
                )
    else:
        nu = 0

    if N < 10:
        warnings.warn(
            f"Very short data (N = {N}). Estimates will be unreliable.",
            stacklevel=2,
        )

    return y, u, N, ny, nu, is_time_series, n_traj
