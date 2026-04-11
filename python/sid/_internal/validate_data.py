# Copyright (c) 2026 Pedro Lourenco. All rights reserved.
# This code is released under the MIT License. See LICENSE file in the
# project root for full license information.
#
# This module is part of the Open Source System Identification Toolbox (SID).
# https://github.com/pdlourenco/sid

"""Validate and orient data for frequency-domain estimation functions."""

from __future__ import annotations

import warnings

import numpy as np

from sid._exceptions import SidError


def validate_data(
    y: np.ndarray | list,
    u: np.ndarray | list | None = None,
    *,
    preserve_lengths: bool = False,
) -> tuple[np.ndarray | list, np.ndarray | list | None, int, int, int, bool, int]:
    """Validate and orient data for ``sid.freq_*`` functions.

    This is the Python port of ``sidValidateData.m``.  It ensures column
    orientation, checks for NaN/Inf and complex data, and verifies size
    consistency.  Multi-trajectory input is supported via 3-D arrays or
    lists of arrays.

    By default, variable-length list inputs are trimmed to the shortest
    trajectory length and a warning is emitted.  When
    ``preserve_lengths=True``, list inputs with different lengths are
    returned as lists of per-trajectory arrays instead of being trimmed;
    callers are then responsible for per-segment filtering (SPEC.md §6.2
    for ``sidFreqMap``).

    Parameters
    ----------
    y : ndarray or list of ndarray
        Output data.  Accepted shapes: ``(N,)``, ``(N, ny)``,
        ``(N, ny, L)``, or a length-*L* list of 1-D / 2-D arrays.
    u : ndarray, list of ndarray, or None, optional
        Input data.  Same shape conventions as *y*.  Pass ``None`` (the
        default) for time-series (output-only) identification.
    preserve_lengths : bool, optional
        When ``True`` and the input is a list of variable-length
        trajectories, the result is returned as a list of per-trajectory
        arrays rather than being trimmed to the shortest length.  This
        is used by :func:`sid.freq_map` to implement the per-segment
        filtering required by SPEC.md §6.2.  Default: ``False``.

    Returns
    -------
    y : ndarray or list of ndarray
        Oriented output data.  By default, a 2-D/3-D ndarray.  When
        ``preserve_lengths=True`` and the input was a variable-length
        list, a list of ``(N_l, ny)`` arrays (one per trajectory).
    u : ndarray, list of ndarray, or None
        Oriented input data, with the same shape convention as *y*, or
        ``None`` for time-series mode.
    N : int
        Number of samples per trajectory.  For preserved variable-length
        lists this is ``max(N_l)``; otherwise the common trimmed length.
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
    2026-04-10 : Add ``preserve_lengths`` for SPEC.md §6.2 compliance.
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

        # Ensure each trajectory is at least 2-D and find the per-trajectory length
        y_list: list[np.ndarray] = []
        lengths = np.empty(L, dtype=int)
        for traj in range(L):
            yt = np.asarray(y[traj])
            if np.iscomplexobj(yt):
                raise SidError(
                    "complex_data",
                    "Complex data is not supported in v1.0. Input y must be real.",
                )
            yt = yt.astype(np.float64)
            if yt.ndim == 1:
                yt = yt[:, np.newaxis]
            if not np.all(np.isfinite(yt)):
                raise SidError("non_finite", f"Data y[{traj}] contains NaN or Inf.")
            y_list.append(yt)
            lengths[traj] = yt.shape[0]

        ny = y_list[0].shape[1]
        for traj in range(1, L):
            if y_list[traj].shape[1] != ny:
                raise SidError(
                    "dim_mismatch",
                    f"y[{traj}] has {y_list[traj].shape[1]} columns, expected {ny}.",
                )

        # Validate/prepare u_list if present
        u_list: list[np.ndarray] | None = None
        nu_raw: int = 0
        if not is_time_series:
            assert isinstance(u, list)  # for type checker
            u_list = []
            for traj in range(L):
                ut = np.asarray(u[traj])
                if np.iscomplexobj(ut):
                    raise SidError(
                        "complex_data",
                        "Complex data is not supported in v1.0. Input u must be real.",
                    )
                ut = ut.astype(np.float64)
                if ut.ndim == 1:
                    ut = ut[:, np.newaxis]
                if not np.all(np.isfinite(ut)):
                    raise SidError("non_finite", f"Data u[{traj}] contains NaN or Inf.")
                if traj == 0:
                    nu_raw = ut.shape[1]
                elif ut.shape[1] != nu_raw:
                    raise SidError(
                        "dim_mismatch",
                        f"u[{traj}] has {ut.shape[1]} columns, expected {nu_raw}.",
                    )
                if ut.shape[0] != lengths[traj]:
                    raise SidError(
                        "size_mismatch",
                        f"u[{traj}] has {ut.shape[0]} samples but y[{traj}] has "
                        f"{lengths[traj]}. Trajectory-wise lengths must match.",
                    )
                u_list.append(ut)

        variable_length = bool(np.any(lengths != lengths[0]))

        if preserve_lengths and variable_length:
            # List-preserving path: return per-trajectory arrays for callers
            # that do per-segment filtering (SPEC.md §6.2 for sidFreqMap).
            N_max = int(np.max(lengths))
            if N_max < 2:
                raise SidError("too_short", "Data must have at least 2 samples.")
            if N_max < 10:
                warnings.warn(
                    f"Very short data (N_max = {N_max}). Estimates will be unreliable.",
                    stacklevel=2,
                )
            return (
                y_list,
                u_list,
                N_max,
                ny,
                (nu_raw if not is_time_series else 0),
                is_time_series,
                L,
            )

        # Trim-to-shortest path (default, and only option when lengths are uniform)
        N_common = int(np.min(lengths))
        y_3d = np.zeros((N_common, ny, L), dtype=np.float64)
        for traj in range(L):
            y_3d[:, :, traj] = y_list[traj][:N_common, :]
        y = y_3d

        if not is_time_series:
            assert u_list is not None
            u_3d = np.zeros((N_common, nu_raw, L), dtype=np.float64)
            for traj in range(L):
                u_3d[:, :, traj] = u_list[traj][:N_common, :]
            u = u_3d

        if variable_length:
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
