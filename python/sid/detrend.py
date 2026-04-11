# Copyright (c) 2026 Pedro Lourenco. All rights reserved.
# This code is released under the MIT License. See LICENSE file in the
# project root for full license information.
#
# This module is part of the Open Source System Identification Toolbox (SID).
# https://github.com/pdlourenco/sid

"""Polynomial detrending for time-domain data."""

from __future__ import annotations

import warnings

import numpy as np

from sid._exceptions import SidError


def detrend(
    x: np.ndarray,
    *,
    order: int = 1,
    segment_length: int | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """Remove a polynomial trend from time-domain data.

    Fits and subtracts a polynomial trend of degree *order* from each
    channel (and each trajectory) of *x*.  Detrending is standard
    preprocessing before spectral estimation: unremoved trends bias
    low-frequency spectral estimates and violate stationarity
    assumptions.

    Parameters
    ----------
    x : ndarray, shape (N,), (N, n_ch), or (N, n_ch, L)
        Real-valued time-domain data.  A 1-D array is treated as a
        single-channel column signal and the returned arrays are
        squeezed back to 1-D.  For multi-trajectory data, pass a 3-D
        array ``(N, n_ch, L)`` where *L* is the number of trajectories.
    order : int, optional
        Polynomial degree to remove.  ``0`` removes the mean, ``1``
        removes a linear trend, ``2`` removes a quadratic trend, etc.
        Must be a non-negative integer.  Default is ``1``.
    segment_length : int, optional
        When set, each non-overlapping segment of this length is
        detrended independently.  Must be a positive integer.  If
        *segment_length* exceeds *N* it is silently clamped to *N*.
        Default is ``N`` (full-record detrend).

    Returns
    -------
    x_detrended : ndarray, same shape as *x*
        The detrended data, ``x - trend``.
    trend : ndarray, same shape as *x*
        The polynomial trend that was removed, so that
        ``x == x_detrended + trend`` (up to floating-point rounding).

    Raises
    ------
    SidError
        If *x* is complex (code: ``'complex_data'``).
    SidError
        If *x* contains NaN or Inf values (code: ``'non_finite'``).
    SidError
        If *order* is not a non-negative integer
        (code: ``'bad_order'``).
    SidError
        If *segment_length* is not a positive integer
        (code: ``'bad_segment_length'``).

    Examples
    --------
    Remove a linear trend (the default):

    >>> import numpy as np
    >>> from sid.detrend import detrend
    >>> t = np.arange(100, dtype=float)
    >>> x = 3.0 * t + np.random.default_rng(0).standard_normal(100)
    >>> x_dt, trend = detrend(x)
    >>> x_dt.shape
    (100,)

    Remove the mean only:

    >>> x_dm, _ = detrend(x, order=0)

    Segment-wise detrending:

    >>> x_ds, _ = detrend(x, segment_length=50)

    Multi-channel data:

    >>> x2d = np.column_stack([t, -2.0 * t])
    >>> x_dt2, trend2 = detrend(x2d)
    >>> x_dt2.shape
    (100, 2)

    Typical preprocessing workflow:

    >>> import sid  # doctest: +SKIP
    >>> y_dt, _ = detrend(y)  # doctest: +SKIP
    >>> u_dt, _ = detrend(u)  # doctest: +SKIP
    >>> result = sid.freq_bt(y_dt, u_dt)  # doctest: +SKIP

    Notes
    -----
    **Algorithm:**

    1. Validate inputs (real, finite, correct parameter types).
    2. For each trajectory and channel, split the signal into
       non-overlapping segments of length *segment_length*.
    3. In each segment, fit a polynomial of degree
       ``min(order, segment_length - 1)`` using least-squares
       (``numpy.polyfit``) and evaluate it (``numpy.polyval``) to obtain
       the local trend.
    4. Subtract the assembled trend from the original data.

    If a segment is too short for the requested polynomial order, the
    order is automatically reduced and a warning is issued.

    **Specification:** Data preprocessing -- not yet in SPEC.md.

    See Also
    --------
    sid.freq_bt : Blackman-Tukey spectral analysis.
    sid.freq_etfe : Empirical Transfer Function Estimate.
    sid.freq_map : Maximum a-posteriori frequency response estimation.

    Changelog
    ---------
    2026-04-08 : First version (Python port) by Pedro Lourenco.
    """
    x = np.asarray(x)

    # ---- Validate data ----
    if np.iscomplexobj(x):
        raise SidError("complex_data", "Input x must be real.")

    x = x.astype(np.float64, copy=False)

    # ---- Handle 1-D (vector) input ----
    was_1d = x.ndim == 1
    if was_1d:
        x = x[:, np.newaxis]
    if not np.all(np.isfinite(x)):
        raise SidError("non_finite", "Input x contains NaN or Inf values.")

    # ---- Dimensions ----
    N = x.shape[0]
    n_ch = x.shape[1]
    if x.ndim == 3:
        n_traj = x.shape[2]
    else:
        n_traj = 1

    # ---- Validate parameters ----
    if not isinstance(order, (int, np.integer)) or order < 0:
        raise SidError("bad_order", "Order must be a non-negative integer.")

    seg_len: int
    if segment_length is None:
        seg_len = N
    else:
        if not isinstance(segment_length, (int, np.integer)) or segment_length < 1:
            raise SidError(
                "bad_segment_length",
                "SegmentLength must be a positive integer.",
            )
        seg_len = int(segment_length)

    if seg_len > N:
        seg_len = N

    # ---- Detrend ----
    trend = np.zeros_like(x)

    for lt in range(n_traj):
        for ch in range(n_ch):
            # Extract column
            if n_traj > 1:
                col = x[:, ch, lt]
            else:
                col = x[:, ch]

            trend_col = np.zeros(N, dtype=np.float64)

            # Process each segment
            idx = 0
            while idx < N:
                seg_end = min(idx + seg_len, N)
                seg_n = seg_end - idx
                t = np.arange(seg_n, dtype=np.float64)

                seg = col[idx:seg_end]
                actual_order = min(order, seg_n - 1)
                if actual_order < order:
                    warnings.warn(
                        f"Segment of length {seg_n} is too short for "
                        f"polynomial order {order}. Reduced to order "
                        f"{actual_order}.",
                        stacklevel=2,
                    )
                coeffs = np.polyfit(t, seg, actual_order)
                trend_col[idx:seg_end] = np.polyval(coeffs, t)

                idx = seg_end

            if n_traj > 1:
                trend[:, ch, lt] = trend_col
            else:
                trend[:, ch] = trend_col

    x_detrended = x - trend

    # ---- Squeeze back to 1-D if input was 1-D ----
    if was_1d:
        x_detrended = x_detrended.ravel()
        trend = trend.ravel()

    return x_detrended, trend
