# Copyright (c) 2026 Pedro Lourenco. All rights reserved.
# This code is released under the MIT License. See LICENSE file in the
# project root for full license information.
#
# This module is part of the Open Source System Identification Toolbox (SID).
# https://github.com/pdlourenco/sid

"""Spring-mass-damper (SMD) plants for the sid example notebooks.

This module provides physically-interpretable plants used by the Jupyter
notebooks in this directory. Three variants are supported:

- :func:`util_msd` — LTI n-mass chain (wall--k1--m1--k2--m2--...--kn--mn)
  discretized via exact zero-order hold.
- :func:`util_msd_ltv` — LTV n-mass chain with per-time-step parameter
  trajectories (continuous ramps, step changes, sinusoidal variations, etc.).
- :func:`util_msd_nl` — Nonlinear simulation with Duffing-style cubic
  stiffness, integrated via fixed-step RK4 with zero-order-hold inputs.

The chain topology is the classic "wall + n masses in series" arrangement:

    wall --k1,c1-- m1 --k2,c2-- m2 --k3,c3-- ... --kn,cn-- mn

with state vector ``[x1, x2, ..., xn, v1, v2, ..., vn]`` (n positions
followed by n velocities) and base units kg / N·m⁻¹ / N·s·m⁻¹ / s.

These functions are sibling utilities to the notebooks — import them with a
plain ``from util_msd import ...`` when the current working directory is
``python/examples/`` (the default under ``pytest --nbmake`` and when Jupyter
is launched from this directory).
"""

from __future__ import annotations

import numpy as np
from scipy.linalg import expm


# ----------------------------------------------------------------------
# Internal helpers
# ----------------------------------------------------------------------


def _build_K_C(k_spring: np.ndarray, c_damp: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Build the stiffness and damping tridiagonal matrices for an n-mass chain.

    Spring ``i`` connects the wall (``i=0``) or mass ``i-1`` (``i>=1``) to
    mass ``i``. The resulting K and C matrices have identical structure.
    """
    n = len(k_spring)
    K = np.zeros((n, n))
    C = np.zeros((n, n))
    for i in range(n):
        K[i, i] = k_spring[i] + (k_spring[i + 1] if i + 1 < n else 0.0)
        C[i, i] = c_damp[i] + (c_damp[i + 1] if i + 1 < n else 0.0)
        if i + 1 < n:
            K[i, i + 1] = -k_spring[i + 1]
            K[i + 1, i] = -k_spring[i + 1]
            C[i, i + 1] = -c_damp[i + 1]
            C[i + 1, i] = -c_damp[i + 1]
    return K, C


def _continuous_ss(
    m: np.ndarray,
    k_spring: np.ndarray,
    c_damp: np.ndarray,
    F: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Assemble the continuous-time state-space matrices ``(Ac, Bc)``.

    ``Ac`` has shape ``(2n, 2n)`` and ``Bc`` has shape ``(2n, q)`` where
    ``q = F.shape[1]``. The state layout is ``[positions; velocities]``.
    """
    n = len(m)
    K, C = _build_K_C(k_spring, c_damp)
    M_mat = np.diag(m)
    Minv = np.linalg.solve(M_mat, np.eye(n))
    Ac = np.block([
        [np.zeros((n, n)), np.eye(n)],
        [-Minv @ K, -Minv @ C],
    ])
    Bc = np.block([
        [np.zeros((n, F.shape[1]))],
        [Minv @ F],
    ])
    return Ac, Bc


# ----------------------------------------------------------------------
# Public API
# ----------------------------------------------------------------------


def util_msd(
    m: np.ndarray,
    k_spring: np.ndarray,
    c_damp: np.ndarray,
    F: np.ndarray,
    Ts: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Discretize an LTI n-mass spring-damper chain with exact zero-order hold.

    Parameters
    ----------
    m : ndarray, shape ``(n,)``
        Masses (kg). ``n`` is inferred from ``len(m)``.
    k_spring : ndarray, shape ``(n,)``
        Spring constants (N/m). ``k_spring[0]`` is the wall-to-mass-1 spring,
        ``k_spring[i]`` (``i>=1``) connects mass ``i`` to mass ``i+1``.
    c_damp : ndarray, shape ``(n,)``
        Damping coefficients (N·s/m), same index convention as ``k_spring``.
    F : ndarray, shape ``(n, q)``
        Force input distribution matrix. Column ``j`` specifies how input
        ``u_j`` is distributed across the masses.
    Ts : float
        Sample time (s).

    Returns
    -------
    Ad : ndarray, shape ``(2n, 2n)``
        Discrete dynamics matrix.
    Bd : ndarray, shape ``(2n, q)``
        Discrete input matrix.

    Examples
    --------
    Single-DoF oscillator (resonance at ``sqrt(k/m) = 10`` rad/s):

    >>> import numpy as np
    >>> Ad, Bd = util_msd(
    ...     m=np.array([1.0]),
    ...     k_spring=np.array([100.0]),
    ...     c_damp=np.array([0.5]),
    ...     F=np.array([[1.0]]),
    ...     Ts=0.01,
    ... )
    >>> Ad.shape, Bd.shape
    ((2, 2), (2, 1))
    """
    m = np.asarray(m, dtype=np.float64)
    k_spring = np.asarray(k_spring, dtype=np.float64)
    c_damp = np.asarray(c_damp, dtype=np.float64)
    F = np.asarray(F, dtype=np.float64)

    if m.ndim != 1:
        raise ValueError(f"m must be 1-D, got shape {m.shape}")
    if k_spring.shape != m.shape:
        raise ValueError(f"k_spring shape {k_spring.shape} must match m shape {m.shape}")
    if c_damp.shape != m.shape:
        raise ValueError(f"c_damp shape {c_damp.shape} must match m shape {m.shape}")
    if F.ndim == 1:
        F = F.reshape(-1, 1)
    if F.shape[0] != m.shape[0]:
        raise ValueError(f"F first dim {F.shape[0]} must equal n = {m.shape[0]}")

    Ac, Bc = _continuous_ss(m, k_spring, c_damp, F)
    ns = Ac.shape[0]
    Ad = expm(Ac * Ts)
    Bd = np.linalg.solve(Ac, (Ad - np.eye(ns))) @ Bc
    return Ad, Bd


def util_msd_ltv(
    m: np.ndarray,
    k_spring: np.ndarray,
    c_damp: np.ndarray,
    F: np.ndarray,
    Ts: float,
    N: int | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """Discretize an LTV n-mass spring-damper chain with per-step parameters.

    Any of ``m``, ``k_spring``, ``c_damp`` may be 1-D (constant in time) or
    2-D with shape ``(n, N)`` (time-varying). ``F`` may be ``(n, q)`` or
    ``(n, q, N)``. ``N`` is inferred from the first 2-D input; if all inputs
    are 1-D, ``N`` must be passed explicitly and the LTI result is replicated
    ``N`` times.

    Parameters
    ----------
    m : ndarray, shape ``(n,)`` or ``(n, N)``
        Masses (kg).
    k_spring : ndarray, shape ``(n,)`` or ``(n, N)``
        Spring constants (N/m).
    c_damp : ndarray, shape ``(n,)`` or ``(n, N)``
        Damping coefficients (N·s/m).
    F : ndarray, shape ``(n, q)`` or ``(n, q, N)``
        Force input distribution matrix (optionally time-varying).
    Ts : float
        Sample time (s).
    N : int, optional
        Number of time steps. Required only when all other inputs are
        time-invariant (1-D / 2-D).

    Returns
    -------
    Ad : ndarray, shape ``(2n, 2n, N)``
    Bd : ndarray, shape ``(2n, q, N)``

    Examples
    --------
    2-mass chain with linearly ramping first-spring stiffness over 200 steps:

    >>> import numpy as np
    >>> N = 200
    >>> m = np.array([1.0, 1.0])
    >>> c = np.array([0.5, 0.5])
    >>> k = np.zeros((2, N))
    >>> k[0, :] = np.linspace(100, 150, N)
    >>> k[1, :] = 80.0
    >>> F = np.array([[1.0], [0.0]])
    >>> Ad, Bd = util_msd_ltv(m, k, c, F, Ts=0.01)
    >>> Ad.shape, Bd.shape
    ((4, 4, 200), (4, 1, 200))
    """
    m = np.asarray(m, dtype=np.float64)
    k_spring = np.asarray(k_spring, dtype=np.float64)
    c_damp = np.asarray(c_damp, dtype=np.float64)
    F = np.asarray(F, dtype=np.float64)

    def _n_from(arr: np.ndarray, name: str) -> int:
        if arr.ndim == 1:
            return arr.shape[0]
        if arr.ndim == 2:
            return arr.shape[0]
        raise ValueError(f"{name} must be 1-D or 2-D, got shape {arr.shape}")

    n = _n_from(m, "m")
    if _n_from(k_spring, "k_spring") != n:
        raise ValueError("k_spring and m must have matching first dimension")
    if _n_from(c_damp, "c_damp") != n:
        raise ValueError("c_damp and m must have matching first dimension")
    if F.ndim == 2:
        if F.shape[0] != n:
            raise ValueError(f"F first dim {F.shape[0]} must equal n = {n}")
        F_is_tv = False
    elif F.ndim == 3:
        if F.shape[0] != n:
            raise ValueError(f"F first dim {F.shape[0]} must equal n = {n}")
        F_is_tv = True
    else:
        raise ValueError(f"F must be 2-D or 3-D, got shape {F.shape}")

    # Infer N from the first 2-D / 3-D input
    inferred_N: int | None = None

    def _check_N(candidate: int, source: str) -> None:
        nonlocal inferred_N
        if inferred_N is None:
            inferred_N = candidate
        elif candidate != inferred_N:
            raise ValueError(
                f"Inconsistent N: {source}={candidate} vs earlier={inferred_N}"
            )

    if m.ndim == 2:
        _check_N(m.shape[1], "m")
    if k_spring.ndim == 2:
        _check_N(k_spring.shape[1], "k_spring")
    if c_damp.ndim == 2:
        _check_N(c_damp.shape[1], "c_damp")
    if F_is_tv:
        _check_N(F.shape[2], "F")

    if inferred_N is None:
        if N is None:
            raise ValueError(
                "All inputs are time-invariant; pass N explicitly to replicate "
                "the LTI result across time steps"
            )
        inferred_N = N
    elif N is not None and N != inferred_N:
        raise ValueError(f"Explicit N={N} conflicts with inferred N={inferred_N}")

    q = F.shape[1]
    Ad_seq = np.zeros((2 * n, 2 * n, inferred_N))
    Bd_seq = np.zeros((2 * n, q, inferred_N))

    # Fast path when nothing actually varies: compute once, broadcast.
    if m.ndim == 1 and k_spring.ndim == 1 and c_damp.ndim == 1 and not F_is_tv:
        Ad, Bd = util_msd(m, k_spring, c_damp, F, Ts)
        Ad_seq[:] = Ad[:, :, None]
        Bd_seq[:] = Bd[:, :, None]
        return Ad_seq, Bd_seq

    for k in range(inferred_N):
        m_k = m if m.ndim == 1 else m[:, k]
        k_k = k_spring if k_spring.ndim == 1 else k_spring[:, k]
        c_k = c_damp if c_damp.ndim == 1 else c_damp[:, k]
        F_k = F if not F_is_tv else F[:, :, k]
        Ad_seq[:, :, k], Bd_seq[:, :, k] = util_msd(m_k, k_k, c_k, F_k, Ts)

    return Ad_seq, Bd_seq


def util_msd_nl(
    m: np.ndarray,
    k_lin: np.ndarray,
    k_cubic: np.ndarray,
    c_damp: np.ndarray,
    F: np.ndarray,
    Ts: float,
    u: np.ndarray,
    x0: np.ndarray | None = None,
    substeps: int = 1,
) -> np.ndarray:
    """Simulate an n-mass chain with Duffing-style cubic stiffness via RK4.

    Each spring's restoring force is

        f_spring_i = k_lin[i] * delta_i + k_cubic[i] * delta_i**3

    where ``delta_0 = x[0]`` (wall-to-mass-1 spring) and
    ``delta_i = x[i] - x[i-1]`` for ``i >= 1``. Damping is linear:
    ``f_damp_i = c_damp[i] * d(delta_i)/dt``. The input ``u`` is held
    constant over each sample interval (zero-order hold).

    Integration uses classic fixed-step RK4 at step ``Ts / substeps``.
    Setting ``substeps > 1`` tightens the local truncation error at the cost
    of proportionally more RHS evaluations per sample.

    Parameters
    ----------
    m : ndarray, shape ``(n,)``
        Masses (kg).
    k_lin : ndarray, shape ``(n,)``
        Linear spring constants (N/m).
    k_cubic : ndarray, shape ``(n,)``
        Cubic stiffness coefficients (N/m³). Pass zeros for a linear plant;
        use ``k_cubic[0] > 0`` for a hardening Duffing oscillator on the
        wall-to-mass-1 spring.
    c_damp : ndarray, shape ``(n,)``
        Damping coefficients (N·s/m).
    F : ndarray, shape ``(n, q)``
        Force input distribution matrix.
    Ts : float
        Sample time (s).
    u : ndarray, shape ``(N, q)`` or ``(N,)``
        Input signal (``q`` inputs, ``N`` samples).
    x0 : ndarray, shape ``(2n,)``, optional
        Initial state ``[positions; velocities]``. Defaults to zeros.
    substeps : int, default 1
        RK4 sub-steps per sample ``Ts``. Use ``>=4`` for strong
        nonlinearities or large ``Ts``.

    Returns
    -------
    x : ndarray, shape ``(N + 1, 2n)``
        State trajectory. Row ``k`` is the state at time ``k * Ts``.

    Examples
    --------
    SDOF Duffing oscillator driven by white noise:

    >>> import numpy as np
    >>> rng = np.random.default_rng(0)
    >>> N = 500
    >>> u = rng.standard_normal((N, 1))
    >>> x = util_msd_nl(
    ...     m=np.array([1.0]),
    ...     k_lin=np.array([100.0]),
    ...     k_cubic=np.array([1000.0]),
    ...     c_damp=np.array([0.5]),
    ...     F=np.array([[1.0]]),
    ...     Ts=0.01,
    ...     u=u,
    ...     substeps=4,
    ... )
    >>> x.shape
    (501, 2)
    """
    m = np.asarray(m, dtype=np.float64)
    k_lin = np.asarray(k_lin, dtype=np.float64)
    k_cubic = np.asarray(k_cubic, dtype=np.float64)
    c_damp = np.asarray(c_damp, dtype=np.float64)
    F = np.asarray(F, dtype=np.float64)
    u = np.asarray(u, dtype=np.float64)

    n = len(m)
    if k_lin.shape != m.shape:
        raise ValueError("k_lin shape must match m")
    if k_cubic.shape != m.shape:
        raise ValueError("k_cubic shape must match m")
    if c_damp.shape != m.shape:
        raise ValueError("c_damp shape must match m")
    if F.ndim == 1:
        F = F.reshape(-1, 1)
    if F.shape[0] != n:
        raise ValueError(f"F first dim {F.shape[0]} must equal n = {n}")
    if u.ndim == 1:
        u = u.reshape(-1, 1)
    if u.shape[1] != F.shape[1]:
        raise ValueError(
            f"u has {u.shape[1]} inputs but F has {F.shape[1]} columns"
        )
    if substeps < 1:
        raise ValueError(f"substeps must be >= 1, got {substeps}")

    N_samples = u.shape[0]
    ns = 2 * n

    if x0 is None:
        x = np.zeros(ns)
    else:
        x0 = np.asarray(x0, dtype=np.float64)
        if x0.shape != (ns,):
            raise ValueError(f"x0 shape must be ({ns},), got {x0.shape}")
        x = x0.copy()

    inv_m = 1.0 / m

    def rhs(state: np.ndarray, force_vec: np.ndarray) -> np.ndarray:
        pos = state[:n]
        vel = state[n:]
        net_force = np.zeros(n)
        # Spring 0: wall to mass 0
        delta0 = pos[0]
        dvel0 = vel[0]
        f0 = k_lin[0] * delta0 + k_cubic[0] * delta0**3 + c_damp[0] * dvel0
        net_force[0] -= f0
        # Springs i=1..n-1: between mass i-1 and mass i
        for i in range(1, n):
            delta = pos[i] - pos[i - 1]
            dvel = vel[i] - vel[i - 1]
            fi = k_lin[i] * delta + k_cubic[i] * delta**3 + c_damp[i] * dvel
            net_force[i - 1] += fi
            net_force[i] -= fi
        # External force (ZOH: constant over the sample interval)
        net_force += F @ force_vec
        acc = inv_m * net_force
        return np.concatenate([vel, acc])

    h = Ts / substeps
    x_traj = np.zeros((N_samples + 1, ns))
    x_traj[0] = x
    for k in range(N_samples):
        uk = u[k]
        for _ in range(substeps):
            k1 = rhs(x, uk)
            k2 = rhs(x + 0.5 * h * k1, uk)
            k3 = rhs(x + 0.5 * h * k2, uk)
            k4 = rhs(x + h * k3, uk)
            x = x + (h / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)
        x_traj[k + 1] = x

    return x_traj
