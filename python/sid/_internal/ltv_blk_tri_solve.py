# Copyright (c) 2026 Pedro Lourenco. All rights reserved.
# This code is released under the MIT License. See LICENSE file in the
# project root for full license information.
#
# This module is part of the Open Source System Identification Toolbox (SID).
# https://github.com/pdlourenco/sid-matlab

"""Generic block tridiagonal forward-backward solver."""

from __future__ import annotations

import warnings

import numpy as np


def blk_tri_solve(
    S: list[np.ndarray],
    U_blk: list[np.ndarray],
    Theta: list[np.ndarray],
) -> tuple[list[np.ndarray], list[np.ndarray]]:
    """Solve a symmetric block tridiagonal linear system.

    This is the Python port of ``sidLTVblkTriSolve.m``.

    Solves the system using Gaussian elimination (forward pass) followed
    by back-substitution (backward pass).  Supports non-uniform block
    sizes via lists of arrays.

    The system has the form::

        [ S[0]    U[0]                          ] [ w[0]   ]   [ Theta[0]   ]
        [ U[0]'   S[1]    U[1]                  ] [ w[1]   ]   [ Theta[1]   ]
        [         U[1]'   S[2]   ...             ] [ ...    ] = [ ...        ]
        [                  ...   ...   U[K-2]    ] [ ...    ]   [ ...        ]
        [                        U[K-2]'  S[K-1] ] [ w[K-1] ]   [ Theta[K-1] ]

    where the sub-diagonal blocks are ``U[k].T`` (symmetric Hessian).

    Parameters
    ----------
    S : list of ndarray
        Length-*K* list.  ``S[k]`` is a ``(m_k, m_k)`` symmetric
        positive-definite diagonal block.
    U_blk : list of ndarray
        Length-*(K-1)* list.  ``U_blk[k]`` is a ``(m_k, m_{k+1})``
        super-diagonal block coupling step *k* to step *k+1*.
        The sub-diagonal is ``U_blk[k].T`` by symmetry.
    Theta : list of ndarray
        Length-*K* list.  ``Theta[k]`` is a ``(m_k, n_rhs)``
        right-hand side.

    Returns
    -------
    w : list of ndarray
        Length-*K* list.  ``w[k]`` is a ``(m_k, n_rhs)`` solution.
    Lbd : list of ndarray
        Length-*K* list.  ``Lbd[k]`` is a ``(m_k, m_k)`` forward
        Schur complement (stored for potential reuse in uncertainty).

    Examples
    --------
    >>> import numpy as np
    >>> K = 5; m = 3
    >>> S = [np.eye(m) * (k + 2) for k in range(K)]
    >>> U_blk = [0.1 * np.eye(m) for _ in range(K - 1)]
    >>> Theta = [np.ones((m, 1)) for _ in range(K)]
    >>> w, Lbd = blk_tri_solve(S, U_blk, Theta)
    >>> len(w)
    5

    Notes
    -----
    **Specification:** docs/cosmic_output.md -- Appendices A and B

    **Algorithm:**

    Forward pass (Gaussian elimination)::

        Lbd[0] = S[0],  Y[0] = Lbd[0] \\ Theta[0]
        Lbd[k] = S[k] - U[k-1].T @ (Lbd[k-1] \\ U[k-1])
        Y[k]   = Lbd[k] \\ (Theta[k] - U[k-1].T @ Y[k-1])

    Backward pass (back-substitution)::

        w[K-1] = Y[K-1]
        w[k]   = Y[k] - Lbd[k] \\ (U[k] @ w[k+1])

    Complexity: O(K * m^3) where m is the typical block size.

    See Also
    --------
    sid._internal.ltv_cosmic_solve.cosmic_solve :
        Specialised COSMIC block tridiagonal solver.
    sid.ltv_state_est.ltv_state_est :
        Uses this solver for RTS smoothing.

    Changelog
    ---------
    2026-04-09 : First version (Python port) by Pedro Lourenco.
    """

    K = len(S)

    Lbd: list[np.ndarray | None] = [None] * K
    Y: list[np.ndarray | None] = [None] * K
    w: list[np.ndarray | None] = [None] * K

    # ---- Forward pass: Gaussian elimination ------------------------------
    Lbd[0] = S[0].copy()
    Y[0] = np.linalg.solve(Lbd[0], Theta[0])

    for k in range(1, K):
        rc = np.linalg.cond(Lbd[k - 1])
        if 1.0 / rc < np.finfo(float).eps:
            warnings.warn(
                f"Block tridiagonal forward pass: Lbd[{k - 1}] is "
                f"near-singular (rcond={1.0 / rc:.2e}). "
                "Results may be unreliable.",
                stacklevel=2,
            )
        LbdInvU = np.linalg.solve(Lbd[k - 1], U_blk[k - 1])
        Lbd[k] = S[k] - U_blk[k - 1].T @ LbdInvU
        Y[k] = np.linalg.solve(Lbd[k], Theta[k] - U_blk[k - 1].T @ Y[k - 1])

    # ---- Backward pass: back-substitution --------------------------------
    w[K - 1] = Y[K - 1].copy()

    for k in range(K - 2, -1, -1):
        w[k] = Y[k] - np.linalg.solve(Lbd[k], U_blk[k] @ w[k + 1])

    return w, Lbd  # type: ignore[return-value]
