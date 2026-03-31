# Output-COSMIC: Extension of COSMIC to Partial State Observation

## 1. Problem Statement

Consider a discrete-time linear time-varying (LTV) system with partial state observation:

$$
x(k+1) = A(k)\,x(k) + B(k)\,u(k), \qquad k = 0, \ldots, N-1
$$

$$
y(k) = H\,x(k)
$$

where $x(k) \in \mathbb{R}^n$ is the state, $u(k) \in \mathbb{R}^q$ is a known input, $y(k) \in \mathbb{R}^{p_y}$ is the measurement, and $H \in \mathbb{R}^{p_y \times n}$ is a known, time-invariant observation matrix. The system matrices $A(k) \in \mathbb{R}^{n \times n}$ and $B(k) \in \mathbb{R}^{n \times q}$ are unknown and to be estimated.

Standard COSMIC assumes $H = I$ (full state observation). This document extends COSMIC to the case $p_y \leq n$, where the state is only partially observed.

## 2. Joint Objective

Define the combined dynamics matrix $C(k) = \begin{bmatrix} A^\top(k) \\ B^\top(k) \end{bmatrix} \in \mathbb{R}^{(n+q) \times n}$ and the augmented data vector $d(k) = \begin{bmatrix} x(k) \\ u(k) \end{bmatrix} \in \mathbb{R}^{n+q}$.

The joint objective over the state sequence $\mathbf{x} = \{x(k)\}_{k=0}^{N}$ and the dynamics sequence $\mathbf{C} = \{C(k)\}_{k=0}^{N-1}$ is:

$$
J(\mathbf{x}, \mathbf{C}) = \sum_{k=0}^{N} \|y(k) - H\,x(k)\|^2_{R^{-1}} + \sum_{k=0}^{N-1} \|x(k+1) - A(k)\,x(k) - B(k)\,u(k)\|^2 + \lambda \sum_{k=1}^{N-1} \|C(k) - C(k-1)\|_F^2
$$

where $\|v\|^2_{R^{-1}} = v^\top R^{-1} v$ denotes the Mahalanobis norm with respect to the measurement noise covariance $R \in \mathbb{R}^{p_y \times p_y}$ (symmetric positive definite), and:

- **Observation fidelity** ($R^{-1}$ term): penalises deviation of the estimated state from measurements, weighted by the measurement information matrix $R^{-1}$. This naturally downweights noisy measurement channels and upweights precise ones. If $R$ is unknown, setting $R = I$ recovers unweighted least squares.
- **Dynamics fidelity** (middle term): penalises deviation from the estimated dynamics model, coupling $\mathbf{x}$ and $\mathbf{C}$.
- **Dynamics smoothness** ($\lambda$ term): the standard COSMIC regulariser, penalising rapid variation of the system matrices between consecutive time steps.

### 2.1 Recovery of Standard COSMIC

When $H = I$ (full state observation), taking $R \to 0$ (infinite measurement precision) forces $x(k) = y(k)$ exactly. The observation fidelity term vanishes (it is zero at the constraint), the state is no longer a free variable, and $J$ reduces to the standard COSMIC cost:

$$
J(\mathbf{C}) = \sum_{k=0}^{N-1} \|y(k+1) - C^\top(k)\,d(k)\|^2 + \lambda \sum_{k=1}^{N-1} \|C(k) - C(k-1)\|_F^2
$$

which is exactly the original COSMIC formulation of Carvalho et al. (2022). No additional hyperparameters are introduced in the fully-observed case.

## 3. Multi-Trajectory Extension

Assume $L$ trajectories are collected, each of length $N+1$ time steps, under the same LTV dynamics $\{A(k), B(k)\}_{k=0}^{N-1}$. Denote the state and output of trajectory $l$ at time $k$ as $x_l(k)$ and $y_l(k)$ respectively. The joint objective generalises to:

$$
J(\mathbf{X}, \mathbf{C}) = \sum_{l=1}^{L} \sum_{k=0}^{N} \|y_l(k) - H\,x_l(k)\|^2_{R^{-1}} + \sum_{l=1}^{L} \sum_{k=0}^{N-1} \|x_l(k+1) - A(k)\,x_l(k) - B(k)\,u_l(k)\|^2 + \lambda \sum_{k=1}^{N-1} \|C(k) - C(k-1)\|_F^2
$$

The dynamics smoothness term is shared across trajectories (single $\mathbf{C}$), while the observation and dynamics fidelity terms are summed over trajectories. Each trajectory has its own state sequence $\{x_l(k)\}$ but shares the same $\{C(k)\}$.

## 4. Alternating Minimisation Algorithm

The objective $J(\mathbf{x}, \mathbf{C})$ is non-convex jointly (due to the bilinear coupling $A(k)\,x(k)$), but is strictly convex in each block separately when the other is held fixed (given $\lambda > 0$ and $R$ positive definite). This motivates a two-block alternating minimisation scheme.

### 4.1 Initialisation: State and Input Matrix Estimation with $A = I$

Before entering the alternating loop, an initialisation step estimates the state sequences and the input matrices $B(k)$ jointly, without requiring knowledge of $A(k)$. This is achieved by evaluating the joint objective $J$ at $A(k) = I$ for all $k$ and minimising over the remaining unknowns $\{\mathbf{x}, \mathbf{B}\}$.

Setting $A(k) = I$ in the multi-trajectory objective (Section 3), the smoothness term becomes $\lambda \sum_{k=1}^{N-1} \|C(k) - C(k-1)\|_F^2 = \lambda \sum_{k=1}^{N-1} \|B(k) - B(k-1)\|_F^2$ (since the $A = I$ block contributes zero). The initialisation cost is therefore:

$$
J_{\text{init}}(\mathbf{X}, \mathbf{B}) = J(\mathbf{X}, \mathbf{C})\big|_{A(k)=I} = \sum_{l=1}^{L} \sum_{k=0}^{N} \|y_l(k) - H\,x_l(k)\|^2_{R^{-1}} + \sum_{l=1}^{L} \sum_{k=0}^{N-1} \|x_l(k+1) - x_l(k) - B(k)\,u_l(k)\|^2 + \lambda \sum_{k=1}^{N-1} \|B(k) - B(k-1)\|_F^2
$$

This is jointly convex in $\{x_l(k)\}$ and $\{B(k)\}$. All terms are quadratic with no bilinear coupling between the two sets of unknowns: $B(k)\,u_l(k)$ is linear in $B(k)$ since $u_l(k)$ is known data, and $x_l(k)$ appears linearly. The $B(k)$ matrices are shared across trajectories (same LTV dynamics), while each trajectory has its own state sequence. The minimiser is unique and can be obtained in a single linear solve.

**Structure of the linear system.** The unknowns are the $n$-vectors $x_l(0), \ldots, x_l(N)$ for each trajectory $l$, and the $n \times q$ matrices $B(0), \ldots, B(N-1)$ shared across trajectories. The per-trajectory state unknowns are coupled temporally through the dynamics fidelity term and coupled to $B(k)$ through the input terms. The normal equations exhibit block-tridiagonal-plus-arrow structure.

**Physical interpretation.** The initialisation models the autonomous state evolution as a random walk ($A = I$) with input-driven corrections. The smoothness prior on $B(k)$ prevents the input matrix from absorbing dynamics that should be attributed to $A(k)$. This provides a state estimate that is robust to ignorance of $A(k)$, at the cost of a weaker dynamical model. Crucially, since $J_{\text{init}} = J\big|_{A=I}$, the initialisation is not a separate heuristic — it is the exact minimisation of the global objective over a restricted subspace.

### 4.2 Alternating Loop

Starting from the initialisation, the algorithm alternates between two steps:

**State Step.** Fix $\mathbf{C} = \{C(k)\}$, solve for $\mathbf{x} = \{x(k)\}$:

$$
\min_{\mathbf{x}} \; \sum_{k=0}^{N} \|y(k) - H\,x(k)\|^2_{R^{-1}} + \sum_{k=0}^{N-1} \|x(k+1) - A(k)\,x(k) - B(k)\,u(k)\|^2
$$

This is a linear least-squares problem in $\{x(k)\}$. It is exactly the Rauch–Tung–Striebel (RTS) smoother with measurement noise covariance $R$ and process noise covariance $Q = I$, conditioned on the full observation sequence $\{y(k)\}_{k=0}^{N}$. The solution is unique and can be computed in $O(N\,n^3)$ time via the standard forward–backward recursion.

For multi-trajectory: the state step decouples across trajectories (each trajectory has its own state sequence but uses the shared $\mathbf{C}$), so $L$ independent smoothers are run.

**COSMIC Step.** Fix $\mathbf{x} = \{x(k)\}$ (or $\mathbf{X} = \{X(k)\}$ for multi-trajectory), solve for $\mathbf{C}$:

$$
\min_{\mathbf{C}} \; \sum_{k=0}^{N-1} \|x(k+1) - C^\top(k)\,d(k)\|^2 + \lambda \sum_{k=1}^{N-1} \|C(k) - C(k-1)\|_F^2
$$

This is exactly the standard COSMIC problem with $x(k)$ replacing the measured state. The observation fidelity term does not involve $\mathbf{C}$ and is therefore constant in this step. The solution is obtained via the COSMIC tridiagonal closed-form algorithm in $O(N\,(n+q)^3)$ time.

For multi-trajectory: the data matrices $D(k) = \begin{bmatrix} X(k) \\ U(k) \end{bmatrix}$ and targets $X'(k) = X(k+1)$ pool across trajectories exactly as in the standard multi-trajectory COSMIC formulation.

### 4.3 Trust-Region Interpolation (Optional)

In cases where the transition from the initialisation ($A = I$) to the first COSMIC estimate of $A(k)$ is too abrupt — for instance with high noise, long trajectories, or poorly conditioned data — a trust-region mechanism can be employed.

Define the interpolated dynamics matrix:

$$
\tilde{A}(k) = (1 - \mu)\,A(k) + \mu\,I
$$

where $\mu \in [0, 1]$ is the trust-region parameter. The state step uses $\tilde{A}(k)$ in place of $A(k)$:

$$
\min_{\mathbf{x}} \; \sum_{k} \|y(k) - H\,x(k)\|^2_{R^{-1}} + \sum_{k} \|x(k+1) - \tilde{A}(k)\,x(k) - B(k)\,u(k)\|^2
$$

The COSMIC step is unaffected (it always solves for $A(k)$ and $B(k)$ freely).

**Adaptive schedule.** The trust-region parameter $\mu$ is managed by an outer loop wrapping the alternating minimisation:

1. Initialise $\mu = 1$ (which reduces the first state step to the $A = I$ initialisation).
2. Run the alternating state–COSMIC loop to convergence for the current $\mu$, yielding $J^*(\mu)$.
3. Reduce $\mu$: set $\mu \leftarrow \mu / 2$.
4. Run the alternating loop to convergence with the new $\mu$, yielding $J^*(\mu/2)$.
5. **Accept/reject:** If $J^*(\mu/2) \leq J^*(\mu)$, accept the reduction and continue from step 3. If $J^*(\mu/2) > J^*(\mu)$, increase $\mu$ back (set $\mu \leftarrow \min(2\mu, 1)$) and terminate.
6. Terminate when $\mu < \varepsilon$ (e.g., $\varepsilon = 10^{-6}$) and set $\mu = 0$ for a final pass.

At each value of $\mu$, the inner alternating loop minimises a well-defined objective (with $\tilde{A}(k)$ fixed for the duration of that inner loop), so the monotone decrease guarantee holds within each inner loop. The outer loop provides a monotone homotopy from the robust random-walk model ($\mu = 1$) to the fully dynamics-informed model ($\mu = 0$).

**When $\mu = 0$ throughout:** the trust-region mechanism is inactive and the algorithm reduces to the base two-block alternation described in Section 4.2, with a single initialisation step at $A = I$. This is expected to be sufficient for most practical cases.

## 5. Convergence Analysis

### 5.1 Monotone Decrease

**Proposition 1.** *Let $(\mathbf{x}^{(t)}, \mathbf{C}^{(t)})$ denote the iterates of the alternating minimisation (Section 4.2) at iteration $t$. Then the sequence $\{J(\mathbf{x}^{(t)}, \mathbf{C}^{(t)})\}$ is monotonically non-increasing.*

*Proof.* At each iteration:

1. The state step sets $\mathbf{x}^{(t+1)} = \arg\min_{\mathbf{x}} J(\mathbf{x}, \mathbf{C}^{(t)})$, so $J(\mathbf{x}^{(t+1)}, \mathbf{C}^{(t)}) \leq J(\mathbf{x}^{(t)}, \mathbf{C}^{(t)})$.
2. The COSMIC step sets $\mathbf{C}^{(t+1)} = \arg\min_{\mathbf{C}} J(\mathbf{x}^{(t+1)}, \mathbf{C})$, so $J(\mathbf{x}^{(t+1)}, \mathbf{C}^{(t+1)}) \leq J(\mathbf{x}^{(t+1)}, \mathbf{C}^{(t)})$.

Combining: $J(\mathbf{x}^{(t+1)}, \mathbf{C}^{(t+1)}) \leq J(\mathbf{x}^{(t)}, \mathbf{C}^{(t)})$. Since $J \geq 0$, the sequence converges. $\square$

### 5.2 Convergence to a Stationary Point

**Proposition 2.** *Every limit point of the sequence $\{(\mathbf{x}^{(t)}, \mathbf{C}^{(t)})\}$ is a stationary point of $J$.*

*Proof.* Both subproblems are strictly convex (the state step has $R^{-1} \succ 0$ ensuring strict convexity in $\mathbf{x}$; the COSMIC step has $\lambda > 0$ ensuring strict convexity in $\mathbf{C}$), so each block minimiser is unique. This satisfies the conditions of the two-block alternating minimisation convergence theorem (Grippo and Sciandrone, 2000, Theorem 2.1): for two-block problems where each block minimiser is unique, every limit point of the alternating sequence is a stationary point of the joint objective. $\square$

### 5.3 Non-Convexity and Global Optimality

The joint objective $J(\mathbf{x}, \mathbf{C})$ is non-convex due to the bilinear term $A(k)\,x(k)$ in the dynamics fidelity. Multiple stationary points may exist. In particular, a similarity transformation $T \in \mathbb{R}^{n \times n}$ (invertible) maps any solution $(x(k), A(k), B(k))$ to an equivalent solution $(T\,x(k),\; T\,A(k)\,T^{-1},\; T\,B(k))$ with identical cost. The observation term $\|y(k) - H\,x(k)\|^2_{R^{-1}}$ breaks this symmetry partially (requiring $H\,T^{-1}$ to produce the same outputs), but does not eliminate it unless $H$ has full column rank.

Global optimality is therefore not guaranteed. The initialisation (Section 4.1) and optional trust-region mechanism (Section 4.3) serve to place the iterates in a favourable basin of attraction.

### 5.4 Trust-Region Convergence

**Proposition 3.** *The outer trust-region loop (Section 4.3) produces a monotonically non-increasing sequence of converged objectives $\{J^*(\mu_s)\}$ as long as reductions are accepted.*

*Proof.* At each outer iteration $s$, the inner alternating loop converges to a stationary point for fixed $\mu_s$ (Propositions 1–2). A reduction $\mu_{s+1} = \mu_s / 2$ is accepted only if $J^*(\mu_{s+1}) \leq J^*(\mu_s)$. If the reduction is rejected, the algorithm reverts to $\mu_s$ and terminates. Therefore, the accepted sequence of $J^*$ values is monotonically non-increasing, and the algorithm terminates in a finite number of outer iterations (since $\mu$ halves at each accepted step and terminates below $\varepsilon$). $\square$

## 6. Similarity Transformation Ambiguity

As noted in Section 5.3, the output-COSMIC problem admits equivalent solutions under state-space similarity transformations. For any invertible $T$:

$$
\tilde{x}(k) = T\,x(k), \quad \tilde{A}(k) = T\,A(k)\,T^{-1}, \quad \tilde{B}(k) = T\,B(k), \quad \tilde{H} = H\,T^{-1}
$$

produces identical input-output behaviour. Since $H$ is fixed and known, the transformation is constrained by $H = \tilde{H}\,T = H\,T^{-1}\,T = H$, which is automatically satisfied — the point is that different initialisations may converge to different (but equivalent) state-space realisations.

This ambiguity is inherent to all output-based system identification methods and does not affect prediction or control design (which depend only on input-output behaviour). If a canonical form is desired, it can be imposed as a post-processing step (e.g., balanced realisation, observable canonical form).

## 7. Algorithm Summary

**Algorithm: Output-COSMIC**

**Input:** Measurements $\{y_l(k)\}$ for $l = 1, \ldots, L$, $k = 0, \ldots, N$; inputs $\{u_l(k)\}$ for $l = 1, \ldots, L$, $k = 0, \ldots, N-1$; observation matrix $H$; measurement noise covariance $R$ (or $R = I$ if unknown); regularisation weight $\lambda$; convergence tolerance $\varepsilon_J$; optional trust-region parameters $(\mu_0, \varepsilon_\mu)$.

**Output:** Estimated dynamics $\{A(k), B(k)\}_{k=0}^{N-1}$; estimated states $\{x_l(k)\}$.

**Initialisation:**
1. Solve $J_{\text{init}} = J\big|_{A=I}$ (Section 4.1) for $\{x_l(k)\}$ and $\{B(k)\}$ jointly, pooling across all $L$ trajectories. Store $\mathbf{X}^{(0)}$, $\mathbf{B}^{(0)}$. Set $A^{(0)}(k) = I$ for all $k$.

**Main loop** (iteration $t = 0, 1, 2, \ldots$):
2. **COSMIC step:** Using the current state estimates $\mathbf{X}^{(t)}$, solve the standard COSMIC problem (Section 4.2) pooling across trajectories. Obtain $\mathbf{C}^{(t+1)} = \{A^{(t+1)}(k), B^{(t+1)}(k)\}$.
3. **State step:** Using $\mathbf{C}^{(t+1)}$ (with trust-region interpolation $\tilde{A}(k) = (1-\mu)\,A^{(t+1)}(k) + \mu\,I$ if enabled), run the RTS smoother for each trajectory independently. Obtain $\mathbf{X}^{(t+1)}$.
4. Evaluate $J^{(t+1)} = J(\mathbf{X}^{(t+1)}, \mathbf{C}^{(t+1)})$.
5. If $|J^{(t+1)} - J^{(t)}| / |J^{(t)}| < \varepsilon_J$, declare convergence for the current $\mu$.

**Trust-region outer loop** (if enabled):
6. If $\mu > \varepsilon_\mu$: set $\mu \leftarrow \mu / 2$, warm-start from current iterates, return to step 2.
7. If $\mu \leq \varepsilon_\mu$: set $\mu = 0$, run steps 2–5 one final time.

**When $H = I$:** Take $R \to 0$ (or equivalently fix $x_l(k) = y_l(k)$). The state step is trivial, the initialisation is unnecessary, and the algorithm reduces to a single COSMIC solve — the original algorithm.

## 8. Computational Complexity

- **Initialisation:** Single linear solve, block-structured, $O(N\,(n + n\,q)^2)$ per trajectory.
- **State step:** RTS smoother, $O(N\,n^3)$ per trajectory, $O(L\,N\,n^3)$ total.
- **COSMIC step:** Standard COSMIC tridiagonal solve, $O(N\,(n+q)^3)$, independent of $L$ (trajectories are pooled into the data matrices).
- **Per outer iteration:** $O(L\,N\,n^3 + N\,(n+q)^3)$.
- **Total:** Proportional to the number of inner iterations times the number of trust-region steps. In practice, convergence is expected within a small number of alternating iterations (typically 5–20), and trust-region steps (typically 0–5 if enabled).

The linear scaling in $N$ — the hallmark of COSMIC — is preserved.

## 9. Hyperparameter Selection

**$\lambda$ (dynamics smoothness):** Same role and selection criteria as in standard COSMIC. Controls the trade-off between data fidelity and temporal smoothness of the estimated system matrices.

**$R$ (measurement noise covariance):** The measurement noise covariance matrix, weighting the observation fidelity term via $R^{-1}$. When $R$ is known from sensor specifications or calibration, it is used directly — no tuning is required for the observation term. When $R$ is unknown, setting $R = I$ recovers unweighted least squares. The relative scaling between $R^{-1}$ and the dynamics fidelity term (which implicitly assumes unit process noise covariance) determines the balance between trusting measurements and trusting the dynamics model.

**$\mu$ (trust-region):** Start at $\mu = 1$ if used, halve adaptively. For well-conditioned problems, $\mu = 0$ from iteration 2 onward (i.e., trust-region disabled) is expected to suffice.

## 10. Relationship to Existing Methods

**Expectation-Maximisation (EM) for LTV systems.** The alternating scheme is structurally identical to EM for linear dynamical systems: the state step is the E-step (state smoothing given parameters), the COSMIC step is the M-step (parameter estimation given states). The key difference is that the M-step uses COSMIC's smoothness regularisation rather than maximum likelihood, making this a MAP-EM algorithm.

**Subspace identification (N4SID, MOESP).** These methods can provide alternative initialisations by estimating a local LTI model within sliding windows. The $A = I$ initialisation proposed here avoids the stationarity assumption inherent in windowed subspace methods and is better suited for systems with continuous parameter variation.

**Kalman smoother.** The state step is exactly an RTS smoother. Output-COSMIC can be viewed as wrapping a Kalman smoother with a COSMIC parameter estimator in an EM loop, with the COSMIC smoothness prior replacing the standard ML parameter update.
