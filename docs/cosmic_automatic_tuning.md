# Automatic Tuning for Discrete LTV Regularised System Identification

## 1. Problem Statement

### 1.1 The Identification Problem

Consider the discrete-time linear time-varying (LTV) system

$$x(k+1) = A(k)x(k) + B(k)u(k), \quad k = 0, \ldots, N-1,$$

where $x(k) \in \mathbb{R}^p$ is the state, $u(k) \in \mathbb{R}^q$ is the control input, and the system matrices $A(k) \in \mathbb{R}^{p \times p}$, $B(k) \in \mathbb{R}^{p \times q}$ vary with time.

The COSMIC algorithm [1] identifies these matrices from $L$ trajectories of state-input data by solving a regularised least-squares problem. Defining $\mathbf{C}(k) = [A^T(k) \; B^T(k)]^T$ and the data matrix $\mathbf{D}(k) = [X^T(k) \; U^T(k)]$, the problem is

$$\min_{\mathbf{C}} \; f(\mathbf{C}) := \underbrace{\frac{1}{2} \sum_{k=0}^{N-1} \|\mathbf{D}(k)\mathbf{C}(k) - \mathbf{X}'^T(k)\|_F^2}_{h(\mathbf{C}) \;\text{(data fidelity)}} + \underbrace{\frac{1}{2} \sum_{k=1}^{N-1} \lambda_k \|\mathbf{C}(k) - \mathbf{C}(k-1)\|_F^2}_{g(\mathbf{C}) \;\text{(temporal smoothness)}}.$$

The regularisation parameters $\lambda_k > 0$ penalise variation in system dynamics between consecutive time steps. COSMIC admits a closed-form solution via LU factorisation of a block tridiagonal system, with computational cost linear in $N$ and cubic in $p + q$ [1].

### 1.2 The Tuning Problem

The quality of the identified model depends critically on the choice of $\lambda_k$. This is a standard bias-variance trade-off:

- **$\lambda_k$ too large:** The solution is over-smoothed — it cannot track genuine dynamics changes and the identified model is biased toward time-invariance.
- **$\lambda_k$ too small:** The solution overfits measurement noise, producing spurious variation in $A(k)$ and $B(k)$.

In existing work, $\lambda_k$ is either set to a single constant $\lambda$ across all time steps [1, §3.2], or split into a small number of manually chosen regimes (e.g., $\lambda_\text{middle}$ and $\lambda_\text{borders}$ in the Comet Interceptor application [1, §5]). Tuning is performed via grid search over orders of magnitude, evaluating each candidate against a validation dataset with noise-free state measurements [2].

This document presents three increasingly sophisticated approaches to automating this process.


## 2. Method 1: L-Curve Selection (Scalar $\lambda$, No Validation Data)

### 2.1 Rationale

For a constant $\lambda_k = \lambda$, the COSMIC solution $\mathbf{C}^*(\lambda)$ traces a one-parameter family as $\lambda$ varies from $0$ to $\infty$. At each $\lambda$, the solution achieves specific values of the two competing objectives:

$$R(\lambda) = h(\mathbf{C}^*(\lambda)) = \frac{1}{2} \sum_k \|\mathbf{D}(k)\mathbf{C}^*(k) - \mathbf{X}'^T(k)\|_F^2 \quad \text{(residual)},$$

$$S(\lambda) = \frac{1}{2} \sum_{k=1}^{N-1} \|\mathbf{C}^*(k) - \mathbf{C}^*(k-1)\|_F^2 \quad \text{(unweighted smoothness penalty)}.$$

Note that $S(\lambda)$ is the *unweighted* variation — the total squared change in system matrices — not the regularisation term $g(\mathbf{C})$ which includes the $\lambda$ factor. This is essential: plotting the regularisation-weighted term against the residual does not produce a useful L-curve because $\lambda$ appears in both the solution and the penalty.

The parametric curve $(\log S(\lambda), \log R(\lambda))$ as $\lambda$ traverses $(0, \infty)$ has a characteristic L-shape:

- **Small $\lambda$** (bottom-right): Low residual $R$, high variation $S$. The model fits the data closely but allows large inter-step changes.
- **Large $\lambda$** (top-left): High residual $R$, low variation $S$. The model is nearly time-invariant.
- **Corner** (the "knee"): The point of maximum curvature, where increasing $\lambda$ begins to degrade data fidelity without substantial further smoothing. This represents the best unsupervised compromise.

### 2.2 Efficient Corner Search via Bisection

Because COSMIC is a closed-form solver, evaluating $R(\lambda)$ and $S(\lambda)$ for a given $\lambda$ costs $O(N(p+q)^3)$ — the same as a single identification run. This makes repeated evaluation feasible.

Rather than gridding $\lambda$ over a predefined range, the corner can be located efficiently. The curvature $\kappa(\lambda)$ of the L-curve in the $(\log S, \log R)$ plane is unimodal in well-behaved regularisation problems — it has a single peak at the corner. This allows a golden section search or bisection on the curvature:

**Algorithm (L-Curve Corner Search):**

1. Set $\lambda_{\min}$, $\lambda_{\max}$ to bracket the expected range (e.g., $10^{-3}$ to $10^{10}$). Work in log-space: $\ell = \log_{10} \lambda$.
2. Evaluate $R(\lambda)$ and $S(\lambda)$ at three initial points $\ell_1 < \ell_2 < \ell_3$ (e.g., endpoints and midpoint).
3. Approximate the curvature $\kappa$ at $\ell_2$ using finite differences on the parametric curve $(\log S(\ell), \log R(\ell))$.
4. Trisect or golden-section the interval to find the $\ell$ that maximises $\kappa$.
5. Terminate when $|\ell_{\text{high}} - \ell_{\text{low}}| < \varepsilon$ (e.g., $\varepsilon = 0.5$, i.e., half an order of magnitude).

The curvature of a parametric curve $(\xi(\ell), \rho(\ell))$ with $\xi = \log S$, $\rho = \log R$ is

$$\kappa(\ell) = \frac{|\xi' \rho'' - \rho' \xi''|}{(\xi'^2 + \rho'^2)^{3/2}},$$

where derivatives are with respect to $\ell$ and can be approximated by finite differences over three evaluations. A golden-section search over an initial bracket spanning 13 orders of magnitude converges in roughly 20 COSMIC evaluations.

### 2.3 Properties

**Strengths:**

- Requires no validation data — only the training data used for identification.
- Fully automatic given a bracket $[\lambda_{\min}, \lambda_{\max}]$.
- Computationally modest: roughly 20 COSMIC runs, each with linear complexity.
- Well-studied in the inverse problems literature [3, 4], with known theoretical properties for Tikhonov-type regularisation.

**Limitations:**

- Restricted to a single scalar $\lambda$ across the entire trajectory. For systems whose rate of change varies over time, a constant $\lambda$ is suboptimal — it will under-smooth during rapid transitions and over-smooth during quiescent periods.
- The L-curve corner is a heuristic: it does not directly minimise any prediction-relevant criterion. In some problems the corner is poorly defined (e.g., when the curve is nearly straight).
- Assumes the L-curve is well-behaved (unimodal curvature). Pathological data or model mismatch can produce multi-modal curvature profiles, causing the search to converge to a local maximum.


## 3. Method 2: Validation-Based Tuning

### 3.1 Rationale

When a separate validation dataset is available — containing trajectories with noise-free (or low-noise) state measurements — the regularisation parameter can be selected to minimise a prediction-relevant loss function directly.

### 3.2 The Trajectory Prediction Loss

Following [2], define the trajectory prediction loss:

$$L(\lambda) = \frac{1}{|S|} \sum_{l \in S} \sqrt{\frac{1}{N} \sum_{k=1}^{N} \sum_{m=1}^{p} \left(\hat{x}_{k,m}^{(l)}(\lambda) - x_{k,m}^{(l)}\right)^2},$$

where $\hat{x}^{(l)}(\lambda)$ is the state trajectory predicted by the model identified with regularisation $\lambda$, propagated from the initial condition $x^{(l)}(0)$ of validation trajectory $l$, and $x^{(l)}$ is the true state.

The loss is evaluated by:

1. Run COSMIC with a candidate $\lambda$ to obtain $\{\hat{A}(k), \hat{B}(k)\}$.
2. For each validation trajectory $l$, propagate: $\hat{x}^{(l)}(k+1) = \hat{A}(k)\hat{x}^{(l)}(k) + \hat{B}(k)u^{(l)}(k)$.
3. Compare with the true trajectory $x^{(l)}(k)$.

### 3.3 Efficient Search

Since $L(\lambda)$ involves propagation through $N$ time steps per trajectory, small errors in $A(k)$ and $B(k)$ accumulate. The loss is typically a smooth, unimodal function of $\log \lambda$, making it amenable to the same golden-section search described in Section 2.2.

For a scalar $\lambda$, the search cost is again roughly 20 COSMIC evaluations, each followed by $|S|$ trajectory propagations. The propagation cost is $O(|S| \cdot N \cdot p^2)$ per evaluation, which is negligible compared to the COSMIC solve for most practical system sizes.

### 3.4 Extension to Piecewise-Constant $\lambda_k$

When prior knowledge suggests the system has distinct regimes (e.g., slow dynamics far from closest approach, fast dynamics near closest approach [1, §5]), the trajectory can be partitioned into $R$ regions with separate regularisation constants $\lambda^{(1)}, \ldots, \lambda^{(R)}$. The validation loss becomes $L(\lambda^{(1)}, \ldots, \lambda^{(R)})$, and minimisation is performed over an $R$-dimensional space.

For $R = 2$ (as in the Comet Interceptor application), grid search over a $20 \times 20$ log-spaced grid requires 400 COSMIC runs — still feasible given linear complexity, but significantly more expensive than the scalar case. For $R \geq 3$, grid search becomes impractical and gradient-free optimisers (Nelder-Mead, Bayesian optimisation) are preferable.

### 3.5 Properties

**Strengths:**

- Directly optimises a criterion relevant to the end use (state prediction accuracy).
- Can tune piecewise-constant $\lambda_k$ schedules with a small number of regions.
- Theoretically sound: minimising prediction error on held-out data is the standard approach to hyperparameter selection.

**Limitations:**

- Requires validation trajectories with clean state measurements — a strong requirement in many practical scenarios.
- The partition into regions (and the number of regions $R$) must be specified manually, typically based on physical knowledge of the system.
- Does not scale to a fully time-varying $\lambda_k$ schedule (one value per time step): the search space is $N$-dimensional and no amount of grid search or gradient-free optimisation can navigate it.
- The validation loss is not always unimodal in $\lambda$, particularly when the system exhibits qualitatively different dynamics in different trajectory phases. Piecewise-constant schedules mitigate this but require knowing the partition.


## 4. Method 3: Spectral Regularisation (Time-Varying $\lambda_k$)

### 4.1 Motivation

Neither the L-curve nor validation-based tuning produces a fully time-varying $\lambda_k$ schedule that adapts to local dynamics. The L-curve is restricted to scalar $\lambda$; validation-based tuning can handle piecewise-constant schedules but requires both clean data and manual partitioning.

The key observation is that the rate of system variation — the quantity that $\lambda_k$ should track — is directly observable in the frequency domain. By running a non-parametric spectral analysis on sliding windows of the input-output data, we can estimate *where* and *how fast* the system's dynamics are changing, and use this to construct $\lambda_k$ automatically.

### 4.2 Sliding-Window Transfer Function Estimation

Given input $u(t)$ and output $y(t)$ (see Section 4.3 for what "output" means when COSMIC uses state measurements), divide the data into $K$ overlapping segments of length $L_\text{seg}$ with overlap $P$, centered at times

$$t_k = \left((k-1)(L_\text{seg} - P) + \frac{L_\text{seg}}{2}\right) T_s, \quad k = 1, \ldots, K.$$

For each segment, compute the H1 transfer function estimate

$$\hat{G}_k(\omega) = \frac{\hat{\Phi}_{yu}^{(k)}(\omega)}{\hat{\Phi}_{uu}^{(k)}(\omega)}$$

and the squared coherence

$$\hat{\gamma}_k^2(\omega) = \frac{|\hat{\Phi}_{yu}^{(k)}(\omega)|^2}{\hat{\Phi}_{yy}^{(k)}(\omega) \, \hat{\Phi}_{uu}^{(k)}(\omega)},$$

using the Welch averaged periodogram method within each segment. Any standard non-parametric spectral estimator can be substituted (e.g., Blackman-Tukey).

### 4.3 Choosing the Output Signal

COSMIC operates on state measurements $x(k)$ and inputs $u(k)$, while the spectral estimator requires an input-output pair. These are different objects, and the connection must be made explicitly.

- If the system has a natural output equation $y = Cx$, use $y$ and $u$.
- If only state measurements are available (the common case for COSMIC), treat each state component $x_i(k)$ as an output: compute $\hat{G}_i(\omega) = \hat{\Phi}_{x_i u}(\omega) / \hat{\Phi}_{uu}(\omega)$ for each $i = 1, \ldots, p$, then aggregate the variation across channels.
- For MIMO inputs, the transfer function matrix is estimated column by column, requiring sufficiently uncorrelated inputs.

### 4.4 Coherence-Weighted Variation Metric

Define the spectral variation between consecutive segments:

$$\Delta_k = \left(\frac{\sum_{i=1}^{n_f} w_k(\omega_i) \, |\hat{G}_{k+1}(\omega_i) - \hat{G}_k(\omega_i)|^2}{\sum_{i=1}^{n_f} w_k(\omega_i)}\right)^{1/2},$$

with the coherence-based weight

$$w_k(\omega_i) = \hat{\gamma}_k^2(\omega_i) \cdot \hat{\gamma}_{k+1}^2(\omega_i).$$

The product of coherences ensures that only frequency bins where *both* adjacent segments have reliable estimates contribute to $\Delta_k$. This provides two critical properties:

1. **Noise rejection.** Frequency bins dominated by noise have low coherence and are down-weighted, preventing the metric from mistaking estimation noise for genuine system change.
2. **Excitation disambiguation.** A frequency band with poor input excitation produces low coherence regardless of whether the system is changing at that frequency. The weighting correctly prevents such bins from generating false positives. Conversely, it prevents them from generating false negatives in the other direction: if coherence is low, the spectral estimate is unreliable, and we should not draw conclusions either way.

For systems where the transfer function magnitude varies greatly across frequency, a normalised variant may be preferable:

$$\Delta_k^{\text{norm}} = \left(\frac{\sum_{i=1}^{n_f} w_k(\omega_i) \left|\frac{\hat{G}_{k+1}(\omega_i) - \hat{G}_k(\omega_i)}{\hat{G}_k(\omega_i)}\right|^2}{\sum_{i=1}^{n_f} w_k(\omega_i)}\right)^{1/2}.$$

### 4.5 Mapping $\Delta_k$ to $\lambda_k$

The variation sequence $\{\Delta_k\}$ must be mapped to a regularisation schedule $\{\lambda_k\}$ satisfying: $\Delta_k$ large $\Rightarrow$ $\lambda_k$ small (allow variation), $\Delta_k$ small $\Rightarrow$ $\lambda_k$ large (enforce smoothness). A sigmoid mapping provides this:

$$\lambda(\Delta) = \lambda_{\min} + \frac{\lambda_{\max} - \lambda_{\min}}{1 + (\Delta / \Delta_{\text{ref}})^\alpha},$$

where $\lambda_{\min}$, $\lambda_{\max}$ define the admissible range, $\Delta_{\text{ref}}$ sets the transition midpoint, and $\alpha > 0$ controls the steepness (default: $\alpha = 2$).

**Setting $\Delta_{\text{ref}}$.** A data-driven default is the median of $\{\Delta_1, \ldots, \Delta_{K-1}\}$. This ensures roughly half the segments receive above-midpoint regularisation and half below. See Section 4.8 for failure modes of this heuristic.

**Setting $\lambda_{\min}$ and $\lambda_{\max}$.** These can be initialised from the L-curve (Section 2): use the L-curve corner as a baseline $\lambda^*$, then set $\lambda_{\min} = \lambda^* / \beta$ and $\lambda_{\max} = \lambda^* \cdot \beta$ for some spread factor $\beta$ (e.g., $\beta = 10^3$). This anchors the spectral schedule to a data-driven scale rather than requiring the user to guess the order of magnitude.

### 4.6 Temporal Interpolation

The spectral segments (typically $K = 10$–$20$) are much coarser than the COSMIC time grid ($N$ samples). The schedule $\lambda_k^\text{spec}$ at the segment centers must be interpolated onto the COSMIC grid:

$$\lambda(t_j) = \text{interp}\left(\{t_k^\text{spec}\}, \{\lambda_k^\text{spec}\}, t_j\right),$$

using linear or shape-preserving piecewise cubic interpolation, with constant extrapolation at the boundaries.

### 4.7 The Relationship Between $\Delta G$ and $\Delta C$

A caveat must be noted. The spectral variation $\Delta_k$ measures changes in the input-output transfer function, while COSMIC's regularisation penalises changes in the state-space matrices $\mathbf{C}(k)$. The mapping between the two is nonlinear: $G(z) = C_o(zI - A)^{-1}B$ (with output matrix $C_o$), so the sensitivity of $G$ to perturbations in $A$ and $B$ depends on the current operating point and frequency.

This means $\Delta_k$ is a *qualitative indicator* of system variation, not a quantitative proxy for $\|\mathbf{C}(k) - \mathbf{C}(k-1)\|_F$. The method reliably detects *where* the system is changing and whether it is changing *fast or slowly*, but the absolute magnitude of $\Delta_k$ does not map linearly to the magnitude of state-space change. The sigmoid mapping in Section 4.5 accommodates this: it only needs $\Delta_k$ to be monotonically related to the true variation rate, not proportional to it.

### 4.8 Failure Modes and Robustness

**Uniform variation.** If all $\Delta_k$ are similar (the system changes at a roughly constant rate), the median heuristic places $\Delta_\text{ref}$ at the common value and the sigmoid produces roughly uniform $\lambda_k \approx (\lambda_{\min} + \lambda_{\max})/2$. This degenerates gracefully to a near-constant schedule, in which case the L-curve (Section 2) is equally appropriate.

**Sparse reconfiguration.** If the system is LTI for most of the trajectory with a brief reconfiguration event, the median reflects the LTI regime and $\Delta_\text{ref}$ is small. The reconfiguration segments have $\Delta_k \gg \Delta_\text{ref}$, correctly driving their $\lambda_k$ toward $\lambda_{\min}$. This is the ideal operating regime for the method.

**Noisy data with moderate coherence.** If coherence is uniformly moderate ($\gamma^2 \approx 0.5$), noise contributes to $\Delta_k$ and may cause systematic under-regularisation. In this regime, the coherence gating attenuates but does not eliminate the noise contribution.

**Unobservable mode changes.** Changes in $A(k)$ that affect only unobservable or weakly observable modes will not manifest in the transfer function and cannot be detected by the spectral pre-scan. The regularisation for these time steps will default to the ambient level, which may be too high.

### 4.9 Multiple Trajectories

COSMIC can exploit $L$ trajectories; the spectral pre-scan operates on one. When multiple trajectories are available:

- **Averaged variation.** Compute $\Delta_k^{(l)}$ for each trajectory, then average: $\bar{\Delta}_k = \frac{1}{L}\sum_l \Delta_k^{(l)}$. This reduces noise in the variation estimate. It assumes all trajectories experience the same time-varying dynamics, which is true by construction for a deterministic LTV system.
- **Conservative (minimum).** Use $\Delta_k = \min_l \Delta_k^{(l)}$, biasing toward higher regularisation. Appropriate when some trajectories have poor excitation in certain intervals.

### 4.10 Properties

**Strengths:**

- Produces a fully time-varying $\lambda_k$ schedule from training data alone — no validation data required.
- Captures the temporal structure of system variation without manual partitioning.
- Computationally cheap: the spectral pre-scan costs $O(K \cdot L_\text{seg} \log L_\text{seg})$, negligible relative to COSMIC.
- The coherence provides a built-in reliability indicator that prevents noise-driven false positives.

**Limitations:**

- Still requires setting $\lambda_{\min}$, $\lambda_{\max}$, $\alpha$, and $L_\text{seg}$ — four free parameters. The L-curve can anchor $\lambda_{\min}$ and $\lambda_{\max}$, but $\alpha$ and $L_\text{seg}$ require defaults or domain knowledge.
- The relationship between spectral variation and state-space variation is qualitative, not quantitative (Section 4.7).
- Cannot detect changes in unobservable modes.
- Inherits the time-frequency trade-off: the segment length limits both the temporal resolution of the $\lambda_k$ schedule and the frequency resolution of the spectral estimates.
- The method is unvalidated numerically. Its performance relative to scalar L-curve or validation-based tuning has not been demonstrated on benchmark systems.


## 5. Combined Pipeline

The three methods are not mutually exclusive. They can be composed into a pipeline where each stage refines the output of the previous one:

**Stage 1 — L-curve (Section 2).** Find the scalar $\lambda^*$ that best trades off data fidelity against smoothness. This requires only training data and roughly 20 COSMIC evaluations. The result anchors the scale of regularisation.

**Stage 2 — Spectral pre-scan (Section 4).** Using $\lambda^*$ to set $\lambda_{\min} = \lambda^* / \beta$ and $\lambda_{\max} = \lambda^* \cdot \beta$, construct a time-varying $\lambda_k$ schedule from the coherence-weighted spectral variation. This requires one pass over the data to compute the sliding-window spectra and one COSMIC evaluation with the resulting schedule.

**Stage 3 — Validation refinement (Section 3).** If validation data is available, refine $(\lambda_{\min}, \lambda_{\max}, \alpha)$ by minimising the trajectory prediction loss — a 3-dimensional search rather than the $N$-dimensional or scalar search of the standalone methods. Alternatively, apply the trajectory prediction loss to select among a small family of schedules (e.g., different $\alpha$ or $L_\text{seg}$ values).

The pipeline degrades gracefully depending on available resources:

- **Training data only:** Stages 1 + 2 produce a time-varying schedule with no manual tuning beyond bracket selection.
- **Training + validation data:** All three stages produce a refined schedule that combines spectral structure with prediction-optimal calibration.
- **Training data, LTI system suspected:** Stage 1 alone gives the scalar $\lambda^*$.


## 6. Summary

| | L-Curve | Validation | Spectral |
|---|---|---|---|
| **$\lambda_k$ type** | Scalar | Scalar or piecewise | Fully time-varying |
| **Validation data** | No | Yes (clean) | No |
| **Prior knowledge** | Bracket only | Regime partition | Bracket + segment length |
| **COSMIC evaluations** | ~20 | ~20 (scalar), ~400 ($R=2$) | 1 (+ pre-scan) |
| **Captures temporal structure** | No | Manually, with $R$ regions | Automatically |
| **Free parameters** | 2 ($\lambda_{\min}$, $\lambda_{\max}$) | $R$ + partition | 4 ($\lambda_{\min}$, $\lambda_{\max}$, $\alpha$, $L_\text{seg}$) |
| **Theoretical grounding** | Classical [3, 4] | Cross-validation | Novel (unvalidated) |

The three methods form a hierarchy of increasing sophistication: the L-curve provides a robust scalar baseline, validation-based tuning optimises a prediction-relevant criterion, and spectral regularisation introduces time-varying structure without requiring clean validation data or manual partitioning. The combined pipeline (Section 5) uses each method to initialise or constrain the next, reducing the total number of free parameters at each stage.

The spectral regularisation method remains a conjecture: it is conceptually motivated but has not been validated numerically. The next step is to test the full pipeline on the spring-mass-damper benchmarks in [1, 2], comparing the spectral-derived $\lambda_k$ schedule against scalar L-curve, grid-searched scalar $\lambda$, and grid-searched piecewise $\lambda_k$.


## References

1. M. Carvalho, C. Soares, P. Lourenço, R. Ventura. "COSMIC: fast closed-form identification from large-scale data for LTV systems." arXiv:2112.04355v2, 2022.

2. P. Łaszkiewicz, M. Carvalho, C. Soares, P. Lourenço. "The impact of modeling approaches on controlling safety-critical, highly perturbed systems: the case for data-driven models." arXiv:2509.13531v1, 2025.

3. P. C. Hansen. "The L-curve and its use in the numerical treatment of inverse problems." In *Computational Inverse Problems in Electrocardiology*, WIT Press, pp. 119–142, 2001.

4. P. C. Hansen, D. P. O'Leary. "The use of the L-curve in the regularization of discrete ill-posed problems." SIAM J. Sci. Comput., 14(6):1487–1503, 1993.

5. L. Ljung. *System Identification: Theory for the User*, 2nd ed. Prentice Hall, 1999.

6. P. Welch. "The use of fast Fourier transform for the estimation of power spectra." IEEE Trans. Audio Electroacoust., 15(2):70–73, 1967.
