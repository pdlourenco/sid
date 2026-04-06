# COSMIC (Closed-form Optimal data-driven linear time-varying SysteM IdentifiCation) as a Kalman Smoother: Online LTV Identification

**Purpose:** Show that the COSMIC algorithm, under the Bayesian interpretation
developed in the companion derivation, is exactly a Rauch-Tung-Striebel (RTS)
smoother in parameter space. This equivalence yields a natural online/recursive
formulation and connects λ selection to Kalman filter consistency diagnostics.

**Prerequisites:** The [Bayesian interpretation document](uncertainty_derivation.md),
which establishes the probabilistic model and the posterior covariance Σ ⊗ P(k).

---

## 1. The Parameter-Space State-Space Model

The Bayesian model from the companion derivation defines two equations:

**Process model (random walk prior on system matrices):**

```
C(k+1) = C(k) + η(k),     η(k) ~ MN(0, (1/λ_{k+1}) I_{p+q}, Σ)
```

This says the system matrices evolve as a random walk, with step size controlled
by 1/λ. Larger λ → smaller steps → smoother dynamics.

**Measurement model (data likelihood):**

```
X'(k)ᵀ = D(k) C(k) + W(k),     rows of W(k) i.i.d. ~ N(0, Σ)
```

This says the observed next-states are linear in C(k), corrupted by noise.

Together, these form a linear Gaussian state-space model where:
- The "state" is C(k) ∈ ℝ^{(p+q)×p} (the system matrices at time k)
- The "state transition" is the identity plus noise (random walk)
- The "measurement" at time k is the batch of L trajectory observations

This is a standard Kalman filtering/smoothing setup, with one important feature:
the Kronecker structure Σ ⊗ (·) means the p columns of C(k) share identical
dynamics and are coupled only through the measurement noise covariance Σ. As
shown in the companion derivation, Σ cancels from the normal equations, so
all filtering/smoothing operations can be performed on the "row covariance"
P(k) ∈ ℝ^{(p+q)×(p+q)} alone, without knowing Σ.


## 2. The Kalman Filter (Online/Forward Pass)

### 2.1 Algorithm

Starting from a diffuse (uninformative) prior on C(0), the Kalman filter
processes data sequentially:

**Initialization:**

```
Ĉ_filt(-1) = 0          (or any initial guess)
P_filt(-1) = αI          (large α → diffuse prior; α → ∞ recovers improper flat prior)
```

**For k = 0, 1, ..., N-1:**

*Predict (propagate through random walk):*

```
Ĉ_pred(k) = Ĉ_filt(k-1)
P_pred(k) = P_filt(k-1) + (1/λ_k) I
```

The covariance inflates by (1/λ_k)I: the system might have changed since
the last time step, so we are less certain.

(Convention: λ_0 is the prior uncertainty on the first step. If using a
diffuse prior, this is absorbed into the large initial P_filt(-1).)

*Update (incorporate data at time k):*

The data D(k), X'(k) provide a linear observation of C(k). In information
(precision) form, the update is additive:

```
P_filt(k)⁻¹ = P_pred(k)⁻¹ + D(k)ᵀ D(k)

Ĉ_filt(k) = P_filt(k) [P_pred(k)⁻¹ Ĉ_pred(k) + D(k)ᵀ X'(k)ᵀ]
```

Equivalently, in covariance form (standard Kalman gain):

```
K(k) = P_pred(k) D(k)ᵀ (D(k) P_pred(k) D(k)ᵀ + I_L)⁻¹

Ĉ_filt(k) = Ĉ_pred(k) + K(k) [X'(k)ᵀ - D(k) Ĉ_pred(k)]

P_filt(k) = (I - K(k) D(k)) P_pred(k)
```

**But the information form is preferred** because:
- D(k)ᵀD(k) is (p+q) × (p+q), already computed for COSMIC
- The covariance form requires inverting an L×L matrix (D P D' + I), which is
  expensive when L (number of trajectories) is large
- The information form inverts only (p+q) × (p+q) matrices

### 2.2 Complexity

Each time step requires:
- One (p+q) × (p+q) matrix inversion (for P_filt(k))
- One (p+q) × (p+q) matrix-matrix multiply (for the mean update)
- The Gram matrix D(k)ᵀD(k) and the right-hand side D(k)ᵀX'(k)ᵀ

**Total per step: O((p+q)³ + L(p+q)²)**

The L(p+q)² term is for assembling D(k)ᵀD(k) from the L trajectory data
vectors. Once assembled, the filter step itself is O((p+q)³), independent
of L. For the full sequence of N steps: O(N(p+q)³ + NL(p+q)²).

### 2.3 What the Filter Provides

At each time step k, the filter produces:
- Ĉ_filt(k): the best estimate of C(k) using only data from times 0,...,k
- P_filt(k): the posterior (row) covariance using only past data

This is a causal estimate — it does not use future data. It is therefore
suitable for real-time control, online monitoring, and streaming applications.

The filter estimate has **higher uncertainty** than the batch COSMIC estimate
at any interior time step (0 < k < N-1), because it lacks the information
from future data. At k = N-1 (the last step), the filter and the forward pass
of COSMIC produce identical results.


## 3. The RTS Smoother (Batch = COSMIC)

### 3.1 Algorithm

The Rauch-Tung-Striebel (RTS) smoother takes the forward filter output and
runs a backward pass to incorporate future data:

**Initialization:**

```
Ĉ_smooth(N-1) = Ĉ_filt(N-1)
P_smooth(N-1) = P_filt(N-1)
```

**For k = N-2, ..., 0:**

```
G(k) = P_filt(k) × P_pred(k+1)⁻¹

Ĉ_smooth(k) = Ĉ_filt(k) + G(k) [Ĉ_smooth(k+1) - Ĉ_pred(k+1)]

P_smooth(k) = P_filt(k) + G(k) [P_smooth(k+1) - P_pred(k+1)] G(k)ᵀ
```

The gain G(k) is a (p+q) × (p+q) matrix that blends the forward (filter)
estimate with the backward (smoother) correction.

### 3.2 Equivalence to COSMIC

**Claim:** The RTS smoother output (Ĉ_smooth(k), P_smooth(k)) is identical to
the batch COSMIC solution (Ĉ(k), P(k)) derived in the companion document.

**Proof sketch:** Both solve the same optimization problem — minimizing the
negative log-posterior of the Bayesian model. The block tridiagonal LU
factorization (COSMIC) and the Kalman smoother (RTS) are algebraically
equivalent decompositions of the same system of normal equations.

Specifically:
- COSMIC's forward pass computes Λ_k and Y_k. These are related to the
  Kalman filter quantities by:

  ```
  Λ_k = P_filt(k)⁻¹
  Y_k = P_filt(k)⁻¹ Ĉ_filt(k) = Λ_k Ĉ_filt(k)
  ```

  (In COSMIC notation, Y_k is the "information-weighted mean," not the
  mean itself. The mean is Ĉ_filt(k) = Λ_k⁻¹ Y_k.)

- COSMIC's backward pass computes C(k) from Y_k and C(k+1). This is
  algebraically identical to the RTS backward pass.

- COSMIC's P(k) recursion (from the companion derivation) gives:

  ```
  P(k) = (Λ_k - λ_{k+1}² P(k+1))⁻¹
  ```

  This is algebraically equivalent to the RTS covariance recursion.

The equivalence is exact, not approximate. COSMIC is the RTS smoother written
in information form with block tridiagonal algebra.

### 3.3 Mapping Between Notations

```
COSMIC                          Kalman (information form)
──────                          ────────────────────────
Λ_k                             P_filt(k)⁻¹
Y_k                             P_filt(k)⁻¹ Ĉ_filt(k)
C(k) (backward pass output)     Ĉ_smooth(k)
P(k) (diagonal of A⁻¹)          P_smooth(k)
S_{kk}                          P_pred(k)⁻¹ + D(k)ᵀD(k)
                                 = (P_filt(k-1) + (1/λ_k)I)⁻¹ + D(k)ᵀD(k)
λ_k                              1 / process_noise_variance_at_step_k
```


## 4. Verification: Deriving COSMIC from the Kalman Filter

To make the equivalence concrete, we establish the exact relationship between
COSMIC's forward Schur complements and the Kalman filter's information matrix.

### 4.1 The Offset Identity

**Claim.** The COSMIC forward Schur complement Λ_k and the Kalman filter
posterior precision P_filt(k)⁻¹ differ by a constant offset:

```
Λ_k = P_filt(k)⁻¹ + λ_{k+1} I       for k = 0, ..., N-2
Λ_{N-1} = P_filt(N-1)⁻¹              (no next-step term at the boundary)
```

This offset arises because COSMIC's block diagonal S_{kk} includes the
regularization coupling to step k+1 (the λ_{k+1} I term), while the Kalman
filter defers this contribution to the predict step at k+1. The smoother
outputs are identical because the offset cancels in the backward pass.

### 4.2 Proof by Induction

**Base case (k = 0).** With a diffuse prior (P_filt(-1) → ∞I), the Kalman
predict step gives P_pred(0)⁻¹ → 0. The information-form update adds the
data contribution:

```
P_filt(0)⁻¹ = P_pred(0)⁻¹ + D(0)ᵀD(0) = D(0)ᵀD(0)
```

COSMIC defines:

```
Λ_0 = S_{00} = D(0)ᵀD(0) + λ_1 I
```

Difference: Λ_0 - P_filt(0)⁻¹ = λ_1 I. ✓

**Inductive step.** Assume Λ_{k-1} = P_filt(k-1)⁻¹ + λ_k I for some k ≥ 1.

*Kalman predict-update.* The predicted precision at step k is:

```
P_pred(k)⁻¹ = (P_filt(k-1) + (1/λ_k) I)⁻¹
```

Applying the matrix inversion lemma with M = P_filt(k-1) and N = (1/λ_k)I:

```
P_pred(k)⁻¹ = λ_k I - λ_k² (P_filt(k-1)⁻¹ + λ_k I)⁻¹
```

By the inductive hypothesis, P_filt(k-1)⁻¹ + λ_k I = Λ_{k-1}, so:

```
P_pred(k)⁻¹ = λ_k I - λ_k² Λ_{k-1}⁻¹
```

After the data update:

```
P_filt(k)⁻¹ = P_pred(k)⁻¹ + D(k)ᵀD(k)
             = D(k)ᵀD(k) + λ_k I - λ_k² Λ_{k-1}⁻¹
```

*COSMIC recursion:*

```
Λ_k = S_{kk} - λ_k² Λ_{k-1}⁻¹
     = D(k)ᵀD(k) + (λ_k + λ_{k+1}) I - λ_k² Λ_{k-1}⁻¹
```

Taking the difference:

```
Λ_k - P_filt(k)⁻¹ = λ_{k+1} I     ✓
```

At the boundary k = N-1, S_{N-1,N-1} = D(N-1)ᵀD(N-1) + λ_{N-1} I (no
λ_N term), and the Kalman update gives P_filt(N-1)⁻¹ = D(N-1)ᵀD(N-1) +
λ_{N-1} I - λ_{N-1}² Λ_{N-2}⁻¹ = Λ_{N-1}. So the offset vanishes at the
last step.  ∎

### 4.3 Consequences for the Smoother

The offset λ_{k+1} I cancels exactly in the backward pass. To see this, note
that the COSMIC backward recursion:

```
C(k) = Y_k + λ_{k+1} Λ_k⁻¹ C(k+1)
```

and the RTS backward recursion:

```
Ĉ_smooth(k) = Ĉ_filt(k) + G(k) [Ĉ_smooth(k+1) - Ĉ_pred(k+1)]
```

both combine the forward estimate at k with a correction proportional to the
backward innovation at k+1. The smoother gain G(k) = P_filt(k) P_pred(k+1)⁻¹
absorbs the offset, so the smoother outputs (Ĉ_smooth(k), P_smooth(k)) are
algebraically identical to the COSMIC outputs (C(k), P(k)). The equivalence
is exact, not approximate.


## 5. Online Operation

### 5.1 Streaming Filter

For real-time applications, the Kalman filter (§2) runs without the backward
pass. At each new time step:

1. Receive new data: D(k), X'(k) from the L trajectories (or a subset)
2. Predict: inflate covariance by (1/λ_k)I
3. Update: incorporate data via information-form update
4. Output: Ĉ_filt(k), P_filt(k) — the current best estimate and uncertainty

**Cost per step:** O((p+q)³) after assembling D(k)ᵀD(k).

**Memory:** O((p+q)²) — only the current Ĉ_filt(k) and P_filt(k).

### 5.2 Relationship to Batch COSMIC

The streaming filter gives noisier estimates than batch COSMIC because it
lacks future data. The relationship is:

| Quantity | Filter (online) | Smoother (batch COSMIC) |
|----------|-----------------|------------------------|
| Mean at k | Ĉ_filt(k) | Ĉ_smooth(k) |
| Covariance at k | P_filt(k) | P_smooth(k) |
| Data used | 0,...,k | 0,...,N-1 |
| Causal? | Yes | No |

For slowly varying systems (large λ), the filter and smoother are close at
all interior points. For rapidly varying systems (small λ), the filter may
lag behind true variations while the smoother tracks them in hindsight.

### 5.3 Warm-Start Batch Updates

A middle ground between pure online and full batch: periodically run the
backward pass over a recent window to refine the estimates.

**Growing-horizon smoother:**
1. Run the filter forward as data arrives
2. When a batch of W new time steps has accumulated, run the RTS backward
   pass over the most recent W steps
3. The smoothed estimates for those W steps are now optimal (given data up
   to the current time)
4. Continue filtering forward from the smoothed state

**Sliding-window smoother:**
1. Maintain a window of the most recent W time steps
2. When new data arrives, drop the oldest step and add the new one
3. Run full COSMIC (forward + backward) on the W-step window
4. Output the smoothed estimate at the window center

Both approaches cost O(W(p+q)³) per update, with W controlling the tradeoff
between latency and estimation quality.


## 6. Innovation-Based Diagnostics

### 6.1 Innovations

The Kalman filter naturally produces **innovations** — the prediction error
before incorporating new data:

```
ε(k) = X'(k)ᵀ - D(k) Ĉ_pred(k) ∈ ℝ^{L × p}
```

This is the difference between the observed next-states and what the model
predicted. Under the assumed model (correct λ, Gaussian noise), the innovations
have known statistics:

```
rows of ε(k) are i.i.d. ~ N(0, Σ_innov(k))
```

where

```
Σ_innov(k) = D(k) P_pred(k) D(k)ᵀ + I_L       (L × L)
```

In practice, we don't need the full L×L matrix. The key diagnostic is the
normalized innovation squared:

```
NIS(k) = (1/Lp) tr(Σ̂⁻¹ ε(k)ᵀ ε(k))
```

Under the correct model, E[NIS(k)] ≈ 1 (plus a correction from P_pred(k) that
vanishes when L is large relative to p+q).

### 6.2 Using Innovations for λ Selection

**If NIS(k) is systematically > 1:** The model under-predicts the innovation
magnitude. This means either:
- λ is too large (the model is too smooth and can't track real system changes),
  so the prediction errors are larger than expected
- The noise model is wrong (Σ is underestimated)

**If NIS(k) is systematically < 1:** The model over-predicts the innovation
magnitude. This means either:
- λ is too small (the model is too flexible and is fitting noise, making
  overly cautious predictions)
- The noise model is wrong (Σ is overestimated)

**If NIS(k) ≈ 1 on average:** The model is consistent.

The optimal λ can be selected by:

```
λ* = argmin_λ | (1/N) Σ_k NIS(k) - 1 |
```

or, more robustly, by a chi-squared consistency test on the innovations.

### 6.3 Cross-Validation Against sidFreqMap

The innovations provide a time-domain diagnostic. For a frequency-domain
diagnostic, compare the frozen transfer function at each time step against
the non-parametric `sidFreqMap` estimate.

At each time step k and frequency ω:

```
G_cosmic(k, ω) = C_out (e^{jω}I - A_filt(k))⁻¹ B_filt(k)
G_BT(k, ω) = sidFreqMap estimate at time k
```

Under a consistent model:

```
|G_cosmic(k,ω) - G_BT(k,ω)|² ≤ Var_cosmic(k,ω) + Var_BT(k,ω)
```

at most frequency-time grid points (e.g., ≥ 90% of grid points at 95%
individual coverage). This gives a frequency-domain λ consistency criterion
that is completely independent of the innovation-based criterion.

**Combined λ selection:**
1. Use trajectory prediction loss (time-domain, held-out data) for initial λ
2. Verify against innovation consistency (time-domain, filter-based)
3. Verify against `sidFreqMap` consistency (frequency-domain, independent method)

If all three agree, confidence in λ is high. If they disagree, the
discrepancy identifies the nature of the model deficiency.


## 7. Variable-Length Trajectories in Online Mode

The online filter handles variable-length trajectories naturally. At each time
step k, D(k) is assembled from whichever trajectories L(k) have data at k.

Trajectories can:
- Start at different times (a new trajectory appears → D(k) gains a row)
- End at different times (a trajectory terminates → D(k) loses a row)
- Have gaps (a trajectory is missing at some time steps → excluded from D(k))

The filter equations are unchanged. The only difference is that D(k)ᵀD(k) and
D(k)ᵀX'(k)ᵀ are computed from the available trajectories at each k.

**Online trajectory addition:** When a new trajectory ℓ begins at time k_start:
- For k < k_start: trajectory ℓ is absent from D(k)
- For k ≥ k_start: trajectory ℓ contributes to D(k)

No reprocessing of past data is needed. The filter simply incorporates the new
trajectory's data as it arrives.


## 8. Implementation Plan for sid

### 8.1 New Functions

```
sidLtvCosmicOnline.m              % Streaming Kalman filter in parameter space
private/sidLtvCosmicPredict.m    % One predict step: P̄ = P + (1/λ)I
private/sidLtvCosmicUpdate.m     % One update step: incorporate D(k), X'(k)
private/sidLtvCosmicSmooth.m     % RTS backward pass over stored filter output
private/sidLtvCosmicInnovation.m % Innovation computation and NIS diagnostic
```

### 8.2 API

**Streaming mode:**

```matlab
% Initialize filter
state = sidLtvCosmicOnline('init', p, q, 'Lambda', 1e5);

% Process data as it arrives
for k = 1:N
    % Get data for this time step (variable number of trajectories)
    Dk = buildDataMatrix(trajectories, k);   % |L(k)| × (p+q)
    Xpk = buildNextState(trajectories, k);    % p × |L(k)|

    % One filter step
    [Ck, Pk, state] = sidLtvCosmicOnline('step', state, Dk, Xpk);

    % Ck is the current estimate of [A(k)'; B(k)']
    % Pk is the current (row) covariance
end

% Optional: smooth over all stored data
[C_smooth, P_smooth] = sidLtvCosmicOnline('smooth', state);
```

**Batch mode (equivalent to sidLtvCosmic):**

```matlab
state = sidLtvCosmicOnline('init', p, q, 'Lambda', 1e5);
for k = 1:N
    state = sidLtvCosmicOnline('step', state, Dk, Xpk);
end
[C_smooth, P_smooth] = sidLtvCosmicOnline('smooth', state);
% C_smooth is identical to sidLtvCosmic output
```

### 8.3 Storage Requirements

| Mode | Memory | Notes |
|------|--------|-------|
| Filter only | O((p+q)²) | Current state only; no smoothing possible |
| Filter + deferred smoothing | O(N(p+q)²) | Store all P_filt(k), Ĉ_filt(k) for backward pass |
| Sliding window | O(W(p+q)²) | Fixed memory, smoothing over window of size W |

### 8.4 Validation Tests

1. **Filter-smoother agreement at k=N-1:** The filter output at the last time
   step must equal the smoother (batch COSMIC) output. Numerical tolerance: eps.

2. **Smoother = batch COSMIC:** After running filter + RTS backward pass, the
   output must be identical (to numerical precision) to `sidLtvCosmic`.

3. **Innovation whiteness:** On synthetic data with known λ and Σ, the
   innovations ε(k) should be white (uncorrelated across k) and have the
   predicted covariance. Test via autocorrelation and NIS ≈ 1.

4. **λ recovery:** On synthetic data generated with a known λ, the
   innovation-based λ selection should recover the true value.

5. **Online = batch for LTI:** On data from an LTI system, the filter should
   converge to the batch estimate after sufficient data, with the filter
   covariance shrinking monotonically.


## 9. Connections and Implications

### 9.1 COSMIC is Not New (But the Packaging Is)

The equivalence to the RTS smoother means COSMIC is, at its core, a well-known
algorithm applied in a specific context. What the COSMIC paper contributes is:
- The block tridiagonal LU formulation, which is more efficient than the
  standard RTS implementation when processing batch data
- The application to LTV system identification (rather than state estimation)
- The existence/uniqueness conditions in terms of the empirical covariance

### 9.2 Extensions from the Kalman Literature

The Kalman filter equivalence immediately suggests extensions:

- **Square-root filtering:** For numerical stability, work with the Cholesky
  factor of P rather than P itself. This prevents loss of positive definiteness
  due to roundoff.

- **Adaptive λ:** Allow λ_k to vary over time, estimated from the innovations.
  This is the "adaptive Kalman filter" idea: if innovations are large, decrease
  λ (allow more variation); if small, increase λ (enforce smoothness). This
  automates the λ selection problem entirely.

- **Robust filtering:** Replace the Gaussian noise model with a heavy-tailed
  distribution (e.g., Student-t), giving robustness to outliers. This connects
  to the SBCD extension of COSMIC mentioned in the original paper.

- **Multiple-model filtering:** Run several filters with different λ values
  in parallel, weight by their innovation likelihoods. This is the
  Interacting Multiple Model (IMM) approach, which handles systems that switch
  between fast-varying and slow-varying regimes.

### 9.3 When to Use Filter vs. Smoother

| Scenario | Recommended | Why |
|----------|-------------|-----|
| Offline analysis | Smoother (batch COSMIC) | Uses all data; minimum variance |
| Real-time control | Filter (online) | Causal; O(1) memory and compute per step |
| Near-real-time monitoring | Sliding-window smoother | Low latency with some future data |
| Post-experiment analysis | Smoother | No time constraints; best estimates |
| Adaptive λ | Filter with NIS monitoring | Innovation diagnostics require filter |
