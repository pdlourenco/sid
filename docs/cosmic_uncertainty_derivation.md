# Bayesian Interpretation of COSMIC (Closed-form Optimal data-driven linear time-varying System IdentifiCation): Derivation

**Purpose:** Justify the posterior covariance of the COSMIC estimator under a
general (non-isotropic) noise model, and clarify the bias/variance tradeoffs.

---

## 1. Setup and Notation

The discrete linear time-varying (LTV) system with noise is:

```
x_ℓ(k+1) = A(k) x_ℓ(k) + B(k) u_ℓ(k) + w_ℓ(k)
```

for trajectory ℓ = 1,...,L and time step k = 0,...,N-1, where w_ℓ(k) ∈ ℝᵖ is
noise (measurement error, process disturbance, or both).

Define the parameter matrix at time k:

```
C(k) = [A(k)ᵀ; B(k)ᵀ] ∈ ℝ^{(p+q)×p}
```

and the data vector for trajectory ℓ at time k:

```
d_ℓ(k) = [x_ℓ(k); u_ℓ(k)] ∈ ℝ^{p+q}
```

The model becomes:

```
x_ℓ(k+1)ᵀ = d_ℓ(k)ᵀ C(k) + w_ℓ(k)ᵀ
```

Stacking all L trajectories at time k into D(k) ∈ ℝ^{L×(p+q)} and X'(k) ∈ ℝ^{p×L}:

```
X'(k)ᵀ = D(k) C(k) + W(k)       (*)
```

where row ℓ of W(k) is w_ℓ(k)ᵀ.


## 2. The Probabilistic Model

We state the two components of the Bayesian model explicitly.

### 2.1 Likelihood

Each noise vector is independent and Gaussian with a general p × p covariance:

```
w_ℓ(k) ~ N(0, Σ)     i.i.d. across ℓ and k
```

where Σ ∈ ℝ^{p×p} is symmetric positive definite. This allows different noise
levels per state component and cross-correlations between them. For example, in
a spring-mass-damper system, position might be measured with much less noise
than velocity, giving Σ = diag(σ_pos², σ_vel²) with σ_pos ≪ σ_vel.

The log-likelihood for all data, given C = {C(0),...,C(N-1)} and Σ, is:

```
log p(data | C, Σ) = -(1/2) Σ_{k=0}^{N-1} Σ_{ℓ=1}^{L}
                      (x_ℓ(k+1) - C(k)ᵀ d_ℓ(k))ᵀ Σ⁻¹ (x_ℓ(k+1) - C(k)ᵀ d_ℓ(k))
                      + const
```

where the constant absorbs -(NL/2) log|2πΣ|.

Using the trace identity ||M||²_{Σ⁻¹} = tr(Σ⁻¹ MᵀM), this becomes:

```
log p(data | C, Σ) = -(1/2) Σ_{k=0}^{N-1} tr(Σ⁻¹ E(k)ᵀ E(k)) + const
```

where E(k) = X'(k)ᵀ - D(k)C(k) ∈ ℝ^{L×p} is the residual matrix at time k.


### 2.2 Prior: Gaussian Random Walk on C(k)

The COSMIC regularization penalizes ||C(k) - C(k-1)||²_F. The natural Bayesian
interpretation that preserves the structure of the problem is a Gaussian random
walk prior whose covariance scales with the noise covariance:

```
vec(C(k) - C(k-1)) ~ N(0, (1/λ_k) Σ ⊗ I_{p+q})    for k = 1,...,N-1
```

Equivalently, the columns of the matrix increment are:

```
[C(k) - C(k-1)][:,j] ~ N(0, (Σ_{jj}/λ_k) I_{p+q})
```

and cross-column correlations follow Σ.

**Why scale the prior with Σ?** This choice means that the regularization
strength λ_k operates in noise-normalized units: a state component with 10×
more measurement noise is automatically allowed 10× more apparent variation in
its corresponding column of C(k). Without this scaling, the regularization
would over-smooth noisy components and under-smooth clean ones.

The prior log-density is:

```
log p(C | Σ) = -(1/2) Σ_{k=1}^{N-1} λ_k tr(Σ⁻¹ (C(k) - C(k-1))ᵀ (C(k) - C(k-1))) + const
```

**Note on C(0):** The random walk prior constrains only differences. There is no
prior on the absolute value of C(0). This is an improper (flat) prior on C(0),
which is acceptable: the posterior will be proper as long as the data at k=0
provides sufficient information, i.e., D(0)ᵀD(0) ≻ 0.


### 2.3 Factoring out Σ⁻¹: Why COSMIC Doesn't Need to Know Σ

The negative log-posterior is proportional to:

```
-log p(C | data, Σ) ∝ (1/2) Σ_k tr(Σ⁻¹ E(k)ᵀ E(k))
                     + (1/2) Σ_k λ_k tr(Σ⁻¹ ΔC(k)ᵀ ΔC(k))
```

where ΔC(k) = C(k) - C(k-1). Both terms have the form tr(Σ⁻¹ × ···), and the
factor Σ⁻¹ is common. Taking the gradient with respect to C(k):

```
∇_{C(k)} [-log p] = D(k)ᵀ E(k) Σ⁻¹ + λ_k ΔC(k) Σ⁻¹ - λ_{k+1} ΔC(k+1) Σ⁻¹
                   = [D(k)ᵀ E(k) + λ_k ΔC(k) - λ_{k+1} ΔC(k+1)] Σ⁻¹
```

Setting this to zero and right-multiplying by Σ (which is invertible):

```
D(k)ᵀ E(k) + λ_k ΔC(k) - λ_{k+1} ΔC(k+1) = 0
```

**The normal equations are identical to the isotropic (σ²I) case.** The matrix Σ
cancels entirely. Therefore:

- The COSMIC algorithm (forward/backward pass) is unchanged.
- The MAP estimate Ĉ does not depend on Σ.
- The implementation does not need to know or estimate Σ to compute Ĉ.

Σ enters only through the posterior covariance (§3) and must be estimated from
the residuals (§6).


## 3. Posterior Distribution

### 3.1 Matrix-Normal Posterior

Because Σ⁻¹ factors out of the entire problem, and the cost is quadratic in C,
the posterior distribution of C(k) is **matrix-normal**:

```
C(k) | data, Σ  ~  MN(Ĉ(k),  P(k),  Σ)
```

where:
- Ĉ(k) is the COSMIC solution (the MAP estimate / posterior mean)
- P(k) ∈ ℝ^{(p+q)×(p+q)} is the "row covariance" (same for all columns)
- Σ ∈ ℝ^{p×p} is the "column covariance" (the noise covariance)

The matrix-normal MN(M, U, V) means that vec(X) ~ N(vec(M), V ⊗ U). So:

```
Cov(vec(C(k))) = Σ ⊗ P(k)
```

This Kronecker structure has a clean interpretation:
- P(k) governs how much uncertainty there is in each row of C(k) — i.e.,
  uncertainty in the (p+q)-dimensional parameter vector — and is determined by
  the data geometry (how informative the data is at time k) and the
  regularization (how much neighboring time steps contribute).
- Σ governs how the p columns of C(k) are correlated — i.e., cross-correlation
  between state components — and is determined by the noise structure.

### 3.2 Covariance of Individual Entries

From the Kronecker structure, the covariance between any two entries of C(k) is:

```
Cov(C(k)_{ab}, C(k)_{cd}) = Σ_{bd} × P(k)_{ac}
```

where a,c index rows (1,...,p+q) and b,d index columns (1,...,p).

**Special cases:**
- Same column (b = d): Cov(C(k)_{ab}, C(k)_{cb}) = Σ_{bb} × P(k)_{ac}.
  This is just the b-th diagonal of Σ scaling the row covariance P(k).
- Diagonal Σ: columns become uncorrelated, and column j has
  Cov(c_j(k)) = Σ_{jj} × P(k).
- Isotropic Σ = σ²I: all columns have the same covariance σ² P(k),
  recovering the isotropic case.

### 3.3 Connection to COSMIC

The posterior mean is the COSMIC solution:

```
Ĉ_COSMIC = A⁻¹ Vᵀ X'ᵀ = posterior mean
```

This is identical to the MAP estimate. For a Gaussian posterior, the MAP and the
mean coincide, so there is no distinction.

The posterior covariance of the stacked parameter vector is:

```
Cov(vec(C)) = Σ ⊗ A⁻¹
```

and the per-timestep marginal is:

```
Cov(vec(C(k))) = Σ ⊗ P(k)     where P(k) = [A⁻¹]_{kk}
```

**This is not a frequentist statement.** It does not say that E[Ĉ] = C_true.
Rather, it says: given the data and the probabilistic model (Gaussian noise with
covariance Σ + Gaussian random walk prior scaled by Σ), the posterior belief
about C is matrix-normal with the above parameters.


### 3.4 Derivation of the Kronecker Structure

For completeness, we derive the posterior covariance from first principles.

Stack all time steps for column j into c_j = [c_j(0); ...; c_j(N-1)] ∈ ℝ^{N(p+q)},
where c_j(k) is column j of C(k). The full stacked parameter is
c = [c_1; ...; c_p] ∈ ℝ^{Np(p+q)}.

The stacked likelihood for all p columns at time k is:

```
X'(k)ᵀ = D(k) C(k) + W(k)
```

where W(k) has i.i.d. rows ~ N(0, Σ). Vectorizing:

```
vec(X'(k)ᵀ) = (Iₚ ⊗ D(k)) vec(C(k)) + vec(W(k))
```

with vec(W(k)) ~ N(0, Σ ⊗ I_L).

Stacking over all k:

```
x' = (Iₚ ⊗ V) c + w,     w ~ N(0, Σ ⊗ I_{NL})
```

The precision of the likelihood contribution is:

```
(Σ ⊗ I_{NL})⁻¹ = Σ⁻¹ ⊗ I_{NL}
```

so the likelihood precision in the c parameterization is:

```
(Iₚ ⊗ V)ᵀ (Σ⁻¹ ⊗ I) (Iₚ ⊗ V) = Σ⁻¹ ⊗ VᵀV
```

The prior on c has precision:

```
Σ⁻¹ ⊗ (FᵀΛF)
```

(from the prior covariance (1/λ_k) Σ ⊗ I on the differences).

The posterior precision is the sum:

```
Σ⁻¹ ⊗ VᵀV + Σ⁻¹ ⊗ FᵀΛF = Σ⁻¹ ⊗ (VᵀV + FᵀΛF) = Σ⁻¹ ⊗ A
```

Inverting:

```
Cov(c | data, Σ) = (Σ⁻¹ ⊗ A)⁻¹ = Σ ⊗ A⁻¹     ∎
```

This confirms the Kronecker structure and shows that P(k) = [A⁻¹]_{kk} is
independent of Σ. The matrix A depends only on the data geometry (DᵀD) and the
regularization (λ_k), not on the noise covariance.


## 4. Why This is Not the Frequentist Covariance

### 4.1 Frequentist Bias

If we instead ask: "For a fixed, deterministic C_true, what is E[Ĉ]?", the
MAP/COSMIC estimate satisfies:

```
E[Ĉ] = A⁻¹ VᵀV C_true = (VᵀV + FᵀΛF)⁻¹ VᵀV C_true
```

(The expectation is over the noise w, which enters only through X'. Since
E[X'ᵀ] = V C_true, the result follows from the linearity of the COSMIC solution
in X'.)

This equals C_true only when FᵀΛF = 0 (no regularization) or F C_true = 0
(the true system is LTI). Otherwise the estimator is biased:

```
Bias = E[Ĉ] - C_true = -A⁻¹ FᵀΛF C_true
```

The regularization pulls the estimate toward smoothness, introducing systematic
error wherever the true system genuinely varies.

### 4.2 Frequentist Variance (Sandwich Form)

The frequentist covariance of vec(Ĉ), conditional on a fixed C_true, is:

```
Cov_freq(vec(Ĉ)) = E[(vec(Ĉ) - E[vec(Ĉ)])(vec(Ĉ) - E[vec(Ĉ)])ᵀ]
```

Since Ĉ = A⁻¹ Vᵀ X'ᵀ and the randomness is in X' = V C_true + noise:

```
Ĉ - E[Ĉ] = A⁻¹ Vᵀ W_stacked
```

where W_stacked has covariance Σ ⊗ I_{NL} when vectorized. Therefore:

```
Cov_freq(vec(Ĉ)) = (Iₚ ⊗ A⁻¹)(Iₚ ⊗ Vᵀ)(Σ ⊗ I)(Iₚ ⊗ V)(Iₚ ⊗ A⁻¹)
                   = Σ ⊗ (A⁻¹ VᵀV A⁻¹)
```

This is the "sandwich" form. It differs from Σ ⊗ A⁻¹ because VᵀV ≠ A.
Specifically, in the positive semidefinite ordering:

```
Σ ⊗ (A⁻¹ VᵀV A⁻¹)  ≼  Σ ⊗ A⁻¹
```

with equality only when FᵀΛF = 0.

**Note:** The Kronecker structure is preserved in the frequentist case too — the
same Σ appears, only the "row" factor changes from A⁻¹ (Bayesian) to
A⁻¹ VᵀV A⁻¹ (frequentist).

**Interpretation:** The Bayesian posterior variance Σ ⊗ A⁻¹ is larger than the
frequentist variance. It accounts for both estimation uncertainty (from noise)
and prior uncertainty (from not knowing C a priori). The frequentist variance
captures only estimation uncertainty, while ignoring the bias entirely.

### 4.3 Frequentist MSE

The mean squared error combines both:

```
MSE = Cov_freq(vec(Ĉ)) + vec(Bias) vec(Bias)ᵀ
    = Σ ⊗ (A⁻¹ VᵀV A⁻¹) + vec(A⁻¹ FᵀΛF C_true) vec(...)ᵀ
```

This depends on the unknown C_true and therefore cannot be computed from data
alone.

### 4.4 Summary of Comparison

```
                          Accounts for   Accounts for    Computable
                          noise?         regularization? from data?
                          
Bayesian  Σ ⊗ A⁻¹         Yes            Yes (as prior)  Yes
Freq. sandwich             Yes            No (ignores     Yes
  Σ ⊗ A⁻¹VᵀVA⁻¹                         bias)
Freq. MSE                  Yes            Yes             No (needs
                                                          C_true)
```

The Bayesian posterior covariance Σ ⊗ A⁻¹ is the only option that is both
computable and accounts for the role of regularization. The tradeoff is that it
requires accepting the Bayesian framing.

### 4.5 Why the Sandwich is Impractical Anyway

Even if we wanted the frequentist sandwich, computing the diagonal blocks of
A⁻¹ VᵀV A⁻¹ requires the full k-th block row of A⁻¹:

```
[A⁻¹ VᵀV A⁻¹]_{kk} = Σ_j [A⁻¹]_{kj} (Dⱼᵀ Dⱼ) [A⁻¹]_{jk}
```

The inverse of a block tridiagonal matrix is generally full, so this sum runs
over all N blocks. Extracting the diagonal blocks costs O(N²(p+q)³), destroying
the linear-in-N complexity that makes COSMIC efficient.


## 5. Posterior Covariance: Block Structure

### 5.1 Structure of A

The matrix A = VᵀV + FᵀΛF is N(p+q) × N(p+q) with block tridiagonal structure.
Its (p+q) × (p+q) blocks are:

**Diagonal blocks:**

```
A_{00}     = D(0)ᵀD(0) + λ₁ I
A_{kk}     = D(k)ᵀD(k) + (λ_k + λ_{k+1}) I     for k = 1,...,N-2
A_{N-1,N-1} = D(N-1)ᵀD(N-1) + λ_{N-1} I
```

**Off-diagonal blocks:**

```
A_{k,k+1} = A_{k+1,k} = -λ_{k+1} I
```

**All other blocks are zero.**

Note that A does not depend on Σ. The noise covariance affects the posterior
only through the Kronecker product Σ ⊗ A⁻¹.

### 5.2 Diagonal Blocks of A⁻¹ via Backward Recursion

The forward pass of COSMIC computes the Schur complements from the left:

```
Λ₀ = A_{00}

Λ_k = A_{kk} - λ_k² Λ_{k-1}⁻¹     for k = 1,...,N-1
```

These are the pivots of the block LU factorization. Denote their inverses as
Λ_k⁻¹ (already computed and stored during the forward pass).

**Claim:** The diagonal blocks P(k) = [A⁻¹]_{kk} satisfy:

```
P(N-1) = Λ_{N-1}⁻¹

P(k) = (Λ_k - λ_{k+1}² P(k+1))⁻¹     for k = N-2,...,0
```

**Proof by induction:**

**Base case (k = N-1):** The matrix A restricted to the last block is
A_{N-1,N-1} = Λ_{N-1} (the forward Schur complement at the last step has no
further correction). So [A⁻¹]_{N-1,N-1} = Λ_{N-1}⁻¹. ✓

**Inductive step:** Assume P(k+1) = [A⁻¹]_{k+1,k+1} is correct. Partition
A at the boundary between blocks k and k+1:

```
A = | A_≤k    B    |     where B has a single nonzero block:
    | Bᵀ      A_≥k+1|        B_{k,k+1} = -λ_{k+1} I
```

By the block matrix inversion formula, the (k,k) block of A⁻¹ is the (k,k)
block of:

```
(A_≤k - B A_≥k+1⁻¹ Bᵀ)⁻¹
```

Since B has nonzero entries only in block-row k and block-column k+1:

```
[B A_≥k+1⁻¹ Bᵀ]_{kk} = λ_{k+1}² [A_≥k+1⁻¹]_{k+1,k+1} = λ_{k+1}² P(k+1)
```

The (k,k) element of A_≤k, after accounting for the Schur complement from
blocks 0,...,k-1, is exactly Λ_k (this is what the forward pass computes).
Therefore:

```
P(k) = (Λ_k - λ_{k+1}² P(k+1))⁻¹     ∎
```

### 5.3 Complexity

The backward recursion for P(k) requires:
- One (p+q) × (p+q) matrix subtraction per time step
- One (p+q) × (p+q) matrix inversion per time step

Total: O(N (p+q)³), identical to the forward pass. The Λ_k values from the
forward pass must be stored (N matrices of size (p+q) × (p+q)).

**Important:** The recursion computes P(k) = [A⁻¹]_{kk}, which is independent
of Σ. The full posterior covariance Σ ⊗ P(k) is formed only at the end, once Σ
has been estimated.


## 6. Estimating Σ

The noise covariance Σ is generally unknown and must be estimated from the
residuals.

### 6.1 Residual-Based Estimate

Given the COSMIC solution Ĉ, the residuals at time step k are:

```
E(k) = X'(k)ᵀ - D(k) Ĉ(k) ∈ ℝ^{|L(k)| × p}
```

Each row of E(k) is a sample of w_ℓ(k)ᵀ (plus the effect of estimation error
in Ĉ). The sample covariance is:

```
Σ̂ = (1/ν) Σ_{k=0}^{N-1} E(k)ᵀ E(k)
```

where ν is the effective degrees of freedom.

This yields a p × p symmetric positive (semi)definite estimate. For diagonal Σ,
only the diagonal elements are needed:

```
Σ̂_{jj} = (1/ν) Σ_{k=0}^{N-1} Σ_{ℓ ∈ L(k)} e_ℓj(k)²
```

where e_ℓj(k) is the residual for trajectory ℓ, state component j, time step k.

### 6.2 Degrees of Freedom

The exact degrees of freedom for the regularized estimator are:

```
ν = Σ_k |L(k)| - trace(VᵀV A⁻¹)
```

where trace(VᵀV A⁻¹) is the effective number of parameters (analogous to the
hat matrix trace in ridge regression). This satisfies:

```
0 < trace(VᵀV A⁻¹) < N(p+q)
```

with the lower bound approached as λ → ∞ and the upper bound as λ → 0.

Since VᵀV is block diagonal with blocks D(k)ᵀD(k), and A⁻¹ has diagonal blocks
P(k), the trace decomposes as:

```
trace(VᵀV A⁻¹) = Σ_{k=0}^{N-1} trace(D(k)ᵀD(k) × P(k))
```

This is O(N(p+q)²) and can be computed after the backward recursion for P(k),
so the exact degrees of freedom are available at modest cost.

**Simpler conservative approximation:** Use the nominal parameter count:

```
ν_approx = Σ_k |L(k)| - N(p+q)
```

This overestimates the effective parameters (since regularization reduces them),
which underestimates ν, which overestimates Σ̂, which gives conservative
(wider) uncertainty intervals.

### 6.3 Cross-Validation Alternative

An alternative that avoids estimating ν: use the validation loss from the
hyperparameter selection procedure. If λ was chosen by minimizing prediction
error on held-out trajectories, the validation prediction error provides a
direct estimate of Σ that is independent of the degrees-of-freedom calculation.

### 6.4 Diagonal vs. Full Σ

The implementation should support three modes:

| Mode | Estimate | When to use |
|------|----------|-------------|
| Full | Σ̂ = (1/ν) Σ_k E(k)ᵀ E(k) | Default. Captures cross-correlations between state components. |
| Diagonal | Σ̂ = diag(σ̂₁², ..., σ̂ₚ²) | When p is large relative to the number of observations, or when off-diagonal terms are not of interest. |
| Isotropic | Σ̂ = σ̂² Iₚ | Simplest. Equivalent to original COSMIC paper's assumption. |

The choice affects only the final posterior covariance Σ̂ ⊗ P(k), not the
COSMIC algorithm or the P(k) recursion.


## 7. Complete Algorithm for Posterior Covariance

**Inputs:**
- Λ_k for k = 0,...,N-1 (from COSMIC forward pass)
- λ_k for k = 1,...,N-1 (regularization parameters)
- Residuals E(k) for k = 0,...,N-1
- Covariance mode: 'full', 'diagonal', or 'isotropic'

**Step 1: Backward recursion for P(k)**
```
P(N-1) = Λ_{N-1}⁻¹

For k = N-2, ..., 0:
    P(k) = (Λ_k - λ_{k+1}² P(k+1))⁻¹
```

**Step 2: Compute degrees of freedom (exact)**
```
ν = Σ_k |L(k)| - Σ_k trace(D(k)ᵀD(k) P(k))
```

**Step 3: Estimate Σ**
```
Σ̂ = (1/ν) Σ_{k=0}^{N-1} E(k)ᵀ E(k)
```
(apply diagonal or isotropic restriction if requested)

**Step 4: Form posterior covariance at each time step**
```
Cov(vec(C(k))) = Σ̂ ⊗ P(k)
```

**Outputs:**
- P(k) ∈ ℝ^{(p+q)×(p+q)} for each k — row covariance (data/regularization structure)
- Σ̂ ∈ ℝ^{p×p} — estimated noise covariance (column covariance)
- ν — effective degrees of freedom

**Cost:** O(N(p+q)³) for the recursion + O(N(p+q)²) for the trace + O(NLp²) for
the covariance estimate. The dominant cost is the recursion, which matches the
cost of COSMIC itself.


## 8. Extracting Useful Uncertainties

### 8.1 Standard Deviation of Individual Entries

The variance of entry (a, b) of C(k) — the coupling from the a-th component
of [x; u] to the b-th state equation — is:

```
Var(C(k)_{ab}) = Σ̂_{bb} × P(k)_{aa}
```

For a diagonal Σ̂, this is simply σ̂_b² × P(k)_{aa}: the noise level of state
component b times the estimation precision in the a-th parameter direction.

### 8.2 Standard Deviation of A(k) and B(k) Entries

Since C(k) = [A(k)ᵀ; B(k)ᵀ], the entries of A(k)ᵀ occupy rows 1,...,p and
the entries of B(k)ᵀ occupy rows p+1,...,p+q. Therefore:

```
Var(A(k)_{ba}) = Σ̂_{bb} × P(k)_{aa}        for a = 1,...,p
Var(B(k)_{ba}) = Σ̂_{bb} × P(k)_{p+a,p+a}   for a = 1,...,q
```

(Here (b,a) indexes A(k) because C(k) stores A(k)ᵀ.)

### 8.3 Correlation Between State Components

If Σ̂ is not diagonal, the off-diagonal elements give the correlation between
estimation errors in different columns of C(k):

```
Corr(C(k)_{ab}, C(k)_{ad}) = Σ̂_{bd} / √(Σ̂_{bb} Σ̂_{dd})
```

This correlation is the same for all rows a and all time steps k — it reflects
the noise structure, not the data geometry.


## 9. Sanity Checks

### 9.1 No-Regularization Limit (λ → 0)

When λ_k = 0 for all k, the prior is flat, and:
- A = VᵀV (block diagonal)
- A⁻¹ is block diagonal with blocks (D(k)ᵀD(k))⁻¹
- P(k) = (D(k)ᵀD(k))⁻¹
- Cov(vec(C(k))) = Σ ⊗ (D(k)ᵀD(k))⁻¹

This is the standard multivariate OLS covariance (matrix-normal posterior
with known covariance Σ), independently at each time step. ✓

### 9.2 Large-Regularization Limit (λ → ∞)

When λ is very large, the prior forces C(k) ≈ C(k-1) ≈ ... ≈ C(0), i.e., the
system is effectively LTI. All P(k) converge to (Σ_k D(k)ᵀD(k))⁻¹, the
covariance from pooling all data. The posterior covariance shrinks as the
prior provides strong information. ✓

### 9.3 Single Time Step (N = 1)

With N=1 there is no regularization (no F term), and:
- A = D(0)ᵀD(0)
- P(0) = (D(0)ᵀD(0))⁻¹
- Cov(vec(C(0))) = Σ ⊗ (D(0)ᵀD(0))⁻¹

This is ordinary multivariate least squares. ✓

### 9.4 Isotropic Noise (Σ = σ²I)

Reduces to:
- Cov(vec(C(k))) = σ² I ⊗ P(k) = σ² (I ⊗ P(k))
- Columns of C(k) become uncorrelated with identical covariance σ² P(k)

Recovers the isotropic derivation as a special case. ✓

### 9.5 Time Step with No Data

If L(k) = ∅, then D(k)ᵀD(k) = 0, and:
- A_{kk} = (λ_k + λ_{k+1}) I  (only prior information)
- P(k) will be large, reflecting high uncertainty

The posterior mean at k is determined by interpolation from neighbors (via the
backward pass), and the posterior covariance is large. Correct behavior:
no data → wide uncertainty. ✓

### 9.6 Diagonal Σ with σ₁ ≫ σ₂

If state 1 is very noisy and state 2 is clean:
- Σ̂ ≈ diag(σ₁², σ₂²) with σ₁ ≫ σ₂
- Column 1 of C(k) (governing state-1 dynamics) has variance σ₁² P(k)_{aa}
- Column 2 of C(k) (governing state-2 dynamics) has variance σ₂² P(k)_{aa}

The uncertainty is correctly larger for the noisy state and smaller for the
clean state, with the same spatial structure P(k) for both. ✓


## 10. What the Uncertainty Means (and Doesn't Mean)

### 10.1 Correct Interpretation

The posterior is:

```
C(k) | data, Σ̂  ~  MN(Ĉ(k), P(k), Σ̂)
```

Marginal credible intervals for entry (a,b) of C(k) are:

```
Ĉ(k)_{ab} ± z_{α/2} × √(Σ̂_{bb} × P(k)_{aa})
```

These are **Bayesian credible intervals** under the model:
- The noise is i.i.d. N(0, Σ) across trajectories and time steps
- The system matrices follow a Gaussian random walk with covariance
  (1/λ_k) Σ ⊗ I per step

If these assumptions hold, the credible intervals have the stated coverage.

### 10.2 What It Doesn't Cover

- **Model misspecification:** If the true system has nonlinear dynamics,
  the intervals are not calibrated.
- **Wrong λ:** If λ is much larger than the true rate of system variation,
  the intervals will be too narrow (overconfident in smoothness). If λ is
  too small, they will be too wide.
- **Non-Gaussian noise:** The MAP estimate is still the best linear estimator
  under the prior, but credible interval coverage may differ.
- **Correlated noise across trajectories:** If w_ℓ(k) is correlated across ℓ
  (e.g., shared external disturbances), the model is misspecified.
- **Correlated noise across time:** If w_ℓ(k) is correlated across k for a
  given ℓ (e.g., colored process noise), the model is misspecified.

### 10.3 Honest Labeling in sid

When plotting confidence bands derived from the COSMIC posterior, sid should label
them as:

```
"95% posterior credible interval (Bayesian, conditional on λ)"
```

rather than "95% confidence interval," which has a specific frequentist meaning
that is not satisfied here.

In the `sidBodePlot` integration, a pragmatic compromise is to use the label
"uncertainty band" without specifying the statistical framework, similar to how
MathWorks labels the `spa` confidence regions.
