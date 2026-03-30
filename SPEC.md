# sid — Algorithm Specification

**Version:** 0.4.0-draft
**Date:** 2026-03-30
**Reference:** Ljung, L. *System Identification: Theory for the User*, 2nd ed., Prentice Hall, 1999.

---

> **Implementation status:** §1–5 (frequency-domain estimation) and §7 (spectrograms) are implemented. §6 (`sidFreqMap`) is implemented for the BT algorithm; Welch algorithm support is planned. §8 base + §8.4 (`sidLTVdisc`, `sidLTVdiscTune`) are implemented. §8.8–8.12 (variable-length trajectories, Bayesian uncertainty, online/recursive COSMIC, frequency-response lambda tuning, output-only estimation) describe planned features not yet implemented. §9 (`sidFreqETFE`, `sidFreqBTFDR`) is implemented.
> Implementation status: §1–5 (frequency-domain estimation), §6 (sidFreqMap BT + Welch), §7 (spectrograms), §8 base + §8.4 (sidLTVdisc, sidLTVdiscTune), §8.8 (variable-length trajectories), §8.9 (Bayesian uncertainty + sidLTVdiscFrozen), and §9 (sidFreqETFE, sidFreqBTFDR) are implemented. §8.10–8.12 (online/recursive COSMIC, frequency-response lambda tuning, output-only estimation) describe planned features not yet implemented. Multiple trajectory support from the spectral methods are also planned throughout the document.

---

## 1. System Model

All frequency-domain estimation in this package assumes the general linear time-invariant model:

```
y(t) = G(q) u(t) + v(t)       t = 1, 2, ..., N
```

where:

- `y(t)` is the output signal, dimension `n_y × 1`
- `u(t)` is the input signal, dimension `n_u × 1`
- `G(q)` is the transfer function (transfer matrix for MIMO), dimension `n_y × n_u`
- `v(t)` is output disturbance noise, dimension `n_y × 1`, assumed independent of `u(t)`
- `q` is the forward shift operator: `q u(t) = u(t+1)`

The noise `v(t)` may optionally be modeled as filtered white noise:

```
v(t) = H(q) e(t)
```

where `e(t)` is white noise with covariance matrix `Λ`.

**Time series mode:** When no input is present (`n_u = 0`), the model reduces to `y(t) = v(t)` and only the output power spectrum is estimated.

**LTV extension:** The `sidFreqMap` function (§6) relaxes the time-invariance assumption by applying spectral analysis (Blackman-Tukey or Welch) to overlapping segments, producing a time-varying frequency response Ĝ(ω, t). Within each segment, local time-invariance is assumed.

**Multi-trajectory support:** All `sid` functions accept multiple independent trajectories (experiments) of the same system. For frequency-domain functions (`sidFreqBT`, `sidFreqETFE`, `sidFreqMap`, `sidSpectrogram`), spectral estimates are ensemble-averaged across trajectories before forming transfer function ratios or power spectra, reducing variance by a factor of `L` without sacrificing frequency resolution. For `sidLTVdisc`, multiple trajectories are aggregated in the data matrices as described in §8. Multi-trajectory data is passed as 3D arrays `(N × n_ch × L)` when all trajectories share the same length, or as cell arrays `{y1, y2, ..., yL}` when lengths differ. See `docs/multi_trajectory_spectral_theory.md` for the mathematical basis.

---

## 2. `sidFreqBT` — Blackman-Tukey Spectral Analysis

### 2.1 Inputs

| Parameter | Symbol | Type | Default |
|-----------|--------|------|---------|
| Output data | `y` | `(N × n_y)` real matrix | required |
| Input data | `u` | `(N × n_u)` real matrix, or `[]` | `[]` (time series) |
| Window size | `M` | positive integer, `M ≥ 2` | `min(floor(N/10), 30)` |
| Frequencies | `ω` | `(n_f × 1)` vector, rad/sample | 128 points, see §2.2 |
| Sample time | `Ts` | positive scalar (seconds) | `1.0` |

All data must be real-valued and uniformly sampled. If `y` or `u` is a column vector, it is treated as a single channel.

**Multi-trajectory input:** When `y` is `(N × n_y × L)` and `u` is `(N × n_u × L)`, the function computes per-trajectory covariances and averages them before windowing and Fourier transformation:

```
R̂_yu^ens(τ) = (1/L) Σ_{l=1}^{L} R̂_yu^(l)(τ)
```

This ensemble averaging reduces variance by a factor of `L` without affecting frequency resolution. When trajectories have different lengths, pass cell arrays: `y = {y1, y2, ..., yL}`, `u = {u1, u2, ..., uL}`.

### 2.2 Default Frequency Grid

When no frequency vector is specified, the default grid is 128 values **linearly** spaced in `(0, π]`:

```
ω_k = k × π / 128       k = 1, 2, ..., 128
```

in units of rad/sample. To convert to rad/s, divide by `Ts`:

```
ω_k (rad/s) = ω_k (rad/sample) / Ts
```

**Note on returned units:** The result struct stores frequencies in **rad/sample** internally. Plotting functions convert to rad/s using `Ts` when labeling axes.

**Rationale for linear spacing:** The FFT fast path (§2.5) produces linearly spaced frequency bins. Linear spacing is therefore the natural default that enables the FFT optimization. Users who want logarithmic spacing should pass an explicit frequency vector, which triggers the direct DFT path.

### 2.3 Covariance Estimation

Compute the biased sample cross-covariance between signals `x` and `z`, each of length `N`:

```
R̂_xz(τ) = (1/N) Σ_{t=1}^{N-|τ|} x(t+|τ|) z(t)       for τ ≥ 0
R̂_xz(τ) = conj(R̂_zx(-τ))                              for τ < 0
```

The biased estimator (dividing by `N` rather than `N-|τ|`) is used because:
1. It guarantees the resulting spectral estimate is non-negative.
2. It has lower mean-squared error than the unbiased estimator.

For the `sidFreqBT` algorithm, the following covariances are needed for lags `τ = 0, 1, ..., M`:

| Covariance | Signals | Dimensions | Used for |
|------------|---------|------------|----------|
| `R̂_y(τ)` | `y, y` | `n_y × n_y` | Output auto-spectrum |
| `R̂_u(τ)` | `u, u` | `n_u × n_u` | Input auto-spectrum |
| `R̂_yu(τ)` | `y, u` | `n_y × n_u` | Cross-spectrum |

**Time series mode** (`u = []`): Only `R̂_y(τ)` is computed.

**Multi-trajectory covariance:** When `L` trajectories are available, the ensemble-averaged covariance is used:

```
R̂_xz^ens(τ) = (1/L) Σ_{l=1}^{L} R̂_xz^(l)(τ)
```

where `R̂_xz^(l)(τ)` is the biased covariance from trajectory `l`. The averaging is performed at the covariance level, before windowing and Fourier transformation. This preserves the H1 estimator structure (ratio of averaged spectra, not average of ratios).

The Hann (Hanning) window of size `M`:

```
W_M(τ) = 0.5 × (1 + cos(π τ / M))       for |τ| ≤ M
W_M(τ) = 0                                for |τ| > M
```

Properties:
- `W_M(0) = 1`
- `W_M(±M) = 0`
- Symmetric: `W_M(τ) = W_M(-τ)`
- Smooth taper to zero at the edges, reducing spectral leakage

The frequency resolution of the estimate is approximately `2π/M` rad/sample. Larger `M` gives finer resolution but higher variance.

### 2.5 Windowed Spectral Estimates

The spectral estimate at frequency `ω` is the Fourier transform of the windowed covariance:

```
Φ̂_xz(ω) = Σ_{τ=-M}^{M} R̂_xz(τ) × W_M(τ) × exp(-j ω τ)
```

This is computed for all three covariance pairs to produce `Φ̂_y(ω)`, `Φ̂_u(ω)`, and `Φ̂_yu(ω)`.

#### 2.5.1 FFT Fast Path

When using the default frequency grid (§2.2), the computation is done via FFT:

1. Construct the full windowed covariance sequence of length `2M+1`:
   ```
   c(τ) = R̂_xz(τ) × W_M(τ)       for τ = -M, ..., 0, ..., M
   ```

2. Arrange into FFT input order. For a length-`L` FFT where `L ≥ 2M+1`:
   ```
   s(k) = c(k)           for k = 0, 1, ..., M
   s(k) = 0              for k = M+1, ..., L-M-1     (zero-padding)
   s(k) = c(k - L)       for k = L-M, ..., L-1       (negative lags wrapped)
   ```
   In practice, `L = 256` (the smallest power of 2 ≥ 2×128+1, used for the default 128-frequency grid).

3. Compute `S = fft(s)`.

4. Extract the desired frequency bins: `Φ̂(ω_k) = S(k+1)` for `k = 1, ..., 128`
   (MATLAB 1-indexed: bin 1 is DC, bin `k+1` corresponds to frequency `k × 2π/L`).

**Scaling:** No additional scaling factor is applied. The FFT computes the sum directly.

#### 2.5.2 Direct DFT Path

When the user supplies a custom frequency vector `ω`, compute the sum explicitly:

```
Φ̂_xz(ω) = R̂_xz(0) × W_M(0) + 2 × Σ_{τ=1}^{M} Re[ R̂_xz(τ) × W_M(τ) × exp(-j ω τ) ]
```

Wait — this shortcut is only valid when `R̂_xz(τ)` is the auto-covariance of a real signal (which is real and symmetric). For the cross-covariance `R̂_yu(τ)`, the full complex sum must be used:

```
Φ̂_yu(ω) = Σ_{τ=-M}^{M} R̂_yu(τ) × W_M(τ) × exp(-j ω τ)
```

where `R̂_yu(-τ) = R̂_uy(τ)' = conj(R̂_uy(τ))` for the scalar case.

**Implementation:** For each frequency `ω_k`, compute:

```
Φ̂_xz(ω_k) = W_M(0) × R̂_xz(0) + Σ_{τ=1}^{M} W_M(τ) × [ R̂_xz(τ) × exp(-j ω_k τ)
                                                             + conj(R̂_xz(τ)) × exp(+j ω_k τ) ]
```

which for real signals and auto-covariances simplifies to:

```
Φ̂_xx(ω_k) = W_M(0) × R̂_xx(0) + 2 × Σ_{τ=1}^{M} W_M(τ) × R̂_xx(τ) × cos(ω_k τ)
```

This form is real-valued and non-negative, as expected for a power spectrum.

### 2.6 Frequency Response Estimate

**SISO case:**

```
Ĝ(ω) = Φ̂_yu(ω) / Φ̂_u(ω)
```

**MIMO case** (`n_u > 1`):

```
Ĝ(ω) = Φ̂_yu(ω) × Φ̂_u(ω)^{-1}
```

where `Φ̂_yu(ω)` is `n_y × n_u` and `Φ̂_u(ω)` is `n_u × n_u`. The matrix inverse is computed independently at each frequency.

**Regularization:** If `Φ̂_u(ω)` is singular or nearly singular at some frequency `ω_k`:
- SISO: if `|Φ̂_u(ω_k)| < ε × max(|Φ̂_u|)` where `ε = 1e-10`, set `Ĝ(ω_k) = NaN + j×NaN`.
- MIMO: if `cond(Φ̂_u(ω_k)) > 1/ε`, set the corresponding row of `Ĝ(ω_k)` to `NaN`.
- Issue a warning when this occurs.

### 2.7 Noise Spectrum Estimate

**SISO case:**

```
Φ̂_v(ω) = Φ̂_y(ω) - |Φ̂_yu(ω)|² / Φ̂_u(ω)
```

**MIMO case:**

```
Φ̂_v(ω) = Φ̂_y(ω) - Φ̂_yu(ω) × Φ̂_u(ω)^{-1} × Φ̂_yu(ω)'
```

where `'` denotes conjugate transpose.

**Non-negativity:** Due to estimation errors, `Φ̂_v(ω)` may become slightly negative at some frequencies. Clamp to zero:

```
Φ̂_v(ω) = max(Φ̂_v(ω), 0)      (SISO)
```

For MIMO, ensure the matrix is positive semi-definite by zeroing any negative eigenvalues.

**Time series mode:** No noise spectrum is computed separately. The output spectrum `Φ̂_y(ω)` is returned in the `NoiseSpectrum` field.

### 2.8 Normalization

The spectral estimates use the following normalization:

```
Φ̂(ω) = Σ_{τ=-M}^{M} R̂(τ) W_M(τ) exp(-jωτ)
```

This matches the System Identification Toolbox convention. It does **not** include:
- A factor of `Ts` (the Signal Processing Toolbox convention includes `Ts`)
- A factor of `1/(2π)`

To convert to the Signal Processing Toolbox convention, multiply by `Ts`:

```
Φ̂_SPT(ω) = Ts × Φ̂_SID(ω)
```

---

## 3. Uncertainty Estimation

### 3.1 Window Norm

Define the squared window norm:

```
C_W = Σ_{τ=-M}^{M} W_M(τ)²
```

For the Hann window, this evaluates to:

```
C_W = 1 + 2 × Σ_{τ=1}^{M} [0.5 × (1 + cos(πτ/M))]²
```

which can be computed in closed form as `C_W = (3/4)×(2M) + 1/2 = (3M + 1)/2`, but the implementation should compute it numerically from the actual window values to avoid any discrepancy.

### 3.2 Coherence

The squared coherence between input and output:

```
γ̂²(ω) = |Φ̂_yu(ω)|² / (Φ̂_y(ω) × Φ̂_u(ω))
```

This is real-valued and satisfies `0 ≤ γ̂²(ω) ≤ 1`. Values near 1 indicate the output is well explained by the input at that frequency; values near 0 indicate noise dominates.

### 3.3 Variance of the Frequency Response

The asymptotic variance of the frequency response estimate (Ljung 1999, p. 184):

```
Var{Ĝ(ω)} ≈ (C_W / N) × |Ĝ(ω)|² × (1 - γ̂²(ω)) / γ̂²(ω)
```

The standard deviation returned in the result struct is:

```
σ_G(ω) = sqrt(Var{Ĝ(ω)})
```

**Regularization:** If `γ̂²(ω_k) < ε` (where `ε = 1e-10`), set `σ_G(ω_k) = Inf`. This corresponds to frequencies where the input has negligible power and the estimate is unreliable.

**Note:** This formula gives the variance of the complex-valued `Ĝ`. The standard deviation `σ_G` applies equally to real and imaginary parts. Confidence intervals for magnitude are constructed as:

```
|Ĝ(ω)| ± p × σ_G(ω)
```

where `p` is the number of standard deviations (default: 3 for ≈99.7% coverage under Gaussian assumptions).

**Multi-trajectory variance:** When `L` trajectories are ensemble-averaged, the variance is reduced by a factor of `L`:

```
Var{Ĝ^ens(ω)} ≈ (C_W / (L × N)) × |Ĝ(ω)|² × (1 - γ̂²(ω)) / γ̂²(ω)
```

The coherence `γ̂²` is now the ensemble coherence, which is generally higher than any single-trajectory coherence because the noise averages out while the signal accumulates.

### 3.4 Variance of the Noise Spectrum

The asymptotic variance of the spectral estimate (Ljung 1999, p. 188):

```
Var{Φ̂_v(ω)} ≈ (2 × C_W / N) × Φ̂_v(ω)²
```

Standard deviation:

```
σ_Φv(ω) = sqrt(Var{Φ̂_v(ω)})
```

### 3.5 Variance of the Output Spectrum (Time Series Mode)

When no input is present:

```
Var{Φ̂_y(ω)} ≈ (2 × C_W / N) × Φ̂_y(ω)²
```

This is the standard asymptotic result for windowed spectral estimates.

---

## 4. `sidFreqETFE` — Empirical Transfer Function Estimate

### 4.1 Algorithm

The ETFE is the ratio of the output and input discrete Fourier transforms:

```
Ĝ_ETFE(ω_k) = Y(ω_k) / U(ω_k)
```

where:

```
Y(ω_k) = Σ_{t=1}^{N} y(t) exp(-j ω_k t)
U(ω_k) = Σ_{t=1}^{N} u(t) exp(-j ω_k t)
```

This is equivalent to `sidFreqBT` with window size `M = N` (rectangular window). It provides the maximum frequency resolution but has high variance.

**Multi-trajectory ETFE:** When `L` trajectories are available, the cross-periodograms are averaged before forming the ratio:

```
Ĝ_ETFE^ens(ω_k) = Φ̂_yu^ens(ω_k) / Φ̂_u^ens(ω_k)
```

where `Φ̂_yu^ens(ω_k) = (1/L) Σ_l Y_l(ω_k) conj(U_l(ω_k))`. This is the multi-trajectory H1 estimator, reducing variance by a factor of `L`.

### 4.2 Optional Smoothing

A smoothing parameter `S` (positive odd integer) may be specified. When given, the raw ETFE is convolved with a length-`S` rectangular (boxcar) frequency-domain window:

```
Ĝ_smooth(ω_k) = (1/S) × Σ_{j=-(S-1)/2}^{(S-1)/2} Ĝ_ETFE(ω_{k+j})
```

with appropriate handling at the boundaries.

### 4.3 Noise Spectrum

For the ETFE, the noise spectrum estimate is the periodogram of the residuals:

```
Φ̂_v(ω_k) = (1/N) × |Y(ω_k) - Ĝ(ω_k) × U(ω_k)|²
```

### 4.4 Time Series Mode

When no input is present, the ETFE reduces to the **periodogram**:

```
Φ̂_y(ω_k) = (1/N) × |Y(ω_k)|²
```

---

## 5. `sidFreqBTFDR` — Frequency-Dependent Resolution

### 5.1 Concept

`sidFreqBTFDR` is identical to `sidFreqBT` except that the window size `M` varies with frequency, allowing different resolution at different frequencies. The user specifies a **resolution** parameter `R(ω)` (in rad/sample) instead of a window size.

### 5.2 Resolution to Window Size Mapping

At each frequency `ω_k`, the local window size is:

```
M_k = ceil(2π / R_k)
```

where `R_k = R(ω_k)` is the desired resolution at that frequency.

If `R` is a scalar, it applies uniformly. If `R` is a vector of the same length as the frequency grid, each entry specifies the local resolution.

### 5.3 Algorithm

For each frequency `ω_k`:

1. Determine `M_k` from the resolution.
2. Compute the Hann window `W_{M_k}(τ)` of size `M_k`.
3. Compute the windowed spectral estimates `Φ̂_y(ω_k)`, `Φ̂_u(ω_k)`, `Φ̂_yu(ω_k)` using the direct DFT formula with window size `M_k`.
4. Form `Ĝ(ω_k)` and `Φ̂_v(ω_k)` as in §2.6 and §2.7.

**Note:** The FFT fast path cannot be used here because the window size varies across frequencies. All computations use the direct DFT.

### 5.4 Default Resolution

If no resolution is specified:

```
R = 2π / min(floor(N/10), 30)
```

This matches the default behavior of `sidFreqBT`.

---

## 6. `sidFreqMap` — Time-Varying Frequency Response Map

### 6.1 Concept

`sidFreqMap` estimates a **time-varying frequency response** Ĝ(ω, t) by applying spectral analysis to overlapping segments of input-output data. This reveals how the system's transfer function, noise spectrum, and coherence evolve over time.

Two algorithms are supported via the `'Algorithm'` parameter:

| Algorithm | Method | Replaces | Within each segment |
|-----------|--------|----------|---------------------|
| `'bt'` (default) | Blackman-Tukey correlogram | `spa` applied per segment | Covariance → lag window → DFT |
| `'welch'` | Welch's averaged periodogram | MathWorks `tfestimate` | Sub-segments → time-domain window → FFT → average → form ratios |

Both produce identical output structures: Ĝ(ω, t), Φ̂_v(ω, t), γ̂²(ω, t). The choice affects the bias-variance tradeoff within each segment, not the user-facing interface.

For an LTI system, the map is constant along the time axis — this serves as a diagnostic check. For an LTV (linear time-varying) system, the map shows modes appearing, disappearing, shifting in frequency, or changing in gain.

This extends the `spectrogram` concept from single-signal time-frequency analysis to **input-output system identification**:

| Tool | Input | Output | Shows |
|------|-------|--------|-------|
| `spectrogram` / `sidSpectrogram` | One signal | \|X(ω,t)\|² | How signal frequency content changes |
| `sidFreqMap` | Input + output pair | Ĝ(ω,t), Φ̂_v(ω,t), γ̂²(ω,t) | How the *system itself* changes |
| `sidFreqMap` | One signal (time series) | Φ̂_y(ω,t) | How signal spectrum changes (≈ spectrogram) |

When used together, `sidSpectrogram` on `u` and `y` alongside `sidFreqMap` on the pair `(y, u)` provides a complete diagnostic picture: the input's spectral content, the output's spectral content, and the system connecting them — all on aligned time axes.

### 6.2 Inputs

| Parameter | Symbol | Type | Default |
|-----------|--------|------|---------|
| Output data | `y` | `(N × n_y)` real matrix | required |
| Input data | `u` | `(N × n_u)` real matrix, or `[]` | `[]` (time series) |
| Segment length | `L` | positive integer | `min(floor(N/4), 256)` |
| Overlap | `P` | integer, `0 ≤ P < L` | `floor(L/2)` (50% overlap) |
| Algorithm | | `'bt'` or `'welch'` | `'bt'` |
| Sample time | `Ts` | positive scalar (seconds) | `1.0` |

**Algorithm-specific parameters:**

| Parameter | Applies to | Type | Default |
|-----------|-----------|------|---------|
| `WindowSize` (M) | `'bt'` only | positive integer | `min(floor(L/10), 30)` |
| `Frequencies` | `'bt'` only | `(n_f × 1)` vector | 128 linearly spaced in (0, π] |
| `SubSegmentLength` | `'welch'` only | positive integer | `floor(L/4.5)` (matches `tfestimate` default) |
| `SubOverlap` | `'welch'` only | non-negative integer | `floor(SubSegmentLength / 2)` |
| `Window` | `'welch'` only | `'hann'`, `'hamming'`, or vector | `'hann'` |
| `NFFT` | `'welch'` only | positive integer | `max(256, 2^nextpow2(SubSegmentLength))` |

**Multi-trajectory input:** When `y` is `(N × n_y × L)` and `u` is `(N × n_u × L)`, spectral estimates within each segment are ensemble-averaged across trajectories before forming transfer function ratios. For variable-length trajectories, pass cell arrays. At each segment `k`, only trajectories that span segment `k` contribute to the ensemble. This directly parallels COSMIC's multi-trajectory aggregation (§8.3.2), ensuring consistent use of the same data across time-domain and frequency-domain analyses.

### 6.3 Outer Segmentation (Common to Both Algorithms)

Both algorithms share the same outer segmentation:

1. Divide the data into `K` overlapping segments, each of length `L` samples, with overlap `P`:
   ```
   Segment k: samples (k-1)(L-P)+1  through  (k-1)(L-P)+L
   for k = 1, 2, ..., K
   where K = floor((N - L) / (L - P)) + 1
   ```

2. For each segment `k`, extract `y_k = y(start:end, :)` and `u_k = u(start:end, :)`.

3. Apply the selected algorithm to estimate `Ĝ(ω)`, `Φ̂_v(ω)`, `γ̂²(ω)` within the segment.

4. Collect the per-segment results into time-frequency arrays.

### 6.4 Inner Estimation: Blackman-Tukey (`'bt'`)

Within each segment of length `L`, apply `sidFreqBT`:

1. Compute biased covariances `R̂_y(τ)`, `R̂_u(τ)`, `R̂_yu(τ)` for lags `0..M`.
2. Apply Hann lag window `W_M(τ)`.
3. Fourier transform to obtain `Φ̂_y(ω)`, `Φ̂_u(ω)`, `Φ̂_yu(ω)`.
4. Form `Ĝ(ω) = Φ̂_yu(ω) / Φ̂_u(ω)`.
5. Form `Φ̂_v(ω) = Φ̂_y(ω) - |Φ̂_yu(ω)|² / Φ̂_u(ω)`.
6. Compute coherence `γ̂²(ω) = |Φ̂_yu(ω)|² / (Φ̂_y(ω) Φ̂_u(ω))`.
7. Compute asymptotic uncertainty via `sidUncertainty`.

**Frequency resolution** within the segment is controlled by the lag window size `M`. The constraint `L > 2M` must hold.

### 6.5 Inner Estimation: Welch (`'welch'`)

Within each segment of length `L`, apply the Welch method (equivalent to `tfestimate` + `mscohere` + `cpsd`):

1. Divide the segment into `J` overlapping sub-segments of length `L_sub` with overlap `P_sub`:
   ```
   J = floor((L - L_sub) / (L_sub - P_sub)) + 1
   ```

2. For each sub-segment `j`:
   a. Apply the time-domain window `w(n)` (Hann by default):
      ```
      y_j(n) = y_segment(n_start + n) × w(n)
      u_j(n) = u_segment(n_start + n) × w(n)
      ```
   b. Compute FFTs: `Y_j(m) = FFT(y_j)`, `U_j(m) = FFT(u_j)`.

3. Average the cross-spectral and auto-spectral periodograms:
   ```
   Φ̂_yu(ω) = (1/J) Σ_j Y_j(ω) conj(U_j(ω)) / S₁
   Φ̂_u(ω)  = (1/J) Σ_j |U_j(ω)|² / S₁
   Φ̂_y(ω)  = (1/J) Σ_j |Y_j(ω)|² / S₁
   ```
   where `S₁ = Σ_n w(n)²` is the window power normalization.

4. Form `Ĝ(ω) = Φ̂_yu(ω) / Φ̂_u(ω)`.
5. Form `Φ̂_v(ω)` and `γ̂²(ω)` as in the BT case.

**Frequency resolution** is determined by the sub-segment length `L_sub` and the NFFT: `Δf = Fs / NFFT`. The sub-segment overlap `P_sub` controls variance reduction — more sub-segments (higher overlap) → lower variance but no change in resolution.

**Uncertainty:** The variance of the Welch spectral estimate is approximately:

```
Var{Φ̂(ω)} ≈ Φ²(ω) / ν
```

where `ν = 2J × (1 - c_overlap)` is the equivalent degrees of freedom, and `c_overlap` is a correction factor depending on the overlap ratio and window shape. For 50% overlap with a Hann window, `ν ≈ 1.8J`.

### 6.6 Comparison of BT and Welch

| Aspect | BT (`sidFreqBT`) | Welch |
|--------|-------------------|-------|
| Resolution control | Lag window size `M` | Sub-segment length `L_sub` |
| Variance control | `M` (smaller M → lower variance) | Number of sub-segments `J` (more → lower variance) |
| Guaranteed non-negative spectrum | Yes (biased covariance estimator) | Yes (averaged periodograms) |
| Custom frequency grid | Yes (direct DFT path) | No (FFT bins only) |
| Normalization | System ID convention (no Ts factor) | PSD convention (includes Ts) |
| Best for | Smooth spectra, custom frequencies | Standard analysis, `tfestimate` compatibility |

**Default choice:** `'bt'` is the default because it matches the `sid` package's primary use case (system identification with `sidFreqBT`-compatible output) and supports custom frequency grids. Users coming from `tfestimate` should use `'welch'`.

### 6.7 Time Vector

The center time of each segment defines the time axis:

```
t_k = ((k-1)(L-P) + L/2) × Ts       for k = 1, ..., K
```

in units of seconds.

### 6.8 Output Struct

`sidFreqMap` returns a struct with fields:

| Field | Type | Description |
|-------|------|-------------|
| `Time` | `(K × 1)` real | Center time of each segment (seconds) |
| `Frequency` | `(n_f × 1)` real | Frequency vector (rad/sample) |
| `FrequencyHz` | `(n_f × 1)` real | Frequency vector (Hz) |
| `Response` | `(n_f × K)` complex | Time-varying frequency response Ĝ(ω, t) |
| `ResponseStd` | `(n_f × K)` real | Standard deviation of Ĝ per segment |
| `NoiseSpectrum` | `(n_f × K)` real | Time-varying noise spectrum Φ̂_v(ω, t) |
| `NoiseSpectrumStd` | `(n_f × K)` real | Standard deviation per segment |
| `Coherence` | `(n_f × K)` real | Time-varying squared coherence γ̂²(ω, t) |
| `SampleTime` | scalar | Sample time Ts |
| `SegmentLength` | scalar | Segment length L |
| `Overlap` | scalar | Overlap P |
| `WindowSize` | scalar | BT lag window size M (BT only) |
| `Algorithm` | char | `'bt'` or `'welch'` |
| `NumTrajectories` | scalar or `(K × 1)` | Number of trajectories used (scalar if constant, vector if variable-length) |
| `Method` | char | `'sidFreqMap'` |

**Dimensions shown are for SISO.** For MIMO, `Response` becomes `(n_f × K × n_y × n_u)`, etc.

The output struct is identical regardless of algorithm, so `sidMapPlot` and downstream tools (including COSMIC lambda cross-validation in §8.11) work transparently with either.

### 6.9 Visualization: `sidMapPlot`

The natural visualization is a **color map** (like a spectrogram):

- **x-axis:** Time (seconds)
- **y-axis:** Frequency (rad/s or Hz, log scale)
- **Color:** Magnitude of Ĝ(ω, t) in dB, or Φ̂_v(ω, t) in dB, or γ̂²(ω, t)

The function `sidMapPlot` provides selectable plot types via a `'PlotType'` option:

| PlotType | Color represents | Use case |
|----------|-----------------|----------|
| `'magnitude'` (default) | `20 log10(\|Ĝ(ω,t)\|)` | Track gain changes |
| `'phase'` | `angle(Ĝ(ω,t))` in degrees | Track phase drift |
| `'noise'` | `10 log10(Φ̂_v(ω,t))` | Track disturbance evolution |
| `'coherence'` | `γ̂²(ω,t)` on [0, 1] | Identify when LTI assumption breaks down |
| `'spectrum'` | `10 log10(Φ̂_y(ω,t))` | Time series mode (equivalent to spectrogram) |

### 6.10 Compatibility with MathWorks `tfestimate`

`sidFreqMap` with `'Algorithm', 'welch'` replicates the core functionality of the Signal Processing Toolbox `tfestimate`, `mscohere`, and `cpsd` functions. Specifically:

```matlab
% MathWorks style (single-window transfer function estimate):
[Txy, F] = tfestimate(u, y, hann(256), 128, 512, Fs);
[Cxy, F] = mscohere(u, y, hann(256), 128, 512, Fs);

% sid equivalent (time-varying, but with segment = full data → single estimate):
result = sidFreqMap(y, u, 'Algorithm', 'welch', ...
                         'SegmentLength', length(y), ...
                         'SubSegmentLength', 256, ...
                         'SubOverlap', 128, ...
                         'NFFT', 512, ...
                         'SampleTime', 1/Fs);
% result.Response ≈ Txy, result.Coherence ≈ Cxy
```

The key difference: `sidFreqMap` always produces time-varying output. Setting `SegmentLength` equal to the data length reduces it to a single-window estimate equivalent to `tfestimate`.

### 6.11 Design Considerations

**Segment length vs. inner parameters:** The outer segment length `L` determines the temporal resolution of the map (how finely you resolve changes in time). The inner parameters (`M` for BT, `L_sub` for Welch) control frequency resolution and variance within each segment. These are independent choices.

**Computational cost:** `K` calls to the inner estimator. For BT, each is O(L×M + M×n_f). For Welch, each is O(J×L_sub×log(L_sub)). Both are fast for typical parameters.

**Edge effects:** The first and last segments may produce less reliable estimates if the system is non-stationary near the boundaries. No special handling is applied — the uncertainty estimates from each segment naturally reflect the reduced confidence.

---

## 7. `sidSpectrogram` — Short-Time Spectral Analysis

### 7.1 Purpose

`sidSpectrogram` computes the short-time Fourier transform (STFT) spectrogram of one or more signals. It replicates the core functionality of the Signal Processing Toolbox `spectrogram` function, with two additional roles in the `sid` workflow:

1. **Diagnostic companion to `sidFreqMap`.** Plotting the spectrograms of `y` and `u` alongside the time-varying transfer function map lets the user distinguish genuine system changes from input-driven effects. If a spectral feature appears in both the `y` spectrogram and the Ĝ(ω,t) map but *not* in the `u` spectrogram, it's likely a real system change. If it appears in `u` too, it's the input driving the output.

2. **Standalone time-frequency analysis** for users who don't have the Signal Processing Toolbox.

### 7.2 Inputs

| Parameter | Symbol | Type | Default |
|-----------|--------|------|---------|
| Signal | `x` | `(N × n_ch)` real matrix | required |
| Window length | `L` | positive integer | `256` |
| Overlap | `P` | integer, `0 ≤ P < L` | `floor(L/2)` |
| NFFT | `nfft` | positive integer | `max(256, 2^nextpow2(L))` |
| Window function | `win` | `'hann'`, `'hamming'`, `'rect'`, or `(L × 1)` vector | `'hann'` |
| Sample time | `Ts` | positive scalar (seconds) | `1.0` |

**Note on window terminology:** The window here is a **time-domain** tapering window applied to each data segment before FFT — this is distinct from the **lag-domain** Hann window used in `sidFreqBT`. The spectrogram window reduces spectral leakage; the BT lag window controls frequency resolution of the correlogram.

**Multi-trajectory input:** When `x` is `(N × n_ch × L)`, the power spectral density within each segment is averaged across trajectories:

```
P̂^ens(ω, t_k) = (1/L) Σ_l P̂^(l)(ω, t_k)
```

This is the event-related spectral perturbation (ERSP) approach, standard in neuroscience and vibration analysis. It reduces noise while preserving time-locked spectral features that are consistent across realizations. For variable-length trajectories, pass cell arrays.

### 7.3 Algorithm

The standard short-time Fourier transform:

1. Divide the signal `x` into `K` overlapping segments of length `L`, with overlap `P`:
   ```
   x_k(n) = x((k-1)(L-P) + n) × w(n)       n = 1, ..., L
   ```
   where `w(n)` is the time-domain window and `K = floor((N - L) / (L - P)) + 1`.

2. Compute the FFT of each windowed segment:
   ```
   X_k(m) = Σ_{n=1}^{L} x_k(n) × exp(-j 2π (m-1) n / nfft)       m = 1, ..., nfft
   ```

3. Compute the one-sided power spectral density for each segment:
   ```
   P_k(m) = (1 / (Fs × S₁)) × |X_k(m)|²
   ```
   where `S₁ = Σ w(n)²` is the window power, and `Fs = 1/Ts`. For one-sided spectra, the positive-frequency bins (excluding DC and Nyquist) are doubled.

4. The spectrogram is the matrix `P(m, k)` for `m = 1, ..., nfft/2+1` and `k = 1, ..., K`.

### 7.4 Output Struct

| Field | Type | Description |
|-------|------|-------------|
| `Time` | `(K × 1)` real | Center time of each segment (seconds) |
| `Frequency` | `(n_bins × 1)` real | Frequency vector (Hz) |
| `FrequencyRad` | `(n_bins × 1)` real | Frequency vector (rad/s) |
| `Power` | `(n_bins × K × n_ch)` real | Power spectral density per segment |
| `PowerDB` | `(n_bins × K × n_ch)` real | `10 × log10(Power)` |
| `Complex` | `(n_bins × K × n_ch)` complex | Complex STFT coefficients (before squaring) |
| `SampleTime` | scalar | Sample time Ts |
| `WindowLength` | scalar | Segment length L |
| `Overlap` | scalar | Overlap P |
| `NFFT` | scalar | FFT length |
| `Method` | char | `'sidSpectrogram'` |

where `n_bins = floor(nfft/2) + 1` (one-sided spectrum).

### 7.5 Visualization

`sidSpectrogram` can be plotted using `sidMapPlot` with `'PlotType', 'spectrum'`, or with a dedicated call:

```matlab
result = sidSpectrogram(y, 'WindowLength', 256, 'Overlap', 128);
sidSpectrogramPlot(result);
```

`sidSpectrogramPlot` produces a standard spectrogram color map:

- **x-axis:** Time (seconds)
- **y-axis:** Frequency (Hz), linear or log scale
- **Color:** Power in dB

### 7.6 Relationship to `sidFreqMap`

The two functions share segmentation conventions (segment length, overlap, time vector computation) so their time axes align when called with the same parameters. A typical diagnostic workflow:

```matlab
% Same segmentation parameters for alignment
L = 256; P = 128; Ts = 0.001;

% Spectrograms of raw signals
specY = sidSpectrogram(y, 'WindowLength', L, 'Overlap', P, 'SampleTime', Ts);
specU = sidSpectrogram(u, 'WindowLength', L, 'Overlap', P, 'SampleTime', Ts);

% Time-varying transfer function
mapG = sidFreqMap(y, u, 'SegmentLength', L, 'Overlap', P, 'SampleTime', Ts);

% Compare side-by-side
figure;
subplot(3,1,1); sidSpectrogramPlot(specU); title('Input u');
subplot(3,1,2); sidSpectrogramPlot(specY); title('Output y');
subplot(3,1,3); sidMapPlot(mapG, 'PlotType', 'magnitude'); title('G(w,t)');
```

This layout immediately reveals whether spectral features in the output are input-driven or system-driven.

### 7.7 Compatibility with MathWorks `spectrogram`

The MathWorks `spectrogram` function uses the calling convention `spectrogram(x, window, noverlap, nfft, fs)`. `sidSpectrogram` supports a compatible positional syntax:

```matlab
% MathWorks style:
[S, F, T, P] = spectrogram(x, hann(256), 128, 512, 1000);

% sid equivalent:
result = sidSpectrogram(x, 'WindowLength', 256, 'Overlap', 128, ...
                           'NFFT', 512, 'SampleTime', 0.001);
% result.Complex ≈ S, result.Frequency ≈ F, result.Time ≈ T, result.Power ≈ P
```

The normalization follows the PSD convention (power per unit frequency), matching the MathWorks default when `spectrogram` is called with the `'psd'` option.

---

## 8. `sidLTVdisc` — Discrete-Time LTV State-Space Identification

### 8.1 Problem Statement

Identify the time-varying system matrices of a discrete linear time-varying system:

```
x(k+1) = A(k) x(k) + B(k) u(k)       k = 0, 1, ..., N-1
```

where `x(k) ∈ ℝᵖ` is the state, `u(k) ∈ ℝᵍ` is the control input, `A(k) ∈ ℝᵖˣᵖ` and `B(k) ∈ ℝᵖˣᵍ` are the unknown time-varying system matrices.

Given measured state trajectories `X` and control inputs `U`, estimate `A(k)` and `B(k)` for all `k`.

### 8.2 Inputs

| Parameter | Symbol | Type | Default |
|-----------|--------|------|---------|
| State data | `X` | `(N+1 × p)` or `(N+1 × p × L)` | required |
| Input data | `U` | `(N × q)` or `(N × q × L)` | required |
| Regularization | `λ` | scalar, `(N-1 × 1)` vector, or `'auto'` | `'auto'` |
| Algorithm | | `'cosmic'` | `'cosmic'` |
| Precondition | | logical | `false` |

Here `L` is the number of trajectories. All trajectories must have the same horizon `N+1`.

### 8.3 COSMIC Algorithm

**Reference:** Carvalho, Soares, Lourenço, Ventura. "COSMIC: fast closed-form identification from large-scale data for LTV systems." arXiv:2112.04355, 2022.

#### 8.3.1 Optimization Variable

Define the stacked optimization variable:

```
C(k) = [A(k)ᵀ; B(k)ᵀ] ∈ ℝ⁽ᵖ⁺ᵍ⁾ˣᵖ       k = 0, ..., N-1
```

#### 8.3.2 Data Matrices

For `L` trajectories at time step `k`:

```
D(k) = [X(k)ᵀ  U(k)ᵀ] ∈ ℝᴸˣ⁽ᵖ⁺ᵍ⁾       (data matrix)
X'(k) = X(k+1)ᵀ ∈ ℝᴸˣᵖ                    (next-state matrix)
```

where `X(k) = [x₁(k), x₂(k), ..., x_L(k)]` collects states from all trajectories.

#### 8.3.3 Cost Function

```
f(C) = (1/2) Σ_{k=0}^{N-1} ||D(k)C(k) - X'(k)||²_F
     + (1/2) Σ_{k=1}^{N-1} ||λ_k^{1/2} (C(k) - C(k-1))||²_F
```

The first term is **data fidelity**: how well the model predicts next states across all trajectories. The second term is **temporal smoothness**: penalizes large changes in system matrices between consecutive time steps.

`λ_k > 0` is the regularization strength at time step `k`. Higher `λ_k` → smoother transitions (system changes slowly). Lower `λ_k` → more freedom for rapid changes.

#### 8.3.4 Closed-Form Solution

Setting ∇f(C) = 0 yields a **block tridiagonal** linear system. Define:

```
S_00         = D(0)ᵀD(0) + λ₁ I
S_{N-1,N-1}  = D(N-1)ᵀD(N-1) + λ_{N-1} I
S_kk         = D(k)ᵀD(k) + (λ_k + λ_{k+1}) I     for k = 1, ..., N-2
Θ_k          = D(k)ᵀ X'(k)ᵀ                         for k = 0, ..., N-1
```

**Forward pass** (k = 0 to N-1):

```
Λ₀ = S_00
Y₀ = Λ₀⁻¹ Θ₀

For k = 1, ..., N-1:
    Λ_k = S_kk - λ_k² Λ_{k-1}⁻¹
    Y_k = Λ_k⁻¹ (Θ_k + λ_k Y_{k-1})
```

**Backward pass** (k = N-2 to 0):

```
C(N-1) = Y_{N-1}

For k = N-2, ..., 0:
    C(k) = Y_k + λ_{k+1} Λ_k⁻¹ C(k+1)
```

**Complexity:** `O(N × (p+q)³)` — linear in the number of time steps, cubic in state+input dimension, independent of the number of trajectories `L` (which only affects the precomputation of `D(k)ᵀD(k)` and `Θ_k`).

#### 8.3.5 Existence and Uniqueness

A unique solution exists if and only if the empirical covariance of the data is positive definite:

```
Σ = Σ₁ + Σ₂ + ... + Σ_L ≻ 0
```

where:

```
Σ_ℓ = (1/N) Σ_{k=0}^{N} [x_ℓ(k); u_ℓ(k)] [x_ℓ(k); u_ℓ(k)]ᵀ
```

Equivalently, the complete set of `[x_ℓ(k)ᵀ  u_ℓ(k)ᵀ]` vectors across all trajectories and time steps must span `ℝᵖ⁺ᵍ`.

#### 8.3.6 Preconditioning

When data matrices `D(k)ᵀD(k)` are ill-conditioned, preconditioning improves numerical stability by redefining:

```
S_kk^PC = I
S_ij^PC = S_kk⁻¹ S_ij         for i ≠ j
Θ_k^PC  = S_kk⁻¹ Θ_k
```

This rescales each block row of the tridiagonal system to have identity on the diagonal, reducing the condition number of the matrices that need to be inverted.

### 8.4 Lambda Selection

#### 8.4.1 Manual

The user provides `λ` as a scalar (applied uniformly) or as an `(N-1 × 1)` vector (per-step).

#### 8.4.2 L-Curve (Automatic)

When `'Lambda', 'auto'` is specified, `sidLTVdisc` selects λ using the L-curve method:

1. Define a grid of candidate values: `λ_grid = logspace(-3, 15, 50)`.
2. For each candidate `λ_j`, run COSMIC and record:
   - Data fidelity: `F_j = ||VC - X'||²_F`
   - Regularization: `R_j = Σ ||λ^{1/2}(C(k) - C(k-1))||²_F`
3. Plot `log(R_j)` vs. `log(F_j)`. This traces an L-shaped curve.
4. Select the λ at the corner of the L — the point of maximum curvature:
   ```
   κ_j = |F''_j R'_j - F'_j R''_j| / (F'_j² + R'_j²)^{3/2}
   ```
   where derivatives are computed by finite differences along the curve.

The L-curve method requires multiple COSMIC runs, but each is O(N(p+q)³), so the total cost is typically under a second for moderate problems.

#### 8.4.3 Validation-Based Tuning (`sidLTVdiscTune`)

A separate function that wraps `sidLTVdisc` in a grid search over λ, evaluating trajectory prediction loss on validation data:

```matlab
function [bestResult, bestLambda, allLosses] = sidLTVdiscTune(X_train, U_train, X_val, U_val, varargin)
```

**Trajectory prediction loss** (from the COSMIC paper):

```
L(λ) = (1/|S|) Σ_{ℓ∈S} sqrt( (1/N) Σ_{k=1}^{N} Σ_{m=1}^{p} (x̂_km^(ℓ)(λ) - x_km^(ℓ))² )
```

where `x̂` is the state predicted by propagating the identified model from initial conditions, and `S` is the set of validation trajectories.

**Inputs:**

| Parameter | Type | Default |
|-----------|------|---------|
| `X_train` | `(N+1 × p × L_train)` | required |
| `U_train` | `(N × q × L_train)` | required |
| `X_val` | `(N+1 × p × L_val)` | required |
| `U_val` | `(N × q × L_val)` | required |
| `'LambdaGrid'` | vector | `logspace(-3, 15, 50)` |
| `'Algorithm'` | char | `'cosmic'` |

**Outputs:**

| Field | Type | Description |
|-------|------|-------------|
| `bestResult` | struct | `sidLTVdisc` result at optimal λ |
| `bestLambda` | scalar | Optimal λ value |
| `allLosses` | `(n_grid × 1)` | Prediction loss at each λ |

### 8.5 Output Struct

| Field | Type | Description |
|-------|------|-------------|
| `A` | `(p × p × N)` | Time-varying dynamics matrices A(0), ..., A(N-1) |
| `B` | `(p × q × N)` | Time-varying input matrices B(0), ..., B(N-1) |
| `AStd` | `(p × p × N)` | Standard deviation of A(k) elements (requires uncertainty) |
| `BStd` | `(p × q × N)` | Standard deviation of B(k) elements (requires uncertainty) |
| `Covariance` | `(p+q × p+q × N)` | Posterior covariance Σ_kk per step (requires uncertainty) |
| `NoiseVariance` | scalar | Estimated σ̂² (requires uncertainty) |
| `Lambda` | scalar or `(N-1 × 1)` | Regularization values used |
| `Cost` | `(1 × 3)` | `[total, data_fidelity, regularization]` |
| `DataLength` | scalar | N (number of time steps) |
| `StateDim` | scalar | p |
| `InputDim` | scalar | q |
| `NumTrajectories` | scalar | L |
| `Algorithm` | char | `'cosmic'` |
| `Preconditioned` | logical | Whether preconditioning was applied |
| `Method` | char | `'sidLTVdisc'` |

### 8.6 Usage Examples

```matlab
% Basic identification with automatic lambda selection
result = sidLTVdisc(X, U, 'Lambda', 'auto');

% Manual lambda, scalar (uniform)
result = sidLTVdisc(X, U, 'Lambda', 1e5);

% Per-step lambda (e.g., lower near a known transient)
lambdaVec = 1e5 * ones(N-1, 1);
lambdaVec(50:60) = 1e2;    % allow more variation during transient
result = sidLTVdisc(X, U, 'Lambda', lambdaVec);

% With preconditioning for ill-conditioned data
result = sidLTVdisc(X, U, 'Lambda', 1e5, 'Precondition', true);

% Validation-based tuning
[best, bestLam, losses] = sidLTVdiscTune(X_train, U_train, X_val, U_val);
semilogx(logspace(-3,15,50), losses); xlabel('\lambda'); ylabel('RMSE');
```

### 8.7 Relationship to `sidFreqMap`

`sidFreqMap` and `sidLTVdisc` answer the same question from different perspectives:

| Aspect | `sidFreqMap` | `sidLTVdisc` |
|--------|---------------|-------------|
| Domain | Frequency × time | Time (state-space) |
| Model type | Non-parametric G(ω,t) | Parametric A(k), B(k) |
| Requires | Input-output data | State measurements |
| State dimension | Not needed | Must be known/chosen |
| Output | Transfer function estimate | Explicit state-space matrices |
| Use case | Diagnosis: *is* the system changing? | Modeling: *what* are the matrices? |
| Downstream | Visual analysis, coherence checking | Controller design (LTV LQR, MPC) |

A recommended workflow:

1. Run `sidSpectrogram` on `u` and `y` to understand signal characteristics.
2. Run `sidFreqMap` to diagnose whether and where the system is time-varying. When multiple trajectories are available, pass all of them — the ensemble-averaged spectral estimates will be more reliable than any single trajectory.
3. Run `sidLTVdisc` to obtain the explicit state-space model for controller design.
4. Validate: propagate the `sidLTVdisc` model and compare predicted states to measured states.

### 8.8 Variable-Length Trajectories

**Reference:** `docs/cosmic_uncertainty_derivation.md` §1.

When trajectories have different horizons, let `L(k) ⊆ {1,...,L}` be the set of trajectories active at time step `k`. The data matrices become:

```
D(k) = [X_{L(k)}(k)^T  U_{L(k)}(k)^T] / sqrt(|L(k)|)
```

Only the `S_kk` and `Θ_k` terms change; the regularization term `F^T Υ F` is unchanged because it couples only consecutive `C(k)` values and does not reference the data. The forward-backward pass structure is completely preserved.

**API change:** `X` and `U` accept cell arrays:

```matlab
X = {X1, X2, X3};   % X1 is (N1+1 x p), X2 is (N2+1 x p), etc.
U = {U1, U2, U3};   % U1 is (N1 x q), etc.
```

The total horizon `N` is `max(N1, N2, ..., N_L)`. Time steps with fewer active trajectories receive more regularization influence, which is the correct behavior.

### 8.9 Bayesian Uncertainty Estimation

**Reference:** `docs/cosmic_uncertainty_derivation.md` §2–4.

#### 8.9.1 Bayesian Interpretation

Under Gaussian noise `w(k) ~ N(0, σ² I)` on the state measurements, the COSMIC cost function is the negative log-posterior of a Bayesian model:

- **Likelihood:** `p(X' | C) ∝ exp(-h(C) / σ²)` — the data fidelity term.
- **Prior:** `p(C) ∝ exp(-g(C) / σ²)` — the smoothness regularizer is a Gaussian prior on consecutive differences of `C(k)` with precision `λ_k / σ²`.

The posterior is Gaussian:

```
p(C | X') = N(C*, H⁻¹ σ²)
```

where `C*` is the COSMIC solution (the MAP estimate) and `H` is the Hessian:

```
H = V^T V + F^T Υ F
```

This is exactly the block tridiagonal matrix `LM` from the COSMIC derivation. The posterior covariance is `Σ = σ² H⁻¹`.

#### 8.9.2 Diagonal Block Extraction via Forward-Backward Pass

The full `H⁻¹` is `N(p+q) × N(p+q)` — too large to store. But we only need the diagonal blocks `Σ_kk = σ² [H⁻¹]_kk`, which give the marginal posterior covariance of `C(k)` at each time step.

The diagonal blocks of a block tridiagonal inverse can be computed by a second backward pass reusing the `Λ_k` matrices from COSMIC's forward pass.

**Algorithm (Uncertainty Backward Pass):**

```
// Λ_k already computed during COSMIC forward pass

// Initialize at last time step
P(N-1) = Λ_{N-1}⁻¹

// Backward pass: k = N-2, ..., 0
For k = N-2 down to 0:
    G_k = λ_{k+1} Λ_k⁻¹                      // gain matrix
    P(k) = Λ_k⁻¹ + G_k P(k+1) G_k^T          // Joseph form
```

where `P(k) = [H⁻¹]_kk` is the `(p+q) × (p+q)` diagonal block of the inverse Hessian at step `k`.

**Complexity:** `O(N(p+q)³)` — identical to COSMIC itself. The `Λ_k⁻¹` are already computed during the forward pass, so the marginal cost is one additional backward sweep of matrix multiplications.

**Connection to Kalman smoothing:** The forward pass computes `Λ_k` (analogous to the Kalman filter's predicted covariance), and the uncertainty backward pass computes `P(k)` (analogous to the Rauch-Tung-Striebel smoother's smoothed covariance). This is not a coincidence — the Bayesian interpretation of COSMIC's regularized least squares *is* a Kalman smoother applied to the parameter evolution model `C(k+1) = C(k) + w_k`.

#### 8.9.3 Noise Variance Estimation

The noise variance `σ²` can be estimated from the data fidelity residuals:

```
σ̂² = (2 / (N × L × p)) × h(C*)
```

where `h(C*)` is the data fidelity term evaluated at the optimal solution. This is the maximum likelihood estimate under the Gaussian assumption.

#### 8.9.4 Output Fields

| Field | Type | Description |
|-------|------|-------------|
| `AStd` | `(p × p × N)` | Standard deviation of each A(k) element |
| `BStd` | `(p × q × N)` | Standard deviation of each B(k) element |
| `Covariance` | `(p+q × p+q × N)` | Posterior covariance `Σ_kk` at each step |
| `NoiseVariance` | scalar | Estimated `σ̂²` |

The standard deviations are extracted from the diagonal of `Σ_kk`:

```
AStd(i, j, k) = σ̂ × sqrt(P(k)_{j, j})    for the (i,j) element of A(k)
BStd(i, j, k) = σ̂ × sqrt(P(k)_{p+j, p+j}) for the (i,j) element of B(k)
```

(Note: `C(k) = [A(k)'; B(k)']`, so the rows of `C` are columns of `A` and `B`.)

### 8.10 Online/Recursive COSMIC

**Reference:** `docs/cosmic_online_recursion.md`.

#### 8.10.1 The Insight: Forward Pass Is Naturally Causal

COSMIC's forward pass computes `Λ_k` and `Y_k` sequentially — step `k` depends only on steps `0..k`. This means the forward pass can run in real time as data arrives. At any point, the "filtered" estimate `Y_k` is available as a causal estimate of `C(k)`, analogous to the Kalman filter's filtered state.

The backward pass touches all time steps and is non-causal — it requires the full trajectory. However, under the Bayesian/Kalman interpretation, the relationship between forward-only and full solution is precise:

| | Forward only (`Y_k`) | Full solution (`C(k)`) |
|---|---|---|
| Kalman analogy | Filtered estimate | Smoothed estimate |
| Uses data from | `0..k` | `0..N-1` |
| Uncertainty | Larger (`Λ_k⁻¹`) | Smaller (`P(k)`) |
| Available | Causally (real-time) | After full trajectory |

#### 8.10.2 Three Operating Modes

**Mode 1: Batch (existing).** Process full trajectory, forward + backward. Best accuracy. Use when all data is available.

**Mode 2: Filtered (real-time).** Run forward pass only. At each new time step `k`, compute `Λ_k` and `Y_k` from the new data `D(k)`, `X'(k)` and the previous `Λ_{k-1}`, `Y_{k-1}`. The estimate `Y_k` is immediately available. Uncertainty is `Λ_k⁻¹` (larger than smoothed, but honest about the causal constraint).

```
// When new measurement arrives at step k:
D_k = [x(k)^T  u(k)^T] / sqrt(L)
X'_k = x(k+1)^T / sqrt(L)
S_kk = D_k^T D_k + (λ_k + λ_{k+1}) I
Θ_k  = D_k^T X'_k

Λ_k = S_kk - λ_k² Λ_{k-1}⁻¹
Y_k = Λ_k⁻¹ (Θ_k + λ_k Y_{k-1})

// Extract filtered estimate:
A_filtered(k) = Y_k(1:p, :)'
B_filtered(k) = Y_k(p+1:end, :)'
// Filtered uncertainty:
P_filtered(k) = Λ_k⁻¹
```

**Cost per step:** One `(p+q) × (p+q)` matrix inversion + one matrix multiply = `O((p+q)³)`. Constant time per step, independent of history length.

**Mode 3: Windowed smoother.** Maintain a sliding window of the last `W` time steps. At each new step:
1. Extend the forward pass by one step (Mode 2).
2. Run the backward pass over only the window `[k-W+1, ..., k]`, using the forward pass quantities `Λ`, `Y` already stored.
3. The smoothed estimates within the window are improved; older estimates are fixed.

This gives a practical middle ground: `O(W(p+q)³)` per step, with smoothed accuracy within the window. The boundary condition at `k-W` uses the filtered estimate, which introduces a small approximation that decays exponentially with `W` if `λ` provides sufficient coupling.

#### 8.10.3 API

```matlab
% Initialize the recursive estimator
rec = sidLTVdiscInit(p, q, 'Lambda', lambda);

% Process measurements one at a time
for k = 1:N
    rec = sidLTVdiscUpdate(rec, x(:,k), x(:,k+1), u(:,k));
    A_now = rec.A_filtered;  % immediately available
    P_now = rec.P_filtered;  % filtered uncertainty
end

% Optional: smooth over recent window
rec = sidLTVdiscSmooth(rec, 'Window', 50);
A_smoothed = rec.A_smoothed;  % improved estimates for last 50 steps
```

### 8.11 Lambda Tuning via Frequency Response

**Reference:** `docs/cosmic_uncertainty_derivation.md` §5.

#### 8.11.1 Concept

`sidFreqMap` produces a non-parametric estimate `Ĝ_BT(ω, t)` with uncertainty, independent of `λ`. For any candidate `λ`, compute the frozen transfer function from COSMIC's `A(k)`, `B(k)`:

```
G_cosmic(ω, k) = (e^{jω} I - A(k))⁻¹ B(k)
```

and propagate the posterior covariance `Σ_kk` to obtain `σ_cosmic(ω, k)` via the Jacobian of the `(A, B) → G(ω)` mapping.

The criterion: **find the largest λ whose COSMIC posterior bands are consistent with the non-parametric bands.**

**Multi-trajectory:** When multiple trajectories are available, `sidFreqMap` should be called with all `L` trajectories to produce ensemble-averaged estimates. This makes the variation metric `Δ_k` in the spectral pre-scan significantly more reliable, since the within-trajectory estimation noise averages out while genuine system variation is preserved. See `docs/multi_trajectory_spectral_theory.md` §3 and §7.

#### 8.11.2 Consistency Score

At each grid point `(ω_j, t_i)`:

```
d²(j,i) = |G_cosmic(ω_j, t_i) - Ĝ_BT(ω_j, t_i)|² / (σ²_cosmic(j,i) + σ²_BT(j,i))
```

This is a Mahalanobis-like distance. Under the null hypothesis (both estimators are estimating the same true G), `d²` is approximately χ² distributed.

Aggregate score:

```
S(λ) = (1 / n_grid) Σ_{j,i} 1[d²(j,i) < χ²_{0.95}]
```

i.e., the fraction of grid points where the two estimates are consistent at 95% level.

Select `λ* = max{λ : S(λ) > 0.90}` — the largest λ for which at least 90% of grid points are consistent.

#### 8.11.3 Depends On

- `sidFreqMap` (§6) for the non-parametric reference.
- Bayesian uncertainty (§8.9) for COSMIC posterior bands.
- `sidLTVdiscFrozen` utility for computing `G_cosmic(ω, k)`.

### 8.12 Output-Only Estimation (`sidLTVdiscIO`)

Two-stage approach for when only outputs `y(k) = C_obs x(k) + D_obs u(k)` are measured:

1. Use an initial LTI subspace method (or user-supplied observer) to estimate state trajectories `x̂(k)` from `y(k)` and `u(k)`.
2. Feed `x̂(k)` into `sidLTVdisc` as if they were true states.

The user is warned that state estimates carry observer error.

### 8.13 Deferred Extensions

The following are out of scope for v1.0:

- **Alternative algorithms:** TVERA, TVOKID, LTVModels (the `'Algorithm'` parameter is ready for this).
- **Alternative regularization norms:** Non-squared L2, L1 (total variation).
- **EM-style output estimation:** Alternating state reconstruction and COSMIC.
- **Direct output equation formulation:** Joint estimation of dynamics and output matrices.
- **GCV lambda selection.**

---

## 9. Output Struct

All `sidFreq*` functions return a struct with these fields:

| Field | Type | Description |
|-------|------|-------------|
| `Frequency` | `(n_f × 1)` real | Frequency vector in rad/sample |
| `FrequencyHz` | `(n_f × 1)` real | Frequency vector in Hz: `ω / (2π Ts)` |
| `Response` | `(n_f × n_y × n_u)` complex | Frequency response `Ĝ(ω)` |
| `ResponseStd` | `(n_f × n_y × n_u)` real | Standard deviation of `Ĝ` |
| `NoiseSpectrum` | `(n_f × n_y × n_y)` real | Noise spectrum `Φ̂_v(ω)` or `Φ̂_y(ω)` |
| `NoiseSpectrumStd` | `(n_f × n_y × n_y)` real | Standard deviation of noise spectrum |
| `Coherence` | `(n_f × 1)` real | Squared coherence `γ̂²(ω)` (SISO only, `[]` for MIMO) |
| `SampleTime` | scalar | Sample time `Ts` in seconds |
| `WindowSize` | scalar or vector | Window size `M` (scalar for BT, vector for BTFDR) |
| `DataLength` | scalar | Number of samples `N` |
| `NumTrajectories` | scalar | Number of trajectories `L` used in estimation |
| `Method` | char | `'sidFreqBT'`, `'sidFreqBTFDR'`, or `'sidFreqETFE'` |

**Dimension conventions:**
- SISO: `Response` is `(n_f × 1)`, `NoiseSpectrum` is `(n_f × 1)`.
- MIMO: Dimensions are `(n_f × n_y × n_u)` for `Response` and `(n_f × n_y × n_y)` for `NoiseSpectrum`.

**Time series mode:** `Response` and `ResponseStd` are empty (`[]`). `Coherence` is empty. `NoiseSpectrum` contains `Φ̂_y(ω)`.

---

## 10. Edge Cases and Validation

### 10.1 Input Validation

| Condition | Action |
|-----------|--------|
| `N < 2 × M` | Reduce `M` to `floor(N/2)` and issue warning |
| `M < 2` | Error: window size must be at least 2 |
| `size(y,1) ~= size(u,1)` | Error: input and output must have same number of samples |
| `N < 10` | Warning: very short data, estimates will be unreliable |
| `y` or `u` contains NaN or Inf | Error: data must be finite |
| `y` or `u` is not real | Error: complex data not supported in v1.0 |
| Any frequency `ω_k ≤ 0` or `ω_k > π` | Error: frequencies must be in (0, π] rad/sample |
| `Ts ≤ 0` | Error: sample time must be positive |

### 10.2 Numerical Edge Cases

| Condition | Action |
|-----------|--------|
| `Φ̂_u(ω_k) ≈ 0` | Set `Ĝ(ω_k) = NaN`, `σ_G(ω_k) = Inf`, issue warning |
| `Φ̂_v(ω_k) < 0` | Clamp to 0 |
| `γ̂²(ω_k) > 1` (numerical error) | Clamp to 1 |
| `γ̂²(ω_k) < 0` (numerical error) | Clamp to 0 |

### 10.3 Degenerate Inputs

| Condition | Action |
|-----------|--------|
| `u` is constant (zero variance) | Same as `Φ̂_u ≈ 0` at all frequencies; `Ĝ = NaN` everywhere, with warning |
| `y` is constant | Valid; `Φ̂_y ≈ 0` at all frequencies |
| `u = y` (perfect coherence) | Valid; `γ̂² ≈ 1`, `Φ̂_v ≈ 0`, very small `σ_G` |

---

## 11. Plotting

### 11.1 `sidBodePlot`

Produces a two-panel figure:
- **Top panel:** Magnitude `20 × log10(|Ĝ(ω)|)` in dB vs. frequency
- **Bottom panel:** Phase `angle(Ĝ(ω)) × 180/π` in degrees vs. frequency

Both panels use logarithmic frequency axis (rad/s by default, Hz if requested).

Confidence bands are shown as a shaded region at `±p` standard deviations (default `p = 3`):
- Magnitude band: `20 × log10(|Ĝ| ± p × σ_G)` — note this is applied to the linear magnitude, then converted to dB.
- Phase band: `±p × σ_G / |Ĝ| × 180/π` — small-angle approximation for phase uncertainty.

### 11.2 `sidSpectrumPlot`

Single panel: `10 × log10(Φ̂_v(ω))` in dB vs. frequency (log axis).

Confidence band: `10 × log10(Φ̂_v ± p × σ_Φv)` — applied in linear scale, converted to dB.

### 11.3 Options

Both plotting functions accept name-value options:

| Option | Default | Description |
|--------|---------|-------------|
| `'Confidence'` | `3` | Number of standard deviations for shaded band |
| `'FrequencyUnit'` | `'rad/s'` | `'rad/s'` or `'Hz'` |
| `'ShowConfidence'` | `true` | Whether to show the confidence band |
| `'Color'` | MATLAB default | Line color |
| `'LineWidth'` | `1.5` | Line width |
| `'Axes'` | `[]` | Axes handle (creates new figure if empty) |

---

## 12. References

1. Ljung, L. *System Identification: Theory for the User*, 2nd ed. Prentice Hall, 1999.
   - §2.3: Spectral analysis fundamentals
   - §6.3–6.4: Non-parametric frequency-domain methods
   - Table 6.1: Default window sizes
   - p. 184: Asymptotic variance of frequency response estimate
   - p. 188: Asymptotic variance of spectral estimate

2. Blackman, R.B. and Tukey, J.W. *The Measurement of Power Spectra*. Dover, 1959.

3. Kay, S.M. *Modern Spectral Estimation: Theory and Application*. Prentice Hall, 1988.

4. Stoica, P. and Moses, R.L. *Spectral Analysis of Signals*. Prentice Hall, 2005.

5. Carvalho, M., Soares, C., Lourenço, P., and Ventura, R. "COSMIC: fast closed-form identification from large-scale data for LTV systems." arXiv:2112.04355, 2022.

6. Łaszkiewicz, P., Carvalho, M., Soares, C., and Lourenço, P. "The impact of modeling approaches on controlling safety-critical, highly perturbed systems: the case for data-driven models." arXiv:2509.13531, 2025.

7. Carlson, F.B., Robertsson, A., and Johansson, R. "Identification of LTV dynamical models with smooth or discontinuous time evolution by means of convex optimization." IEEE ICCA, 2018.

8. Majji, M., Juang, J.-N., and Junkins, J.L. "Time-varying eigensystem realization algorithm." JGCD 33(1), 2010.

9. Majji, M., Juang, J.-N., and Junkins, J.L. "Observer/Kalman-filter time-varying system identification." JGCD 33(3), 2010.

10. Bendat, J.S. and Piersol, A.G. *Random Data: Analysis and Measurement Procedures*, 4th ed. Wiley, 2010. (Ch. 9: Statistical errors in spectral estimates; Ch. 11: Multiple-input/output relationships.)

11. Antoni, J. and Schoukens, J. "A comprehensive study of the bias and variance of frequency-response-function measurements: optimal window selection and overlapping strategies." Automatica, 43(10):1723–1736, 2007.
