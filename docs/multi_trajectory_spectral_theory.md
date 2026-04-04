# Multi-Trajectory Spectral Analysis for System Identification

## 1. Motivation

COSMIC (Closed-form Optimal data-driven linear time-varying System IdentifiCation) identifies discrete LTV systems from multiple trajectories drawn from the same time-varying dynamics. At each time step k, the data matrix D(k) aggregates state-input pairs from all L trajectories:

```
D(k) = [x_1(k)' u_1(k)'; x_2(k)' u_2(k)'; ...; x_L(k)' u_L(k)'] / sqrt(L)
```

This ensemble averaging is fundamental to COSMIC's statistical strength — more trajectories reduce variance without requiring longer time horizons.

The frequency-domain functions in `sid` (`sidFreqBT`, `sidFreqMap`, `sidSpectrogram`) currently accept only a single input-output record. This creates an asymmetry: the time-domain side of the toolbox can exploit multiple realizations, but the frequency-domain side cannot. Extending the spectral functions to accept multiple trajectories would:

1. Improve spectral estimates by ensemble averaging across independent realizations.
2. Enable direct comparison between multi-trajectory COSMIC results and multi-trajectory spectral diagnostics (e.g., for lambda tuning via `sidFreqMap`).
3. Align the `sid` API so that the same multi-trajectory dataset can be passed to any function without reshaping.


## 2. The Ensemble Averaging Principle

### 2.1 Setup

Consider L trajectories of an LTI system (or locally LTI within a segment):

```
y_l(t) = G(q) u_l(t) + v_l(t)       l = 1, ..., L;  t = 1, ..., N
```

Each trajectory has different input `u_l(t)` and noise `v_l(t)`, but the same system `G(q)`. The inputs and noise are assumed independent across trajectories.

### 2.2 Covariance Averaging (Blackman-Tukey)

For the BT method, the per-trajectory biased covariance estimates are:

```
R̂_yu^(l)(τ) = (1/N) Σ_{t=1}^{N-|τ|} y_l(t+|τ|) u_l(t)
```

The multi-trajectory covariance is the ensemble average:

```
R̂_yu^ens(τ) = (1/L) Σ_{l=1}^{L} R̂_yu^(l)(τ)
```

This averaged covariance is then windowed and Fourier-transformed as usual:

```
Φ̂_yu^ens(ω) = Σ_{τ=-M}^{M} R̂_yu^ens(τ) × W_M(τ) × exp(-jωτ)
```

**Why this works:** Each R̂_yu^(l)(τ) is an unbiased estimate of the true cross-covariance R_yu(τ) (up to the 1/N vs 1/(N-|τ|) normalization). Since the trajectories are independent, their estimation errors are independent and zero-mean. Averaging L of them reduces the variance of R̂_yu^ens(τ) by a factor of L, without affecting bias or frequency resolution.

**Key property:** Ensemble averaging reduces variance without sacrificing frequency resolution. This is fundamentally different from increasing the BT window size M (which reduces variance by reducing resolution) or increasing the data length N within one trajectory (which improves both, but requires a longer experiment).

### 2.3 Periodogram Averaging (Welch / ETFE)

For the Welch method or ETFE, the per-trajectory cross-periodograms are:

```
Ŝ_yu^(l)(ω) = Y_l(ω) conj(U_l(ω))
```

The multi-trajectory average is:

```
Ŝ_yu^ens(ω) = (1/L) Σ_{l=1}^{L} Ŝ_yu^(l)(ω)
```

This is equivalent to having L additional sub-segments in the Welch averaging — each trajectory is an independent realization, so the averaged periodogram has L times the effective degrees of freedom:

```
ν_ens = L × ν_single
```

where ν_single is the degrees of freedom from the within-trajectory Welch averaging. The total effective degrees of freedom combine both sources of averaging.

### 2.4 The Result

For any spectral estimator, the multi-trajectory estimates are:

```
Φ̂_yu^ens(ω) = (1/L) Σ_l Φ̂_yu^(l)(ω)
Φ̂_u^ens(ω)  = (1/L) Σ_l Φ̂_u^(l)(ω)
Φ̂_y^ens(ω)  = (1/L) Σ_l Φ̂_y^(l)(ω)
```

The transfer function and coherence are then formed from the averaged spectra:

```
Ĝ^ens(ω) = Φ̂_yu^ens(ω) / Φ̂_u^ens(ω)
γ̂²^ens(ω) = |Φ̂_yu^ens(ω)|² / (Φ̂_y^ens(ω) × Φ̂_u^ens(ω))
```

**Important:** The averaging happens at the cross-spectral level, before forming ratios. This is equivalent to the H1 estimator applied to the concatenated dataset. Averaging the transfer functions directly (i.e., averaging Ĝ^(l)(ω) across l) would give a different and generally inferior result, because the ratio of averages ≠ the average of ratios.


## 3. Extension to Time-Varying Systems

### 3.1 Multi-Trajectory `sidFreqMap`

For an LTV system, each trajectory l provides an independent realization of the same time-varying dynamics:

```
y_l(t) = G(q, t) u_l(t) + v_l(t)
```

Within each time segment k of `sidFreqMap`, we have L independent input-output records of length L_seg. The ensemble-averaged spectral estimates within the segment are:

```
Φ̂_yu^ens(ω, t_k) = (1/L) Σ_l Φ̂_yu^(l)(ω, t_k)
```

where Φ̂_yu^(l)(ω, t_k) is the BT or Welch estimate from trajectory l's data within segment k.

This directly parallels how COSMIC uses multiple trajectories: at each time step k, COSMIC aggregates state-input pairs across trajectories; at each time segment k, `sidFreqMap` aggregates spectral estimates across trajectories.

### 3.2 Variance Reduction

The variance of the ensemble-averaged spectral estimate within each segment is:

**BT:**
```
Var{Φ̂^ens(ω)} ≈ (1/L) × (2 C_W / N_seg) × Φ(ω)²
```

**Welch:**
```
Var{Φ̂^ens(ω)} ≈ Φ(ω)² / (L × ν_single)
```

In both cases, L additional trajectories are equivalent to multiplying the effective data length by L, but without the time-resolution cost of using a longer segment.

### 3.3 Uncertainty Formulas

The asymptotic variance of the multi-trajectory transfer function estimate becomes:

```
Var{Ĝ^ens(ω)} ≈ (C_W / (L × N_seg)) × |G(ω)|² × (1 - γ²(ω)) / γ²(ω)
```

Note the factor of L in the denominator — this is the only change from the single-trajectory formula. The coherence γ² is now the ensemble coherence, which is generally higher than any single-trajectory coherence because the noise averages out while the signal accumulates.

### 3.4 Different Inputs Across Trajectories

An important practical point: the L trajectories need not use the same input signal. In fact, using different (ideally uncorrelated) inputs across trajectories is beneficial:

- It broadens the frequency content of the aggregate excitation.
- It breaks any correlation between input and disturbance that might exist within a single trajectory.
- For MIMO systems, trajectories with different input directions improve the conditioning of Φ̂_u^ens(ω).

The only requirement is that the *system* is the same across trajectories (same G(q, t) for LTV, same G(q) for LTI). The inputs and noise realizations can differ freely.


## 4. Variable-Length Trajectories

When trajectories have different lengths N_l, not all trajectories contribute to every time segment in `sidFreqMap`. Let L(k) be the set of trajectories that span segment k:

```
L(k) = { l : trajectory l has data covering segment k }
```

The ensemble average within segment k uses only the active trajectories:

```
Φ̂_yu^ens(ω, t_k) = (1/|L(k)|) Σ_{l ∈ L(k)} Φ̂_yu^(l)(ω, t_k)
```

This mirrors COSMIC's variable-length trajectory handling (§8.8 of the SPEC): time segments with fewer active trajectories produce noisier estimates, which is correctly reflected in the uncertainty.

The uncertainty at segment k scales as:

```
Var{Ĝ^ens(ω, t_k)} ∝ 1 / |L(k)|
```

Segments covered by all L trajectories have the lowest variance; segments at the extremes (covered by only the longest trajectories) have higher variance.


## 5. Multi-Trajectory `sidSpectrogram`

For spectrograms of a single signal (no system identification), ensemble averaging across trajectories computes the averaged periodogram:

```
P̂^ens(ω, t_k) = (1/L) Σ_l P̂^(l)(ω, t_k)
```

This is the standard approach in fields like neuroscience (event-related spectral perturbation / ERSP), where multiple trials of the same stimulus are averaged in the time-frequency domain to separate evoked responses from background noise.

In the `sid` context, this is useful for:
- Comparing the averaged spectrogram of u and y alongside the averaged `sidFreqMap` output, with consistent variance reduction across all three views.
- Reducing noise in the spectral pre-scan for lambda tuning (§4 of the tuning document), where averaged variation metrics ∆_k are more reliable than single-trajectory estimates.


## 6. Multi-Trajectory `sidFreqBT` (LTI Case)

For a pure LTI system, all time segments see the same transfer function. Multiple trajectories simply provide more data for the same estimation problem. The ensemble-averaged covariance is:

```
R̂_yu^ens(τ) = (1/L) Σ_l R̂_yu^(l)(τ)
```

This is computationally equivalent to concatenating all L trajectories into one long record — except that concatenation would introduce discontinuities at the boundaries, while ensemble averaging avoids this.

The variance reduction is a factor of L, and the implementation is trivial: compute per-trajectory covariances, average, then proceed with the standard BT pipeline.


## 7. Connection to COSMIC Multi-Trajectory

The multi-trajectory extension of the spectral functions creates a unified data flow:

```
                 Multi-trajectory data
                 (y_l, u_l for l = 1..L)
                        │
           ┌────────────┼────────────┐
           │            │            │
    sidFreqMap     sidSpectrogram  sidLTVdisc
    (ensemble-     (ensemble-     (COSMIC,
     averaged       averaged      D(k) from
     spectra)       power)        all L)
           │            │            │
           └────────────┼────────────┘
                        │
              Consistent multi-trajectory
              estimates across all domains
```

All three functions aggregate information from L trajectories at each time step/segment, using the same data and consistent conventions. This enables:

- **Lambda tuning:** The spectral pre-scan (§8d of SPEC) uses ensemble-averaged sidFreqMap output, making the variation metric ∆_k much more reliable than single-trajectory estimates.
- **Cross-validation:** Comparing ensemble-averaged frequency response from sidFreqMap against the frozen transfer function from COSMIC, both using the same multi-trajectory data.
- **Diagnostics:** If the ensemble-averaged coherence is high but a single trajectory's coherence is low, the problem is poor excitation in that trajectory, not system nonlinearity.


## 8. API Convention

All `sid` functions that accept multi-trajectory data should use a consistent convention:

**3D arrays** (when all trajectories have the same length):
```matlab
y = (N x n_y x L)       % L trajectories, each N samples, n_y channels
u = (N x n_u x L)       % L trajectories, each N samples, n_u channels
```

**Cell arrays** (when trajectories have different lengths):
```matlab
y = {y1, y2, ..., yL}   % y_l is (N_l x n_y)
u = {u1, u2, ..., uL}   % u_l is (N_l x n_u)
```

This matches the existing COSMIC convention in `sidLTVdisc` (§8.2 and §8.8 of SPEC).

For `sidFreqBT` (LTI), the function accepts multi-trajectory input and ensemble-averages covariances. For `sidFreqMap` and `sidSpectrogram`, the function accepts multi-trajectory input and ensemble-averages spectral estimates within each segment.

The output struct is unchanged — it always contains the final ensemble-averaged estimates. A new field `NumTrajectories` records L.


## 9. Summary

| Function | Single-trajectory (existing) | Multi-trajectory (extension) |
|----------|------------------------------|------------------------------|
| `sidFreqBT` | Covariance from one record | Average covariances across L records |
| `sidFreqETFE` | Periodogram from one record | Average periodograms across L records |
| `sidFreqMap` | Per-segment BT/Welch from one record | Per-segment ensemble-averaged spectra |
| `sidSpectrogram` | Per-segment PSD from one record | Per-segment averaged PSD (ERSP-like) |
| `sidLTVdisc` | Single trajectory (already supports L) | Already supports L trajectories |

The mathematical operation is the same in every case: average the spectral/covariance estimates across trajectories before forming ratios. This is the correct statistical procedure because cross-spectral averaging preserves the H1 estimator structure (ratio of averages, not average of ratios).

The variance reduction is uniformly a factor of L, and the frequency resolution is unchanged. This is the fundamental advantage of ensemble averaging over any within-trajectory technique for improving spectral estimates.


## References

1. Bendat, J.S. and Piersol, A.G. *Random Data: Analysis and Measurement Procedures*, 4th ed. Wiley, 2010. (Ch. 9: Statistical errors in spectral estimates; Ch. 11: Multiple-input/output relationships.)

2. Ljung, L. *System Identification: Theory for the User*, 2nd ed. Prentice Hall, 1999. (§2.3, §6.3–6.4.)

3. Carvalho, M. et al. "COSMIC: fast closed-form identification from large-scale data for LTV systems." arXiv:2112.04355, 2022.

4. Makeig, S. "Auditory event-related dynamics of the EEG spectrum and effects of exposure to tones." Electroencephalography and Clinical Neurophysiology, 86(4):283–293, 1993. (Event-related spectral perturbation / ensemble averaging in time-frequency analysis.)

5. Antoni, J. and Schoukens, J. "A comprehensive study of the bias and variance of frequency-response-function measurements: optimal window selection and overlapping strategies." Automatica, 43(10):1723–1736, 2007.
