# Open-Source Implementation of MATLAB's `spa` Function

## Algorithm Review & Implementation Plan

---

## 1. What `spa` Does

The `spa` function in the MATLAB System Identification Toolbox performs **non-parametric frequency-response estimation** using the **Blackman-Tukey (correlogram) method**. Given time-domain input/output data from a linear system, it produces three things:

1. **A frequency-response function estimate** Ĝ(e^{jω}) — the estimated transfer function evaluated at discrete frequencies on the unit circle.
2. **A noise spectrum estimate** Φ̂_v(ω) — the power spectral density of the output disturbance.
3. **Uncertainty estimates** (standard deviations) for both quantities, enabling confidence interval plots.

The underlying system model is:

```
y(t) = G(q) u(t) + v(t)
```

where `u(t)` is the input, `y(t)` is the output, `G(q)` is the transfer function, and `v(t)` is output noise assumed independent of `u(t)`.

---

## 2. Algorithm Breakdown

The algorithm proceeds in three stages. Every step references Ljung, *System Identification: Theory for the User*, 2nd ed. (1999), which is the canonical source.

### Stage 1: Covariance Estimation

Compute the sample auto-covariance and cross-covariance functions from the N data points:

```
R̂_y(τ)  = (1/N) Σ_{t=1}^{N-τ} y(t+τ) y(t)        (output auto-covariance)
R̂_u(τ)  = (1/N) Σ_{t=1}^{N-τ} u(t+τ) u(t)        (input auto-covariance)
R̂_yu(τ) = (1/N) Σ_{t=1}^{N-τ} y(t+τ) u(t)        (cross-covariance)
```

These are computed for lags τ = -M, ..., 0, ..., M where M is the **window size** (lag size).

**Key details:**
- The `1/N` normalization (biased estimator) is used, not `1/(N-|τ|)`. The biased estimator guarantees a non-negative spectral estimate and has lower mean-squared error.
- In MATLAB's implementation, the internal function `covf` handles this step efficiently.
- For MIMO systems, these become matrix-valued functions: R̂_yu(τ) is N_y × N_u.

### Stage 2: Windowed Fourier Transform (Spectral Estimation)

The covariance estimates are multiplied by a **Hann (Hanning) window** W_M(τ) of width M, then Fourier-transformed:

```
Φ̂_y(ω)  = Σ_{τ=-M}^{M}  R̂_y(τ)  · W_M(τ) · e^{-jωτ}
Φ̂_u(ω)  = Σ_{τ=-M}^{M}  R̂_u(τ)  · W_M(τ) · e^{-jωτ}
Φ̂_yu(ω) = Σ_{τ=-M}^{M}  R̂_yu(τ) · W_M(τ) · e^{-jωτ}
```

**The Hann window** is defined as:

```
W_M(τ) = 0.5 (1 + cos(π τ / M))    for |τ| ≤ M
W_M(τ) = 0                           for |τ| > M
```

Note: An older version of the documentation references a **Hamming** window; the current version uses **Hann**. The difference is minor (Hamming: 0.54 + 0.46 cos(...); Hann: 0.5 + 0.5 cos(...)). An implementation should default to Hann but could offer Hamming as an option.

**Window size M and frequency resolution:**
- Default: `M = min(N/10, 30)` where N is the data length.
- The frequency resolution is approximately `2π/M` rad/sample.
- Larger M → finer resolution but higher variance (the classic bias-variance tradeoff).

**Frequency grid:**
- Default: 128 equally spaced values in (0, π] rad/sample: `ω_k = k·π/128` for k = 1..128.
- When using the default grid, the Fourier transforms can be computed via **FFT** by zero-padding the windowed covariance sequence to length 128 (or larger power of 2), which is much faster than the direct DFT summation.
- When the user specifies custom frequencies, the summation must be computed directly for each ω.

### Stage 3: Forming the Estimates

**Frequency response:**
```
Ĝ(e^{jω}) = Φ̂_yu(ω) / Φ̂_u(ω)
```

**Noise spectrum** (derived from the assumption that u and v are independent):
```
Φ̂_v(ω) = Φ̂_y(ω) − |Φ̂_yu(ω)|² / Φ̂_u(ω)
```

This is the coherent subtraction formula. It removes the portion of the output spectrum explained by the input, leaving the unexplained noise.

**For time-series data** (no input), only the output power spectrum Φ̂_y(ω) is computed — stages involving input-related quantities are skipped.

### Uncertainty Estimation

The standard deviations of the estimates are computed using asymptotic variance formulas from Ljung (1999), pages 184 and 188. The key results are:

**For the spectral estimate Φ̂_y(ω):**
```
Var{Φ̂_y(ω)} ≈ (2/N) · Φ_y²(ω) · Σ_{τ=-M}^{M} W_M²(τ)
```

The sum `Σ W_M²(τ)` is a constant that depends only on the window and M. For the Hann window, this equals approximately `(3/8)·(2M+1)`.

**For the frequency response Ĝ(e^{jω}):**
The variance depends on the coherence function γ²(ω) between input and output:

```
γ²(ω) = |Φ_yu(ω)|² / (Φ_y(ω) · Φ_u(ω))
```

```
Var{Ĝ(e^{jω})} ≈ (1/N) · |G(e^{jω})|² · (1 - γ²(ω)) / γ²(ω) · Σ W_M²(τ)
```

In practice, the estimated coherence and transfer function are substituted for their true values.

---

## 3. MIMO Extension

For multi-input multi-output systems, the covariances become matrices:

```
R̂_z(τ) = (1/N) Σ_{t=1}^{N-τ} z(t+τ) z(t)ᴴ
```

where `z(t) = [y(t); u(t)]` is the stacked output-input vector. The full spectral matrix is:

```
S(ω) = Σ_{τ=-M}^{M} R̂_z(τ) · W_M(τ) · e^{-jωτ}
```

This N_z × N_z matrix (where N_z = N_y + N_u) is then partitioned to extract Φ̂_yy, Φ̂_uu, and Φ̂_yu sub-blocks. The frequency response becomes:

```
Ĝ(e^{jω}) = Φ̂_yu(ω) · Φ̂_uu(ω)⁻¹
```

---

## 4. Comparison with Related Methods

| Method | Description | Tradeoff |
|--------|-------------|----------|
| **spa (Blackman-Tukey)** | Windowed covariance → DFT. Fixed resolution. | Smooth estimates, easy uncertainty quantification. Resolution limited by M. |
| **etfe** | Ratio of output/input DFTs (periodogram). | Maximum resolution but very noisy. No built-in smoothing. |
| **spafdr** | Like spa but with frequency-dependent resolution. | Allows fine resolution at some frequencies and coarse at others. More flexible but more complex. |
| **Welch's method** | Segmented, overlapping periodograms averaged. | Widely available (scipy.signal.welch), but doesn't naturally produce a transfer function estimate with uncertainty. |

The key differentiator of `spa` versus what's available in scipy/numpy is the **complete system identification workflow**: transfer function estimation with uncertainty, not just power spectrum estimation.

---

## 5. Implementation Plan

### Target: Python library using NumPy/SciPy

### 5.1 Data Structures

```python
@dataclass
class SpaResult:
    """Result of spectral analysis."""
    frequencies: np.ndarray       # ω values (rad/sample or rad/s)
    response: np.ndarray          # Complex frequency response G(e^{jω}), shape (n_freq, n_out, n_in)
    response_std: np.ndarray      # Standard deviation of G, same shape
    noise_spectrum: np.ndarray    # Φ_v(ω), shape (n_freq, n_out, n_out)
    noise_spectrum_std: np.ndarray
    spectrum_matrix: np.ndarray   # Full spectral matrix S(ω), shape (n_freq, n_z, n_z)
    sample_time: float            # Ts
    window_size: int              # M
```

### 5.2 Core Function Signature

```python
def spa(
    y: np.ndarray,
    u: np.ndarray | None = None,
    window_size: int | None = None,
    frequencies: np.ndarray | None = None,
    sample_time: float = 1.0,
    n_frequencies: int = 128,
) -> SpaResult:
    """
    Estimate frequency response and noise spectrum via Blackman-Tukey spectral analysis.
    
    Parameters
    ----------
    y : array_like, shape (N,) or (N, n_out)
        Output time series.
    u : array_like, shape (N,) or (N, n_in), optional
        Input time series. If None, only the output spectrum is estimated.
    window_size : int, optional
        Hann window lag size M. Default: min(N//10, 30).
    frequencies : array_like, optional
        Frequencies (rad/sample) at which to evaluate. Default: n_frequencies
        log-spaced values in (0, π].
    sample_time : float
        Sample time in seconds. Default: 1.0.
    n_frequencies : int
        Number of default frequency points. Default: 128.
    
    Returns
    -------
    SpaResult
        Frequency response, noise spectrum, uncertainties, and full spectral matrix.
    """
```

### 5.3 Implementation Modules

**Module 1: `covariance.py` — Covariance Estimation**

```python
def estimate_covariance(x, y, max_lag):
    """Biased cross-covariance estimate R̂_xy(τ) for τ = 0..max_lag."""
    # Use np.correlate or FFT-based correlation for efficiency
    # Return shape: (max_lag+1, n_x, n_y) for MIMO
```

Key considerations:
- Use FFT-based correlation (`scipy.signal.fftconvolve` or manual FFT multiply) for large N, since direct summation is O(N·M) while FFT is O(N log N).
- Compute for both positive and negative lags, exploiting R̂_xy(-τ) = R̂_yx(τ)ᴴ.

**Module 2: `windows.py` — Window Functions**

```python
def hann_window(M):
    """Hann (Hanning) lag window of size M."""
    tau = np.arange(-M, M + 1)
    return 0.5 * (1 + np.cos(np.pi * tau / M))

def window_norm_squared(M, window='hann'):
    """Compute Σ W²(τ) for variance formulas."""
    w = hann_window(M)
    return np.sum(w**2)
```

**Module 3: `spectral.py` — Spectral Estimation Core**

```python
def windowed_spectrum(cov, window, frequencies=None):
    """
    Compute Φ̂(ω) = Σ R̂(τ) W(τ) e^{-jωτ}
    
    If frequencies is None, use FFT (fast path).
    Otherwise, compute DFT at specified frequencies (slow path).
    """
```

Two code paths:
1. **FFT path** (default frequencies): Construct the windowed covariance sequence, zero-pad, apply `np.fft.fft`, and extract the desired frequency bins.
2. **Direct DFT path** (custom frequencies): For each ω_k, compute the sum explicitly. Can be vectorized with outer products.

**Module 4: `uncertainty.py` — Variance Estimation**

```python
def transfer_function_variance(G_hat, coherence_sq, N, W_norm_sq):
    """Asymptotic variance of frequency response estimate."""
    return (W_norm_sq / N) * np.abs(G_hat)**2 * (1 - coherence_sq) / np.maximum(coherence_sq, 1e-10)

def spectrum_variance(Phi_hat, N, W_norm_sq):
    """Asymptotic variance of spectral estimate."""
    return (2 * W_norm_sq / N) * Phi_hat**2
```

**Module 5: `spa.py` — Main Entry Point**

Orchestrates the pipeline: validate inputs → compute covariances → apply window → FFT/DFT → form ratios → compute uncertainties → return `SpaResult`.

### 5.4 Visualization

```python
def bode_plot(result: SpaResult, confidence=3, ax=None):
    """Plot Bode diagram (magnitude + phase) with confidence bands."""

def spectrum_plot(result: SpaResult, confidence=3, ax=None):
    """Plot noise spectrum with confidence bands."""
```

### 5.5 Validation Strategy

1. **Unit tests against known analytic cases:**
   - White noise input → flat spectrum.
   - Known first-order system G(z) = 1/(1 - 0.9z⁻¹) → compare estimate to true Bode plot.
   - Pure sinusoid → spectral peak at correct frequency.

2. **Cross-validation against MATLAB:**
   - Generate test data in MATLAB, run `spa`, export results.
   - Run Python implementation on same data, compare numerically.
   - Target: agreement within 1e-6 for identical inputs and parameters.

3. **Edge cases:**
   - Very short data (N < 100).
   - Window size M close to N (should warn user).
   - Time-series mode (no input).
   - MIMO systems (2×2, 3×1, etc.).

### 5.6 Development Phases

| Phase | Deliverable | Effort |
|-------|-------------|--------|
| **Phase 1** | SISO time-series (output spectrum only) | 1–2 days |
| **Phase 2** | SISO input-output (full G + Φ_v estimation) | 2–3 days |
| **Phase 3** | Uncertainty estimation and confidence bands | 1–2 days |
| **Phase 4** | MIMO generalization | 2–3 days |
| **Phase 5** | Bode/spectrum plotting with confidence intervals | 1 day |
| **Phase 6** | Frequency-domain input data support | 1–2 days |
| **Phase 7** | Documentation, packaging, MATLAB cross-validation | 2–3 days |

---

## 6. Existing Open-Source Landscape

Before building from scratch, it's worth noting what already exists:

- **`scipy.signal.welch`** — Welch's method for PSD estimation. No transfer function, no uncertainty.
- **`scipy.signal.csd`** — Cross-spectral density. Could be used as a building block, but uses Welch's method (segment-averaging), not Blackman-Tukey (windowed covariance).
- **`spectrum` (PyPI)** — Correlogram PSD estimation exists, but no transfer function estimation or uncertainty.
- **`python-control`** — Control systems library. Has frequency response tools but no non-parametric estimation from data.
- **`sippy` (Systems Identification Package for Python)** — Has some system identification tools but limited spectral analysis.

**Gap:** No existing Python package provides the complete `spa` workflow: Blackman-Tukey-based transfer function estimation with uncertainty from input-output data. This is a genuine gap worth filling.

---

## 7. Design Decisions and Recommendations

**Window choice:** Default to Hann (matching current MATLAB), but expose a `window` parameter accepting `'hann'`, `'hamming'`, or a callable.

**Normalization:** Match MATLAB's System Identification Toolbox normalization, which differs from scipy's conventions. The spectrum is defined as `Φ(ω) = Σ R(τ) e^{-jωτ}` (no `2π` factor, no `Ts` factor in the summation itself). Document this clearly.

**FFT vs direct computation:** Always prefer FFT when using the default frequency grid. For custom frequencies, consider using the Goertzel algorithm if only a few frequencies are needed, or the chirp-z transform for arbitrary frequency grids.

**Numerical stability:** When computing Ĝ = Φ̂_yu / Φ̂_u, protect against division by zero at frequencies where input has no power. Use a regularization floor or return NaN with a warning.

**License:** MIT or BSD-3-Clause recommended for maximum compatibility with both academic and commercial use.

---

## 8. Summary

The `spa` algorithm is conceptually straightforward — it's the Blackman-Tukey correlogram method applied to system identification rather than just spectrum estimation. The core mathematical operations are covariance estimation, windowing, and Fourier transformation. The real value added by a proper implementation lies in:

1. Correctly computing the **transfer function estimate** (not just the power spectrum).
2. Providing **rigorous uncertainty estimates** based on asymptotic theory.
3. Handling the **MIMO case** cleanly.
4. Offering sensible **defaults** (window size, frequency grid) that work well in practice.

All of the underlying mathematics are well-documented in Ljung (1999) and no proprietary algorithms are involved — the method predates MATLAB itself. An open-source implementation is both feasible and fills a real gap in the Python ecosystem.
