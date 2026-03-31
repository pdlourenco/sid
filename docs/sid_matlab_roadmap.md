# sid — Open-Source System Identification Toolbox for MATLAB/Octave

## Naming Convention

All public functions follow the pattern `sid` + `Domain` + `Method`:

```
sid  [Domain]  [Method/Variant]
 │      │          │
 │      │          └── BT, BTFDR, ETFE, ARX, N4SID, AR, ...
 │      └──────────── Freq, TF, SS, TS, LTV, ...
 └─────────────────── system identification (root)
```

### Function Catalog (v1.0 scope in bold)

| Function | Replaces | Description | Status |
|----------|----------|-------------|--------|
| **`sidFreqBT`** | `spa` | Frequency response via Blackman-Tukey | ✅ |
| **`sidFreqBTFDR`** | `spafdr` | Blackman-Tukey, frequency-dependent resolution | ✅ |
| **`sidFreqETFE`** | `etfe` | Empirical transfer function estimate | ✅ |
| **`sidFreqMap`** | `tfestimate`, `mscohere`, `cpsd` | Time-varying frequency response map (BT or Welch) | ✅ |
| **`sidSpectrogram`** | `spectrogram` | Short-time FFT spectrogram | ✅ |
| **`sidLTVdisc`** | — | Discrete LTV state-space identification (COSMIC) | ✅ |
| **`sidLTVdiscTune`** | — | Lambda tuning (validation-based and frequency-response) | ✅ |
| **`sidLTVdiscFrozen`** | — | Frozen transfer function G(ω,k) from A(k), B(k) | ✅ |
| **`sidLTVdiscInit`** | — | Initialize recursive/online COSMIC estimator | ⬜ |
| **`sidLTVdiscUpdate`** | — | Process one time step (filtered estimate) | ⬜ |
| **`sidLTVdiscSmooth`** | — | Backward pass over window (smoothed estimates) | ⬜ |
| **`sidLTVdiscIO`** | — | Partial-observation LTV identification (alternating COSMIC + RTS smoother) | ⬜ |
| **`sidModelOrder`** | — | Model order estimation from frequency response (Hankel SVD) | ⬜ |
| **`sidDetrend`** | `detrend` | Polynomial detrending (preprocessing) | ✅ |
| **`sidResidual`** | `resid` | Residual analysis (whiteness + independence tests) | ✅ |
| **`sidCompare`** | `compare` | Model output comparison with fit metric | ✅ |
| `sidTfARX` | `arx` | Transfer function, ARX model | — |
| `sidTfARMAX` | `armax` | Transfer function, ARMAX model | — |
| `sidSsN4SID` | `n4sid` | State-space, N4SID subspace method | — |
| `sidTsAR` | `ar` | Time series, autoregressive | — |
| `sidTsARMA` | `arma` | Time series, ARMA | — |

### Plotting Functions

| Function | Description | Status |
|----------|-------------|--------|
| **`sidBodePlot`** | Bode diagram with confidence bands | ✅ |
| **`sidSpectrumPlot`** | Power spectrum with confidence bands | ✅ |
| **`sidMapPlot`** | Time-frequency color map (for sidFreqMap results) | ✅ |
| **`sidSpectrogramPlot`** | Spectrogram color map (for sidSpectrogram results) | ✅ |
| `sidNyquistPlot` | Nyquist plot | — |
| `sidPolePlot` | Pole-zero map | — |

### Utility Functions

| Function | Description |
|----------|-------------|
| **`sidCov`** | Biased cross-covariance estimation (single and multi-trajectory) |
| **`sidHannWin`** | Hann lag window |
| **`sidWindowedDFT`** | Windowed Fourier transform (FFT + direct paths) |
| **`sidUncertainty`** | Asymptotic variance formulas |

---

## Package Structure

```
sid-matlab/
├── sidFreqBT.m              % Blackman-Tukey spectral analysis
├── sidFreqBTFDR.m           % BT with frequency-dependent resolution
├── sidFreqETFE.m            % Empirical transfer function estimate
├── sidFreqMap.m              % Time-varying frequency response map (BT or Welch)
├── sidSpectrogram.m         % Short-time FFT spectrogram
├── sidLTVdisc.m             % Discrete LTV state-space identification (COSMIC)
├── sidLTVdiscTune.m         % Lambda tuning via validation or frequency response
├── sidLTVdiscFrozen.m       % Frozen transfer function from A(k), B(k)
├── sidLTVdiscInit.m         % Initialize recursive COSMIC estimator              (planned)
├── sidLTVdiscUpdate.m       % Online: process one time step                      (planned)
├── sidLTVdiscSmooth.m       % Windowed backward pass for smoothed estimates      (planned)
├── sidLTVdiscIO.m           % Partial-observation LTV identification              (planned)
├── sidModelOrder.m          % Model order estimation (Hankel SVD)                 (planned)
├── sidBodePlot.m            % Bode diagram with confidence bands
├── sidSpectrumPlot.m        % Power spectrum plot
├── sidMapPlot.m             % Time-frequency color map
├── sidSpectrogramPlot.m     % Spectrogram color map
├── sidDetrend.m             % Polynomial detrending (preprocessing)
├── sidResidual.m            % Residual analysis (whiteness + independence)
├── sidCompare.m             % Model output comparison with fit metric
├── internal/
│   ├── sidCov.m             % Biased covariance estimation
│   ├── sidDFT.m             % DFT computation (FFT + direct paths)
│   ├── sidHannWin.m         % Hann window generation
│   ├── sidWindowedDFT.m     % Windowed DFT (FFT + direct)
│   ├── sidUncertainty.m     % Asymptotic variance formulas
│   ├── sidValidate.m        % Input parsing and validation
│   └── sidValidateData.m    % Data validation helper
├── tests/
│   ├── runAllTests.m        % Master test runner
│   ├── test_sidFreqBT.m     % SISO + time series + MIMO
│   ├── test_sidFreqBTFDR.m
│   ├── test_sidFreqETFE.m
│   ├── test_sidFreqMap.m
│   ├── test_sidSpectrogram.m
│   ├── test_sidMapPlot.m
│   ├── test_sidSpectrogramPlot.m
│   ├── test_sidHannWin.m
│   ├── test_sidCov.m
│   ├── test_sidDFT.m
│   ├── test_sidWindowedDFT.m
│   ├── test_sidUncertainty.m
│   ├── test_sidValidate.m
│   ├── test_sidPlotting.m
│   ├── test_validation.m
│   ├── test_crossMethod.m
│   ├── test_compareSpa.m    % vs. MathWorks spa
│   ├── test_compareEtfe.m   % vs. MathWorks etfe
│   ├── test_compareSpafdr.m % vs. MathWorks spafdr
│   ├── test_compareWelch.m  % vs. MathWorks tfestimate/mscohere
│   ├── test_compareMultiTraj.m
│   ├── test_multiTrajectory.m
│   ├── test_sidLTVdisc.m
│   ├── test_sidLTVdiscTune.m
│   ├── test_sidLTVdiscUncertainty.m
│   ├── test_sidLTVdiscFrozen.m
│   ├── test_sidLTVdiscVarLen.m
│   ├── test_sidDetrend.m
│   ├── test_sidResidual.m
│   ├── test_sidCompare.m
│   └── test_sidModelOrder.m                                     (planned)
├── examples/
│   ├── README.md
│   ├── exampleSISO.m
│   ├── exampleMIMO.m
│   ├── exampleETFE.m
│   ├── exampleFreqDepRes.m
│   ├── exampleFreqMap.m
│   ├── exampleSpectrogram.m
│   ├── exampleCoherence.m
│   ├── exampleMethodComparison.m
│   ├── exampleLTVdisc.m
│   ├── exampleMultiTrajectory.m
│   ├── exampleOutputCOSMIC.m                                    (planned)
│   └── runAllExamples.m
├── docs/
│   ├── sid_matlab_roadmap.md
│   ├── cosmic_uncertainty_derivation.md
│   ├── cosmic_online_recursion.md
│   ├── cosmic_automatic_tuning.md
│   ├── multi_trajectory_spectral_theory.md
│   └── cosmic_output.md
├── SPEC.md
├── LICENSE                    % MIT
├── README.md
└── sidInstall.m               % Adds sid to MATLAB/Octave path
```

### Usage

```matlab
% Setup (once)
sidInstall

% Basic SISO frequency response estimation
result = sidFreqBT(y, u);
sidBodePlot(result);

% With options
result = sidFreqBT(y, u, 'WindowSize', 50, ...
                         'Frequencies', logspace(-2, pi, 256));

% Time series (output spectrum only)
result = sidFreqBT(y, []);
sidSpectrumPlot(result, 'Confidence', 3);

% Empirical transfer function estimate
result = sidFreqETFE(y, u, 'Smoothing', 5);
```

---

## Result Struct

All `sidFreq*` functions return the same struct:

```matlab
result.Frequency          % (n_freq x 1) rad/sample
result.FrequencyHz        % (n_freq x 1) Hz
result.Response           % (n_freq x n_out x n_in) complex
result.ResponseStd        % (n_freq x n_out x n_in) real
result.NoiseSpectrum      % (n_freq x n_out x n_out) real
result.NoiseSpectrumStd   % (n_freq x n_out x n_out) real
result.Coherence          % (n_freq x 1) squared coherence (SISO only, [] for MIMO)
result.SampleTime         % scalar (seconds)
result.WindowSize         % scalar integer (or vector for BTFDR)
result.DataLength         % N (number of samples used)
result.Method             % 'sidFreqBT', 'sidFreqBTFDR', 'sidFreqETFE', or 'sidFreqMap'
```

---

## Revised Roadmap

### Phase 1 — Spec + Scaffolding (~3 days) ✅

- Write `SPEC.md`: exact formulas, defaults, edge cases, normalization
- Create package folder structure
- Write `sidInstall.m`
- Stub every v1.0 function with full help text
- Set up test framework (MATLAB `runtests` compatible)

### Phase 2 — `sidFreqBT` SISO Core (~4 days) ✅

- `sidCov.m` — biased cross-covariance for lags 0..M
- `sidHannWin.m` — Hann lag window
- `sidWindowedDFT.m` — FFT fast path + direct DFT slow path
- `sidValidate.m` — input parsing (supports both positional and name-value)
- `sidFreqBT.m` — full SISO pipeline:
  covariance → window → DFT → G = Phi_yu/Phi_u → Phi_v
- Tests against known systems and analytic solutions

### Phase 3 — Time Series Mode (~1 day) ✅

- `sidFreqBT.m` handles empty `u`: returns Phi_y only
- Tests: AR(1) known spectrum, white noise, sinusoid in noise

### Phase 4 — Uncertainty (~3 days) ✅

- `sidUncertainty.m` — asymptotic variance for G and Phi_v
- Coherence computation
- Monte Carlo validation (empirical std vs. formula)
- Result struct extended with `ResponseStd`, `NoiseSpectrumStd`

### Phase 5 — Plotting (~2 days) ✅

- `sidBodePlot.m` — magnitude + phase, log-freq, shaded confidence
- `sidSpectrumPlot.m` — power spectrum (dB) with confidence
- Octave-compatible (subplot, not tiledlayout)
- Return figure/axes handles

### Phase 6 — MIMO (~4 days) ✅

- `sidCov.m` extended to matrix-valued covariances
- `sidFreqBT.m` extended: full spectral matrix, matrix inversion
- `sidBodePlot.m` extended: subplot grid per channel pair
- Tests: 2x2 known system, verify channel-by-channel

### Phase 7 — `sidFreqMap` + `sidSpectrogram` — Time-Varying Analysis (~6 days) ✅

- `sidFreqMap.m`:
  - Outer segmentation: sliding overlapping windows (shared by both algorithms)
  - `'Algorithm', 'bt'` (default): calls `sidFreqBT` per segment (BT correlogram) ✅
  - `'Algorithm', 'welch'`: Welch's method per segment (replaces `tfestimate`/`mscohere`/`cpsd`) ⬜
    - Sub-segmentation within each outer segment
    - Time-domain window (Hann default), FFT, averaged cross/auto periodograms
    - Form G, Phi_v, coherence from averaged spectra
  - Identical output struct regardless of algorithm
  - Share segmentation conventions with `sidSpectrogram` for aligned time axes
- `sidSpectrogram.m` ✅
- `sidMapPlot.m` ✅
- `sidSpectrogramPlot.m` ✅
- Tests:
  - `sidFreqMap` BT: LTI system (constant map), step change, chirp ✅
  - `sidFreqMap` Welch: same tests, verify qualitatively similar results to BT ⬜
  - `sidFreqMap` Welch vs. MathWorks `tfestimate`: numerical comparison ⬜
  - `sidSpectrogram`: chirp signal (verify moving peak), white noise (flat), known sinusoid ✅
  - Alignment test: verify time axes match between `sidSpectrogram` and `sidFreqMap` ✅
  - Compare `sidSpectrogram` output to MathWorks `spectrogram` (if available) ✅

### Phase 8 — `sidLTVdisc` Base (~5 days) ✅

- Integrate existing COSMIC MATLAB implementation into `sid` conventions
- `sidLTVdisc.m`:
  - COSMIC forward-backward pass (closed-form block tridiagonal solver)
  - Multi-trajectory support (same horizon)
  - Preconditioning option
  - L-curve automatic lambda selection
  - Manual scalar or per-step lambda
  - Returns struct with A(k), B(k), cost breakdown
- `sidLTVdiscTune.m`:
  - Grid search over lambda candidates
  - Evaluate trajectory prediction loss on validation data
  - Return best result, best lambda, all losses
- Tests:
  - Known LTI system: verify A(k), B(k) are constant
  - Known LTV system: compare to ground truth
  - Multi-trajectory vs. single-trajectory accuracy
  - L-curve lambda selection: verify reasonable choice
  - Preconditioning and uniqueness condition checks

### Phase 8a — Variable-Length Trajectories (~2 days) ✅

- Extend input parsing to accept cell arrays of different-length trajectories
- Modify `buildDataMatrices` to handle per-step active trajectory sets
- Normalization: `1/sqrt(|L(k)|)` per step
- Tests: mix of short and long trajectories, verify identical to uniform when all same length

### Phase 8b — Bayesian Uncertainty (~4 days) ✅

**Theory:** `docs/cosmic_uncertainty_derivation.md`

- Implement uncertainty backward pass reusing stored Λ_k matrices:
  - `P(N-1) = Λ_{N-1}⁻¹`
  - `P(k) = Λ_k⁻¹ + G_k P(k+1) G_k^T` where `G_k = λ_{k+1} Λ_k⁻¹`
  - Same O(N(p+q)³) cost as COSMIC itself
- Noise variance estimation: `σ̂² = 2h(C*) / (NLp)`
- Add `AStd`, `BStd`, `Covariance`, `NoiseVariance` to result struct
- `sidLTVdiscFrozen.m`: compute `G(ω, k) = (zI - A(k))⁻¹ B(k)` with Jacobian-propagated uncertainty
- Monte Carlo validation: 500 realizations, verify empirical std matches posterior
- Integration with `sidBodePlot` for time-varying Bode with confidence bands

### Phase 8c — Online/Recursive COSMIC (~4 days) ⬜

**Theory:** `docs/cosmic_online_recursion.md`

Key insight: COSMIC forward pass = Kalman filter on parameter evolution; backward pass = RTS smoother. Three operating modes:

- `sidLTVdiscInit.m`: initialize recursive estimator (Λ₀, Y₀)
- `sidLTVdiscUpdate.m`: process one time step in O((p+q)³):
  - One forward pass step: `Λ_k`, `Y_k` from new data + previous state
  - Returns filtered estimate `Y_k` and filtered uncertainty `Λ_k⁻¹`
- `sidLTVdiscSmooth.m`: backward pass over stored window for smoothed estimates
- Tests:
  - Filtered estimates converge to smoothed as window grows
  - Filtered uncertainty ≥ smoothed uncertainty at every step
  - Process data one-at-a-time, compare final result to batch COSMIC
  - Timing: O(1) per step regardless of history length

### Phase 8d — Lambda Tuning via Frequency Response (~4 days) ✅

- Extend `sidLTVdiscTune` with `'Method', 'frequency'` option
- Frozen transfer function vs. `sidFreqMap` comparison
- Mahalanobis consistency scoring at each (ω, t) grid point
- Select largest λ where ≥90% of grid points consistent at 95% level
- Depends on: Phase 7 (`sidFreqMap`) and Phase 8b (uncertainty)
- Tests: known LTV system, verify selected lambda is reasonable

### Phase 8e — Output-COSMIC (`sidLTVdiscIO`) (~6 days) ⬜

**Theory:** `docs/cosmic_output.md`

Core: alternating minimisation of joint objective (observation fidelity + dynamics fidelity + COSMIC smoothness). When `H = I`, reduces to standard `sidLTVdisc`.

- `sidModelOrder.m`:
  - Input: `sidFreqBT` (or any `sidFreq*`) result struct
  - Compute impulse response coefficients via IFFT of `Response`
  - Build block Hankel matrix from impulse response
  - SVD → singular values
  - Detect model order from largest singular value ratio gap (default) or user-specified threshold
  - Return `n`, singular values, and Hankel matrix (for inspection)
  - Optional: `'Plot', true` for singular value bar chart with detected gap
- `sidLTVdiscIO.m`:
  - Initialisation: single forward-backward pass for states `{x_l(k)}` and input matrices `{B(k)}` jointly with `A = I` (exact minimisation of `J|_{A=I}`, jointly convex, composite block tridiagonal)
  - COSMIC step: standard `sidLTVdisc` on estimated states (reuses existing implementation)
  - State step: RTS smoother per trajectory with current `A(k)`, `B(k)`, `H`, `R`
  - Alternating loop with convergence monitoring (`|ΔJ/J| < ε_J`)
  - Optional trust-region interpolation `Ã(k) = (1-μ) A(k) + μ I` with adaptive μ schedule
  - Multi-trajectory support: shared `C(k)`, independent state sequences
  - Variable-length trajectories via cell arrays
  - Returns estimated `A(k)`, `B(k)`, `X(k)`, cost history, iteration count
  - Uncertainty: run COSMIC uncertainty backward pass (Phase 8b) at final iteration
- `exampleOutputCOSMIC.m`:
  - End-to-end workflow: `sidFreqBT` → `sidModelOrder` → construct `H` → `sidLTVdiscIO`
  - Demonstrates partially-known H case (some states measured, hidden states estimated)
- Tests:
  - `sidModelOrder`: known 4th-order system, verify n = 4 detected; known 2nd-order, verify n = 2
  - Known LTV system with `H = I`: verify identical to `sidLTVdisc`
  - Known LTV system with `H ≠ I`: compare estimated `A(k)`, `B(k)` to ground truth
  - Convergence: verify monotone cost decrease across iterations
  - Multi-trajectory: verify pooling improves estimates vs single trajectory
  - Noisy measurements: verify `R` weighting improves estimates vs `R = I`
  - Trust-region: verify convergence on difficult case where `μ = 0` diverges
  - State recovery: compare estimated `x̂(k)` to true states

### Phase 9 — `sidFreqETFE` and `sidFreqBTFDR` (~4 days) ✅

- `sidFreqETFE.m` — FFT ratio with optional smoothing
- `sidFreqBTFDR.m` — frequency-dependent window size
- Tests for both

### Phase 9a — Multi-Trajectory Support for Frequency Functions (~3 days) ✅

**Theory:** `docs/multi_trajectory_spectral_theory.md`

- Extend `sidFreqBT`, `sidFreqETFE` input parsing to accept 3D arrays `(N × n_ch × L)` and cell arrays
- Implement ensemble-averaged covariance in `sidCov` (average per-trajectory covariances across 3rd dimension)
- Extend `sidFreqMap` to ensemble-average per-segment spectra across trajectories
- Extend `sidSpectrogram` to ensemble-average per-segment PSDs (ERSP-like)
- Handle variable-length trajectories (per-segment active trajectory count)
- Add `NumTrajectories` field to all output structs
- Update `sidUncertainty` to use `L` in variance formulas (`1/(L×N)` instead of `1/N`)
- Tests:
  - LTI system: verify L-trajectory estimate matches single trajectory with L× data length
  - LTV system: verify ensemble-averaged `sidFreqMap` + `sidLTVdisc` use same data consistently
  - Variable-length: verify graceful handling when not all trajectories span all segments
  - Variance check: Monte Carlo confirm 1/L variance reduction

### Phase 11 — Workflow Utilities (~4 days) ✅

- `sidDetrend.m`:
  - Polynomial detrending: remove mean (d=0), linear trend (d=1, default), or higher order
  - Segment-wise option for long records
  - Returns detrended data and removed trend
  - Multi-channel support (detrend each channel independently)
- `sidResidual.m`:
  - Compute residuals from any sid model + data
  - For frequency-domain models: IFFT filtering of input through Ĝ(ω)
  - For state-space models: forward propagation x̂(k+1) = A(k)x̂(k) + B(k)u(k)
  - Whiteness test: normalised autocorrelation within ±2.58/sqrt(N)
  - Independence test: normalised cross-correlation with input within same bounds
  - Optional diagnostic plot (two-panel: autocorrelation + cross-correlation)
- `sidCompare.m`:
  - Simulate model output and overlay on measured data
  - NRMSE fit percentage: 100 × (1 - ||y - ŷ|| / ||y - mean(y)||)
  - Works with frequency-domain and state-space models
  - Multi-channel and multi-trajectory support
  - Optional comparison plot
- Tests:
  - sidDetrend: known linear trend, verify removal to machine precision
  - sidResidual: good model (white residuals) vs. wrong-order model (coloured residuals)
  - sidCompare: perfect model gives 100% fit, mean predictor gives 0%
  - Cross-method: detrend → estimate → residual → compare end-to-end workflow

### Phase 10 — Validation, Freeze + Release (~4 days) 🔄

- `exampleCompare.m` — head-to-head vs. MathWorks `spa`
- Octave CI on GitHub Actions
- Edge case hardening
- README, examples, MATLAB File Exchange submission

---

## Timeline

| Phase | Effort | Running Total | Status |
|-------|--------|---------------|--------|
| 1. Spec + scaffolding | 3 days | 3 days | ✅ |
| 2. sidFreqBT SISO | 4 days | 7 days | ✅ |
| 3. Time series | 1 day | 8 days | ✅ |
| 4. Uncertainty | 3 days | 11 days | ✅ |
| 5. Plotting | 2 days | 13 days | ✅ |
| 6. MIMO | 4 days | 17 days | ✅ |
| 7. sidFreqMap + sidSpectrogram | 6 days | 23 days | ✅ |
| 8. sidLTVdisc base | 5 days | 28 days | ✅ |
| 8a. Variable-length trajectories | 2 days | 30 days | ✅ |
| 8b. Bayesian uncertainty | 4 days | 34 days | ✅ |
| 8c. Online/recursive COSMIC | 4 days | 38 days | ⬜ |
| 8d. Lambda via frequency response | 4 days | 42 days | ✅ |
| 8e. Output-COSMIC (`sidLTVdiscIO`) | 6 days | 48 days | ⬜ |
| 9. ETFE + BTFDR | 4 days | 52 days | ✅ |
| 9a. Multi-trajectory spectral | 3 days | 55 days | ✅ |
| 11. Workflow utilities | 4 days | 59 days | ✅ |
| 10. Validation, freeze + release | 4 days | 63 days | 🔄 |

---

## Octave Compatibility Rules

- No `"string"` literals — use `'char vectors'`
- No `arguments` blocks — use `inputParser`
- No `tiledlayout`/`nexttile` — use `subplot`
- No `exportgraphics` — use `print`
- No `dictionary` — use `struct` or `containers.Map`
- Test in CI with Octave 8+

---

## Out of Scope for v1.0

- Frequency-domain input data
- Continuous-time models
- `maxSize` data segmentation
- Custom window functions (Hann only for sidFreqBT)
- idfrd-compatible class
- Python / Julia ports (v1.0 freeze → port → release)
- C reference implementation
- Online/recursive COSMIC (Phase 8c — v2)
- Unknown observation matrix estimation (joint H + dynamics)
- Time-varying observation matrix H(k)
- Alternative regularization norms (non-squared L2, L1 total variation)
- Alternative LTV algorithms (TVERA, TVOKID) — `'Algorithm'` parameter is ready
- GCV lambda selection
- Parametric identification: ARX, ARMAX, state-space subspace methods (v2)
- LPV identification: structured parameter-varying models (v2)
