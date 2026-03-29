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

| Function | Replaces | Description |
|----------|----------|-------------|
| **`sidFreqBT`** | `spa` | Frequency response via Blackman-Tukey |
| **`sidFreqBTFDR`** | `spafdr` | Blackman-Tukey, frequency-dependent resolution |
| **`sidFreqETFE`** | `etfe` | Empirical transfer function estimate |
| **`sidFreqMap`** | `tfestimate`, `mscohere`, `cpsd` | Time-varying frequency response map (BT or Welch) |
| **`sidSpectrogram`** | `spectrogram` | Short-time FFT spectrogram |
| **`sidLTVdisc`** | — | Discrete LTV state-space identification (COSMIC) |
| **`sidLTVdiscTune`** | — | Lambda tuning (validation-based and frequency-response) |
| **`sidLTVdiscFrozen`** | — | Frozen transfer function G(ω,k) from A(k), B(k) |
| **`sidLTVdiscInit`** | — | Initialize recursive/online COSMIC estimator |
| **`sidLTVdiscUpdate`** | — | Process one time step (filtered estimate) |
| **`sidLTVdiscSmooth`** | — | Backward pass over window (smoothed estimates) |
| **`sidLTVdiscIO`** | — | Output-only LTV identification (two-stage) |
| `sidTfARX` | `arx` | Transfer function, ARX model |
| `sidTfARMAX` | `armax` | Transfer function, ARMAX model |
| `sidSsN4SID` | `n4sid` | State-space, N4SID subspace method |
| `sidTsAR` | `ar` | Time series, autoregressive |
| `sidTsARMA` | `arma` | Time series, ARMA |
| `sidLtvDisc` | — | Discrete LTV identification |

### Plotting Functions

| Function | Description |
|----------|-------------|
| **`sidBodePlot`** | Bode diagram with confidence bands |
| **`sidSpectrumPlot`** | Power spectrum with confidence bands |
| **`sidMapPlot`** | Time-frequency color map (for sidFreqMap results) |
| **`sidSpectrogramPlot`** | Spectrogram color map (for sidSpectrogram results) |
| `sidNyquistPlot` | Nyquist plot |
| `sidPolePlot` | Pole-zero map |

### Utility Functions

| Function | Description |
|----------|-------------|
| **`sidCov`** | Biased cross-covariance estimation |
| **`sidHannWin`** | Hann lag window |
| **`sidWindowedDFT`** | Windowed Fourier transform (FFT + direct paths) |
| **`sidUncertainty`** | Asymptotic variance formulas |

---

## Package Structure

```
sid/
├── sidFreqBT.m              % Blackman-Tukey spectral analysis
├── sidFreqBTFDR.m           % BT with frequency-dependent resolution
├── sidFreqETFE.m            % Empirical transfer function estimate
├── sidFreqMap.m           % Time-varying frequency response map
├── sidSpectrogram.m         % Short-time FFT spectrogram
├── sidLTVdisc.m             % Discrete LTV state-space identification
├── sidLTVdiscTune.m         % Lambda tuning via validation or frequency response
├── sidLTVdiscFrozen.m       % Frozen transfer function from A(k), B(k)
├── sidLTVdiscInit.m         % Initialize recursive COSMIC estimator
├── sidLTVdiscUpdate.m       % Online: process one time step
├── sidLTVdiscSmooth.m       % Windowed backward pass for smoothed estimates
├── sidLTVdiscIO.m           % Output-only LTV identification (two-stage)
├── sidBodePlot.m            % Bode diagram with confidence bands
├── sidSpectrumPlot.m        % Power spectrum plot
├── sidMapPlot.m             % Time-frequency color map
├── sidSpectrogramPlot.m     % Spectrogram color map
├── internal/
│   ├── sidCov.m             % Biased covariance estimation
│   ├── sidHannWin.m         % Hann window generation
│   ├── sidWindowedDFT.m     % Windowed DFT (FFT + direct)
│   ├── sidUncertainty.m     % Asymptotic variance formulas
│   └── sidValidate.m        % Input parsing and validation
├── test/
│   ├── testSidFreqBT.m      % SISO + time series + MIMO
│   ├── testSidFreqBTFDR.m
│   ├── testSidFreqETFE.m
│   ├── testSidFreqBTMap.m
│   ├── testSidLTVdisc.m
│   ├── testSidUncertainty.m
│   ├── testSidEdgeCases.m
│   ├── testSidOctave.m
│   └── fixtures/
│       └── generateFixtures.m
├── examples/
│   ├── exampleSISO.m
│   ├── exampleTimeSeries.m
│   └── exampleCompare.m      % vs. MathWorks spa if available
├── SPEC.md
├── LICENSE                    % MIT
├── README.md
└── sidInstall.m               % Adds +sid to path
```

### Usage

```matlab
% Setup (once)
sidInstall

% Basic SISO frequency response estimation
result = sid.sidFreqBT(y, u);
sid.sidBodePlot(result);

% With options
result = sid.sidFreqBT(y, u, 'WindowSize', 50, ...
                             'Frequencies', logspace(-2, pi, 256));

% Time series (output spectrum only)
result = sid.sidFreqBT(y, []);
sid.sidSpectrumPlot(result, 'Confidence', 3);

% Empirical transfer function estimate
result = sid.sidFreqETFE(y, u, 'Smoothing', 5);
```

---

## Result Struct

All `sidFreq*` functions return the same struct:

```matlab
result.Frequency          % (n_freq x 1) rad/sample
result.Response           % (n_freq x n_out x n_in) complex
result.ResponseStd        % (n_freq x n_out x n_in) real
result.NoiseSpectrum      % (n_freq x n_out x n_out) real
result.NoiseSpectrumStd   % (n_freq x n_out x n_out) real
result.SampleTime         % scalar (seconds)
result.WindowSize         % scalar integer
result.Method             % 'sidFreqBT', 'sidFreqBTFDR', or 'sidFreqETFE'
result.DataLength          % N (number of samples used)
```

---

## Revised Roadmap

### Phase 1 — Spec + Scaffolding (~3 days)

- Write `SPEC.md`: exact formulas, defaults, edge cases, normalization
- Create package folder structure
- Write `sidInstall.m`
- Stub every v1.0 function with full help text
- Set up test framework (MATLAB `runtests` compatible)

### Phase 2 — `sidFreqBT` SISO Core (~4 days)

- `sidCov.m` — biased cross-covariance for lags 0..M
- `sidHannWin.m` — Hann lag window
- `sidWindowedDFT.m` — FFT fast path + direct DFT slow path
- `sidValidate.m` — input parsing (supports both positional and name-value)
- `sidFreqBT.m` — full SISO pipeline:
  covariance → window → DFT → G = Phi_yu/Phi_u → Phi_v
- Tests against known systems and analytic solutions

### Phase 3 — Time Series Mode (~1 day)

- `sidFreqBT.m` handles empty `u`: returns Phi_y only
- Tests: AR(1) known spectrum, white noise, sinusoid in noise

### Phase 4 — Uncertainty (~3 days)

- `sidUncertainty.m` — asymptotic variance for G and Phi_v
- Coherence computation
- Monte Carlo validation (empirical std vs. formula)
- Result struct extended with `ResponseStd`, `NoiseSpectrumStd`

### Phase 5 — Plotting (~2 days)

- `sidBodePlot.m` — magnitude + phase, log-freq, shaded confidence
- `sidSpectrumPlot.m` — power spectrum (dB) with confidence
- Octave-compatible (subplot, not tiledlayout)
- Return figure/axes handles

### Phase 6 — MIMO (~4 days)

- `sidCov.m` extended to matrix-valued covariances
- `sidFreqBT.m` extended: full spectral matrix, matrix inversion
- `sidBodePlot.m` extended: subplot grid per channel pair
- Tests: 2x2 known system, verify channel-by-channel

### Phase 7 — `sidFreqMap` + `sidSpectrogram` — Time-Varying Analysis (~6 days)

- `sidFreqMap.m`:
  - Outer segmentation: sliding overlapping windows (shared by both algorithms)
  - `'Algorithm', 'bt'` (default): calls `sidFreqBT` per segment (BT correlogram)
  - `'Algorithm', 'welch'`: Welch's method per segment (replaces `tfestimate`/`mscohere`/`cpsd`)
    - Sub-segmentation within each outer segment
    - Time-domain window (Hann default), FFT, averaged cross/auto periodograms
    - Form G, Phi_v, coherence from averaged spectra
  - Identical output struct regardless of algorithm
  - Share segmentation conventions with `sidSpectrogram` for aligned time axes
- `sidSpectrogram.m`:
  - Short-time FFT spectrogram (replaces Signal Processing Toolbox `spectrogram`)
  - Windowed segments → FFT → one-sided PSD
  - Supports Hann, Hamming, rectangular, or custom window vector
  - Returns struct with Time, Frequency, Power, PowerDB, Complex coefficients
- `sidMapPlot.m`:
  - Color map visualization (pcolor/imagesc)
  - Plot types: magnitude, phase, noise, coherence, spectrum
  - Log frequency axis, time on x-axis, colorbar
  - Works identically with BT and Welch results
  - Octave-compatible
- `sidSpectrogramPlot.m`:
  - Standard spectrogram color map (time × frequency × power dB)
  - Shared visual style with `sidMapPlot`
- Tests:
  - `sidFreqMap` BT: LTI system (constant map), step change, chirp
  - `sidFreqMap` Welch: same tests, verify qualitatively similar results to BT
  - `sidFreqMap` Welch vs. MathWorks `tfestimate`: numerical comparison (if available)
  - `sidSpectrogram`: chirp (moving peak), white noise (flat), known sinusoid
  - Alignment test: verify time axes match between `sidSpectrogram` and `sidFreqMap`
  - Compare `sidSpectrogram` output to MathWorks `spectrogram` (if available)

### Phase 8 — `sidLTVdisc` Base (~5 days)

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

### Phase 8a — Variable-Length Trajectories (~2 days)

- Extend input parsing to accept cell arrays of different-length trajectories
- Modify `buildDataMatrices` to handle per-step active trajectory sets
- Normalization: `1/sqrt(|L(k)|)` per step
- Tests: mix of short and long trajectories, verify identical to uniform when all same length

### Phase 8b — Bayesian Uncertainty (~4 days)

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

### Phase 8c — Online/Recursive COSMIC (~4 days)

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

### Phase 8d — Lambda Tuning via Frequency Response (~4 days)

- Extend `sidLTVdiscTune` with `'Method', 'frequency'` option
- Frozen transfer function vs. `sidFreqMap` comparison
- Mahalanobis consistency scoring at each (ω, t) grid point
- Select largest λ where ≥90% of grid points consistent at 95% level
- Depends on: Phase 7 (`sidFreqMap`) and Phase 8b (uncertainty)
- Tests: known LTV system, verify selected lambda is reasonable

### Phase 8e — Output-Only Estimation (~3 days)

- `sidLTVdiscIO.m`: two-stage (initial LTI observer → state reconstruction → COSMIC)
- Warning to user about approximate state estimates
- Tests: known system with measured outputs only, compare to full-state COSMIC

### Phase 9 — `sidFreqETFE` and `sidFreqBTFDR` (~4 days)

- `sidFreqETFE.m` — FFT ratio with optional smoothing
- `sidFreqBTFDR.m` — frequency-dependent window size
- Tests for both

### Phase 10 — Validation + Release (~4 days)

- `exampleCompare.m` — head-to-head vs. MathWorks `spa`
- Octave CI on GitHub Actions
- Edge case hardening
- README, examples, MATLAB File Exchange submission

---

## Timeline

| Phase | Effort | Running Total |
|-------|--------|---------------|
| 1. Spec + scaffolding | 3 days | 3 days |
| 2. sidFreqBT SISO | 4 days | 7 days |
| 3. Time series | 1 day | 8 days |
| 4. Uncertainty | 3 days | 11 days |
| 5. Plotting | 2 days | 13 days |
| 6. MIMO | 4 days | 17 days |
| 7. sidFreqMap + sidSpectrogram | 5 days | 22 days |
| 8. sidLTVdisc base | 5 days | 27 days |
| 8a. Variable-length trajectories | 2 days | 29 days |
| 8b. Bayesian uncertainty | 4 days | 33 days |
| 8c. Online/recursive COSMIC | 4 days | 37 days |
| 8d. Lambda via frequency response | 4 days | 41 days |
| 8e. Output-only (two-stage) | 3 days | 44 days |
| 9. ETFE + BTFDR | 4 days | 48 days |
| 10. Validation + release | 4 days | 52 days |

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
- Python / Julia ports
- C reference implementation
- EM-style or direct output equation LTV identification
- Alternative regularization norms (non-squared L2, L1 total variation)
- Alternative LTV algorithms (TVERA, TVOKID) — `'Algorithm'` parameter is ready
- GCV lambda selection
