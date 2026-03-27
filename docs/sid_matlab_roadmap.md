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
| **`sidFreqBTMap`** | — | Time-varying frequency response map (LTV analysis) |
| **`sidSpectrogram`** | `spectrogram` | Short-time FFT spectrogram |
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
| **`sidMapPlot`** | Time-frequency color map (for sidFreqBTMap results) |
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
├── sidFreqBTMap.m           % Time-varying frequency response map
├── sidSpectrogram.m         % Short-time FFT spectrogram
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

### Phase 7 — `sidFreqBTMap` + `sidSpectrogram` — Time-Varying Analysis (~5 days)

- `sidSpectrogram.m`:
  - Short-time FFT spectrogram (replaces Signal Processing Toolbox `spectrogram`)
  - Windowed segments → FFT → one-sided PSD
  - Supports Hann, Hamming, rectangular, or custom window vector
  - Returns struct with Time, Frequency, Power, PowerDB, Complex coefficients
- `sidFreqBTMap.m`:
  - Segment the data into overlapping windows
  - Run `sidFreqBT` on each segment
  - Collect G(w,t), Phi_v(w,t), coherence(w,t) into 2D arrays
  - Compute time vector from segment centers
  - Share segmentation conventions with `sidSpectrogram` for aligned time axes
- `sidMapPlot.m`:
  - Color map visualization (pcolor/imagesc)
  - Plot types: magnitude, phase, noise, coherence, spectrum
  - Log frequency axis, time on x-axis, colorbar
  - Octave-compatible
- `sidSpectrogramPlot.m`:
  - Standard spectrogram color map (time × frequency × power dB)
  - Shared visual style with `sidMapPlot`
- Tests:
  - `sidSpectrogram`: chirp signal (verify moving peak), white noise (flat), known sinusoid
  - `sidFreqBTMap`: LTI system (constant map), step change in system, chirp
  - Alignment test: verify time axes match between `sidSpectrogram` and `sidFreqBTMap`
  - Compare `sidSpectrogram` output to MathWorks `spectrogram` (if available)

### Phase 8 — `sidFreqETFE` and `sidFreqBTFDR` (~4 days)

- `sidFreqETFE.m` — FFT ratio with optional smoothing
- `sidFreqBTFDR.m` — frequency-dependent window size
- Tests for both

### Phase 9 — Validation + Release (~4 days)

- `exampleCompare.m` — head-to-head vs. MathWorks `spa`
- Octave CI on GitHub Actions
- Edge case hardening (short data, complex, single-precision)
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
| 7. sidFreqBTMap + sidSpectrogram | 5 days | 22 days |
| 8. ETFE + BTFDR | 4 days | 26 days |
| 9. Validation + release | 4 days | 30 days |

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
- Custom window functions (Hann only)
- idfrd-compatible class
- Python / Julia ports
- C reference implementation