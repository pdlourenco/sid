# sid v0.1.0 — matlab/octave — Release Notes

**Release date:** 2026-04-11
**License:** MIT

## Overview

sid is a free, open-source system identification toolbox for MATLAB and
GNU Octave. It provides non-parametric frequency response estimation and
time-varying state-space identification with no dependency on MathWorks'
System Identification or Signal Processing Toolboxes.

The toolbox covers two complementary paths: a frequency-domain path built on
Blackman-Tukey spectral analysis (Ljung, 1999), and a state-space path built
on the COSMIC (Closed-form data-driven linear time-varying SysteM IdentifiCation)
algorithm — a closed-form, O(N)-complexity solver for discrete linear time-varying
system identification (Carvalho et al., 2022). Both paths support SISO, MIMO, time
series, and multi-trajectory data through a unified API.

No toolboxes are required. The entire codebase uses only core MATLAB/Octave
functions and runs on MATLAB R2016b+ and GNU Octave 8.0+, including MATLAB
Online.

## Release Contents

```
sid-v0.1.0-matlab/
├── README.md                           # Overview, features, references
├── LICENSE                             # MIT License
├── spec/
│   ├── SPEC.md                         # Full algorithm specification
│   └── EXAMPLES.md                     # Example-suite specification
└── matlab/
    ├── README.md                       # Installation, quick start, function reference
    ├── sidInstall.m                    # Setup script (adds sid/ to path)
    ├── sid/                            # Public API
    │   ├── sidFreqBT.m                 # Blackman-Tukey spectral estimation
    │   ├── sidFreqBTFDR.m              # BT with frequency-dependent resolution
    │   ├── sidFreqETFE.m               # Empirical transfer function estimate
    │   ├── sidFreqMap.m                # Time-varying frequency response map
    │   ├── sidSpectrogram.m            # Short-time FFT spectrogram
    │   ├── sidLTVdisc.m                # COSMIC: LTV state-space identification
    │   ├── sidLTVdiscTune.m            # Regularisation tuning (validation, frequency)
    │   ├── sidLTVdiscFrozen.m          # Frozen transfer function G(ω,k)
    │   ├── sidLTVdiscIO.m              # Output-COSMIC: partial-observation LTV ID
    │   ├── sidLTVStateEst.m            # Batch LTV state estimation (RTS smoother)
    │   ├── sidLTIfreqIO.m              # LTI realisation from I/O freq. response
    │   ├── sidModelOrder.m             # Model order estimation (Hankel SVD)
    │   ├── sidDetrend.m                # Polynomial detrending
    │   ├── sidResidual.m               # Residual analysis (whiteness + independence)
    │   ├── sidCompare.m                # Model output comparison (NRMSE fit)
    │   ├── sidBodePlot.m               # Bode diagram with confidence bands
    │   ├── sidSpectrumPlot.m           # Power spectrum with confidence bands
    │   ├── sidMapPlot.m                # Time-frequency colour map
    │   ├── sidSpectrogramPlot.m        # Spectrogram colour map
    │   ├── sidResultTypes.m            # Result struct documentation
    │   └── private/                    # Internal helper functions
    │       ├── sidCov.m                #   Biased cross-covariance estimation
    │       ├── sidDFT.m                #   Discrete Fourier transform
    │       ├── sidHannWin.m            #   Hann lag window
    │       ├── sidWindowedDFT.m        #   Windowed FFT (fast + direct paths)
    │       ├── sidUncertainty.m        #   Asymptotic variance formulas
    │       ├── sidValidateData.m       #   Input validation and parsing
    │       ├── sidParseOptions.m       #   Name-value option parsing
    │       ├── sidIsDefaultFreqs.m     #   Default frequency grid detection
    │       ├── sidFreqDomainSim.m      #   Frequency-domain model simulation
    │       ├── sidEstimateNoiseCov.m   #   Noise covariance estimation
    │       ├── sidExtractStd.m         #   Standard deviation extraction
    │       ├── sidTestMSD.m            #   Mass-spring-damper test system
    │       ├── sidLTVbuildDataMatrices.m       # COSMIC data matrix construction
    │       ├── sidLTVbuildDataMatricesVarLen.m # Variable-length trajectory support
    │       ├── sidLTVbuildBlockTerms.m         # Block tridiagonal term assembly
    │       ├── sidLTVcosmicSolve.m             # COSMIC forward-backward solve
    │       ├── sidLTVblkTriSolve.m             # Generic block tridiagonal solver
    │       ├── sidLTVevaluateCost.m            # COSMIC cost function evaluation
    │       └── sidLTVuncertaintyBackwardPass.m # Posterior covariance recursion
    └── examples/
        ├── README.md                   # Examples guide
        ├── runAllExamples.m            # Run all examples
        ├── util_msd.m                  # ZOH discretisation of n-mass SMD chain
        ├── util_msd_ltv.m              # Per-step (Ad(k), Bd(k)) for LTV SMD chain
        ├── util_msd_nl.m               # RK4 simulator for Duffing-style nonlinear SMD
        ├── exampleSISO.m               # Basic SISO frequency response
        ├── exampleETFE.m               # Empirical transfer function
        ├── exampleFreqDepRes.m         # Frequency-dependent resolution
        ├── exampleCoherence.m          # Coherence analysis
        ├── exampleMethodComparison.m   # BT vs. BTFDR vs. ETFE
        ├── exampleMIMO.m                # MIMO estimation
        ├── exampleFreqMap.m             # Time-varying frequency maps
        ├── exampleSpectrogram.m         # Spectrogram
        ├── exampleLTVdisc.m             # LTV identification with COSMIC
        ├── exampleMultiTrajectory.m     # Multi-trajectory ensemble averaging
        └── exampleOutputCOSMIC.m        # Partial-observation identification
```

## Installation

Unpack the release archive, then add sid to the MATLAB/Octave path:

```matlab
run('/path/to/sid-v0.1.0-matlab/matlab/sidInstall.m')
```

Or clone the full repository for update access via `git pull`:

```bash
git clone https://github.com/pdlourenco/sid.git
```

```matlab
run('/path/to/sid/matlab/sidInstall.m')
```

To make the path persistent, add the `sidInstall` line to your `startup.m`
(MATLAB) or `.octaverc` (Octave).

## Quick Start

```matlab
% SISO frequency response estimation
N = 1000; Ts = 0.01;
u = randn(N, 1);
y = filter([1], [1 -0.9], u) + 0.1 * randn(N, 1);
result = sidFreqBT(y, u, 'SampleTime', Ts);
sidBodePlot(result);

% Time series (output spectrum only)
result_ts = sidFreqBT(y, []);
sidSpectrumPlot(result_ts);
```

See the [examples guide](examples/README.md) for ETFE, frequency-dependent
resolution, MIMO, time-varying maps, spectrograms, LTV state-space
identification, multi-trajectory ensemble averaging, and partial-observation
Output-COSMIC.

## Function Reference

### Frequency-Domain Estimation

| Function | Replaces | Description |
|---|---|---|
| `sidFreqBT` | `spa` | Blackman-Tukey frequency response and noise spectrum estimation. Supports SISO, MIMO, time series, and multi-trajectory inputs (3D arrays and cell arrays). |
| `sidFreqBTFDR` | `spafdr` | Blackman-Tukey with frequency-dependent resolution. Adapts the correlation window length across frequency. |
| `sidFreqETFE` | `etfe` | Empirical transfer function estimate via FFT ratio, with optional smoothing. |

### Time-Frequency Analysis

| Function | Replaces | Description |
|---|---|---|
| `sidFreqMap` | `tfestimate`, `mscohere`, `cpsd` | Sliding-window time-varying frequency response map using Blackman-Tukey or Welch windowing. Produces time-varying transfer function, coherence, and noise spectrum. |
| `sidSpectrogram` | `spectrogram` | Short-time FFT spectrogram with configurable window, overlap, and FFT length. |

### State-Space Identification (COSMIC)

| Function | Description |
|---|---|
| `sidLTVdisc` | COSMIC algorithm for discrete LTV system identification from full-state measurements. Identifies A(k), B(k) from single or multiple trajectories (including variable-length via cell arrays) with time-varying or constant regularisation λ. Returns system matrices and optional Bayesian posterior uncertainty. |
| `sidLTVdiscIO` | Partial-observation LTV identification (Output-COSMIC). Alternates between state estimation (RTS smoother) and dynamics identification (COSMIC) to jointly recover A(k), B(k), and hidden states. Includes fast path when H is full-rank. |
| `sidLTVdiscTune` | Regularisation tuning via validation-based grid search, L-curve, or frequency-response consistency scoring (frozen transfer function vs. `sidFreqMap` with Mahalanobis distance). |
| `sidLTVdiscFrozen` | Frozen transfer function G(ω, k) from identified A(k), B(k) with uncertainty propagated from the COSMIC posterior, enabling direct comparison with non-parametric frequency estimates. |
| `sidLTIfreqIO` | LTI realisation from I/O frequency response via Ho-Kalman algorithm with H-basis transform. Used by `sidLTVdiscIO` for initialisation; also usable standalone. |
| `sidLTVStateEst` | Batch LTV state estimation via RTS smoother given A(k), B(k), H, R, and optional Q. Usable standalone or as a building block for `sidLTVdiscIO`. |
| `sidModelOrder` | Model order estimation from any `sidFreq*` result via Hankel SVD. Supports automatic gap detection and user-specified threshold. |

### Analysis and Validation

| Function | Replaces | Description |
|---|---|---|
| `sidDetrend` | `detrend` | Polynomial detrending. Removes constant, linear, or higher-order trends. Multi-channel support. |
| `sidResidual` | `resid` | Residual analysis with whiteness (autocorrelation) and independence (cross-correlation with input) tests, with confidence bounds. Works with both frequency-domain and state-space models. |
| `sidCompare` | `compare` | Model output simulation and comparison. Computes NRMSE fit metric. Supports frequency-domain and state-space models, multi-channel, and multi-trajectory data. |

### Plotting

| Function | Description |
|---|---|
| `sidBodePlot` | Bode diagram with shaded confidence bands. |
| `sidSpectrumPlot` | Power spectrum plot with shaded confidence bands. |
| `sidMapPlot` | Time-frequency colour map for `sidFreqMap` results. |
| `sidSpectrogramPlot` | Spectrogram colour map for `sidSpectrogram` results. |

### Result Types

Every estimation function returns a struct with documented fields. Run
`help sidResultTypes` for the complete reference.

| Result type | Produced by | Key fields |
|---|---|---|
| FreqResult | `sidFreqBT`, `sidFreqBTFDR`, `sidFreqETFE` | `.Response`, `.NoiseSpectrum`, `.Coherence`, `.ResponseStd` |
| FreqMapResult | `sidFreqMap` | `.Time`, `.Response`, `.NoiseSpectrum`, `.Coherence` |
| SpectrogramResult | `sidSpectrogram` | `.Time`, `.Frequency`, `.Power`, `.PowerDB`, `.Complex` |
| LTVResult | `sidLTVdisc`, `sidLTVdiscTune` | `.A`, `.B`, `.Lambda`, `.Cost`, `.AStd`, `.BStd`, `.P` |
| LTVIOResult | `sidLTVdiscIO` | `.A`, `.B`, `.X`, `.H`, `.R`, `.Cost`, `.Iterations` |
| FrozenResult | `sidLTVdiscFrozen` | `.Response`, `.ResponseStd`, `.TimeSteps` |
| CompareResult | `sidCompare` | `.Predicted`, `.Measured`, `.Fit`, `.Residual` |
| ResidualResult | `sidResidual` | `.Residual`, `.AutoCorr`, `.CrossCorr`, `.WhitenessPass` |

## Key Features

- **Closed-form, O(N) LTV identification.** COSMIC solves the regularised
  least-squares problem in a single forward-backward pass — no iterations,
  no convergence tuning. Complexity is linear in the number of time steps
  and cubic in the state dimension, independent of the number of
  trajectories.

- **Partial-observation identification.** Output-COSMIC extends the
  algorithm to systems where only a subset of states is measured. An
  alternating minimisation between dynamics identification and state
  estimation converges to a joint optimum, with a fast path that avoids
  iteration entirely when the observation matrix is full-rank.

- **Multi-trajectory support.** Frequency-domain functions accept 3D arrays
  (N × n_ch × L) for ensemble averaging with 1/L variance reduction. COSMIC
  pools multiple trajectories sharing the same time-varying dynamics, with
  variable-length trajectory support via cell arrays.

- **Bayesian uncertainty quantification.** Per-timestep posterior covariance
  of the identified system matrices via a backward recursion at O(N) cost.
  The posterior has matrix-normal structure and the MAP estimate is
  independent of the noise covariance. Uncertainty propagates through to
  frozen transfer functions for comparison with non-parametric frequency
  estimates.

- **No toolbox dependencies.** Drop-in replacements for `spa`, `spafdr`,
  `etfe`, `resid`, `compare`, and `spectrogram` using only core
  MATLAB/Octave functions.

## Compatibility

| Platform | Version | Status |
|---|---|---|
| MATLAB | R2016b or later | Tested in CI |
| GNU Octave | 8.0 or later | Tested in CI |

## Out of Scope for this version

- Online/recursive COSMIC (deferred to a future release)
- Parametric identification: ARX, ARMAX, state-space subspace methods
- LPV identification
- Unknown or time-varying observation matrix H
- Alternative regularisation norms
- Python and Julia implementations — Python ships as a separate v0.1.0
  release (`v0.1-python`); Julia is planned for a later version.

## References

- Ljung, L. (1999). *System Identification: Theory for the User*, 2nd ed.
  Prentice Hall.
- Blackman, R. B. & Tukey, J. W. (1959). *The Measurement of Power Spectra*.
  Dover.
- Carvalho, M., Soares, C., Lourenco, P., and Ventura, R. (2022). "COSMIC:
  fast closed-form identification from large-scale data for LTV systems."
  [arXiv:2112.04355v2](https://arxiv.org/abs/2112.04355v2)
- Laszkiewicz, P., Carvalho, M., Soares, C., and Lourenco, P. (2025). "The
  impact of modeling approaches on controlling safety-critical, highly
  perturbed systems: the case for data-driven models."
  [arXiv:2509.13531](https://arxiv.org/abs/2509.13531)
