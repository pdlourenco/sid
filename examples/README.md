# sid Examples

Runnable scripts demonstrating the full functionality of the sid toolbox. Each example is self-contained and can be run directly in MATLAB or Octave.

## Example Index

| Example | Description | Functions |
|---|---|---|
| [`exampleSISO`](exampleSISO.m) | Basic SISO frequency response estimation | `sidFreqBT`, `sidBodePlot`, `sidSpectrumPlot` |
| [`exampleETFE`](exampleETFE.m) | Empirical transfer function estimate | `sidFreqETFE`, `sidBodePlot`, `sidSpectrumPlot` |
| [`exampleFreqDepRes`](exampleFreqDepRes.m) | Frequency-dependent resolution | `sidFreqBTFDR`, `sidFreqBT`, `sidBodePlot` |
| [`exampleCoherence`](exampleCoherence.m) | Coherence analysis and signal quality | `sidFreqBT`, `sidBodePlot` |
| [`exampleMethodComparison`](exampleMethodComparison.m) | Comparing BT, BTFDR, and ETFE | `sidFreqBT`, `sidFreqBTFDR`, `sidFreqETFE` |
| [`exampleMIMO`](exampleMIMO.m) | Multi-input multi-output systems | `sidFreqBT` (MIMO mode) |
| [`exampleFreqMap`](exampleFreqMap.m) | Time-varying frequency response maps | `sidFreqMap`, `sidMapPlot` |
| [`exampleSpectrogram`](exampleSpectrogram.m) | Short-time FFT spectrogram | `sidSpectrogram`, `sidSpectrogramPlot` |
| [`exampleLTVdisc`](exampleLTVdisc.m) | LTV state-space identification | `sidLTVdisc`, `sidLTVdiscTune`, `sidLTVdiscFrozen` |

## Descriptions

### exampleSISO

Estimates the frequency response of a first-order SISO system using `sidFreqBT` (Blackman-Tukey). Demonstrates Bode diagram and noise spectrum plotting with confidence bands, the effect of window size on the bias-variance trade-off, and time-series mode for output spectrum estimation.

### exampleETFE

Introduces `sidFreqETFE`, which estimates the frequency response as the ratio of output and input DFTs. Shows the high-resolution but high-variance nature of the raw ETFE, how the `'Smoothing'` parameter reduces variance, exact recovery of a pure delay (noiseless FIR), and time-series periodogram mode. Also demonstrates custom frequency grids and Hz-unit plotting.

### exampleFreqDepRes

Demonstrates `sidFreqBTFDR`, which adapts the smoothing window size at each frequency. Uses a second-order resonant system (poles at `0.9*exp(+/-j*pi/4)`) to show why a fixed window struggles with sharp peaks. Compares scalar and per-frequency resolution vectors, and overlays the true system for validation.

### exampleCoherence

Shows how squared coherence quantifies estimation quality across frequencies. Uses an ARMA system with colored noise to produce frequency-dependent coherence. Demonstrates confidence band customization via the `'Confidence'` option, and compares high-noise vs low-noise scenarios.

### exampleMethodComparison

Head-to-head comparison of `sidFreqBT`, `sidFreqBTFDR`, and `sidFreqETFE` on the same data. Overlays all estimates on one plot with the true system response. Covers logarithmic frequency grids, noise spectrum comparison, and time-series spectrum estimation across methods.

### exampleMIMO

Demonstrates MIMO frequency response estimation with `sidFreqBT`. Shows the 3D `Response` array structure for a 2-output, 1-input system, manual per-channel Bode plotting, the Hermitian noise spectral matrix with positive semi-definiteness checking, and a full 2-output, 2-input transfer matrix. Notes the absence of coherence and uncertainty for MIMO.

### exampleFreqMap

Uses `sidFreqMap` to analyze how a system's frequency response evolves over time. Starts with an LTI baseline (flat map), then simulates a time-varying system with a drifting pole. Compares the Blackman-Tukey and Welch algorithms, explores segment length and overlap tuning, and demonstrates time-series spectrum maps. All visualizations use `sidMapPlot`.

### exampleSpectrogram

Demonstrates `sidSpectrogram` for time-frequency analysis of signals. Generates a chirp signal (frequency sweep) and visualizes it with `sidSpectrogramPlot`. Explores the window length trade-off between time and frequency resolution, compares Hann, Hamming, and rectangular windows, and shows multi-channel and log-frequency-scale plotting.

### exampleLTVdisc

Comprehensive demonstration of `sidLTVdisc` (COSMIC algorithm) for identifying time-varying state-space models `x(k+1) = A(k)x(k) + B(k)u(k)`. Covers LTI system recovery, LTV identification with a ramping pole, automatic vs manual regularization (`'Lambda'`), multi-trajectory benefits, and validation-based tuning via `sidLTVdiscTune`. Also demonstrates Bayesian uncertainty quantification (`'Uncertainty', true`) with confidence bands on recovered parameters, and frozen transfer function computation via `sidLTVdiscFrozen` with propagated uncertainty.

## Running All Examples

To run all examples and verify they execute without error:

```matlab
run('examples/runAllExamples.m')
```

This runner script is also used in CI to validate examples on both MATLAB and GNU Octave.
