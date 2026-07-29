# Function index

Every public function in `sid` exists in both Python and MATLAB/Octave with a
1:1 mapping. The tables below link to the per-function reference page in each
language.

## Frequency-domain estimation

| Python | MATLAB | Description |
|---|---|---|
| [`freq_bt`](python/freq_bt.md) | [`sidFreqBT`](matlab/sidFreqBT.md) | Blackman-Tukey frequency response and noise spectrum |
| [`freq_btfdr`](python/freq_btfdr.md) | [`sidFreqBTFDR`](matlab/sidFreqBTFDR.md) | Blackman-Tukey with frequency-dependent resolution |
| [`freq_etfe`](python/freq_etfe.md) | [`sidFreqETFE`](matlab/sidFreqETFE.md) | Empirical transfer function estimate (FFT ratio) |

## Time-frequency analysis

| Python | MATLAB | Description |
|---|---|---|
| [`freq_map`](python/freq_map.md) | [`sidFreqMap`](matlab/sidFreqMap.md) | Time-varying frequency response map |
| [`spectrogram`](python/spectrogram.md) | [`sidSpectrogram`](matlab/sidSpectrogram.md) | Short-time FFT spectrogram |

## State-space identification (COSMIC)

| Python | MATLAB | Description |
|---|---|---|
| [`ltv_disc`](python/ltv_disc.md) | [`sidLTVdisc`](matlab/sidLTVdisc.md) | Discrete LTV state-space identification |
| [`ltv_disc_io`](python/ltv_disc_io.md) | [`sidLTVdiscIO`](matlab/sidLTVdiscIO.md) | LTV identification from partial observations |
| [`ltv_disc_tune`](python/ltv_disc_tune.md) | [`sidLTVdiscTune`](matlab/sidLTVdiscTune.md) | Regularization tuning |
| [`ltv_disc_frozen`](python/ltv_disc_frozen.md) | [`sidLTVdiscFrozen`](matlab/sidLTVdiscFrozen.md) | Frozen-time transfer function G(ω, k) |
| [`lti_freq_io`](python/lti_freq_io.md) | [`sidLTIfreqIO`](matlab/sidLTIfreqIO.md) | LTI realization from I/O frequency response |
| [`ltv_state_est`](python/ltv_state_est.md) | [`sidLTVStateEst`](matlab/sidLTVStateEst.md) | Batch LTV state estimation (RTS smoother) |
| [`model_order`](python/model_order.md) | [`sidModelOrder`](matlab/sidModelOrder.md) | Model order estimation via Hankel SVD |

## Plotting

| Python | MATLAB | Description |
|---|---|---|
| [`bode_plot`](python/bode_plot.md) | [`sidBodePlot`](matlab/sidBodePlot.md) | Bode diagram with confidence bands |
| [`spectrum_plot`](python/spectrum_plot.md) | [`sidSpectrumPlot`](matlab/sidSpectrumPlot.md) | Power spectrum plot |
| [`map_plot`](python/map_plot.md) | [`sidMapPlot`](matlab/sidMapPlot.md) | Time-frequency map plot |
| [`spectrogram_plot`](python/spectrogram_plot.md) | [`sidSpectrogramPlot`](matlab/sidSpectrogramPlot.md) | Spectrogram color map |

## Analysis and preprocessing

| Python | MATLAB | Description |
|---|---|---|
| [`detrend`](python/detrend.md) | [`sidDetrend`](matlab/sidDetrend.md) | Remove mean / linear / polynomial trends |
| [`residual`](python/residual.md) | [`sidResidual`](matlab/sidResidual.md) | Residual analysis with auto/cross-correlation |
| [`compare`](python/compare.md) | [`sidCompare`](matlab/sidCompare.md) | Simulate model output and compute fit percentage |

## Result types

Every function returns a frozen dataclass (Python) or struct (MATLAB).
Browse the [Python result types](python/results.md) page for field-level
documentation; MATLAB struct layouts are documented inline on each function
page.
