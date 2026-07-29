# Python API

The `sid` package exposes 19 public functions and 8 result dataclasses from
its top-level namespace. Every function maps 1:1 to a MATLAB equivalent —
see the [function index](../index.md) for the pairing table.

## Frequency-domain estimators

- [`freq_bt`](freq_bt.md) — Blackman-Tukey spectral analysis
- [`freq_btfdr`](freq_btfdr.md) — Blackman-Tukey with frequency-dependent resolution
- [`freq_etfe`](freq_etfe.md) — Empirical transfer function estimate
- [`freq_map`](freq_map.md) — Time-varying frequency map
- [`spectrogram`](spectrogram.md) — Sliding-window spectrogram

## State-space identification

- [`ltv_disc`](ltv_disc.md) — COSMIC LTV state-space identification
- [`ltv_disc_io`](ltv_disc_io.md) — Output-COSMIC (output-only identification)
- [`ltv_disc_tune`](ltv_disc_tune.md) — Automatic regularization tuning
- [`ltv_disc_frozen`](ltv_disc_frozen.md) — Frozen-time transfer functions G(ω, k)
- [`lti_freq_io`](lti_freq_io.md) — LTI frequency response from input/output data
- [`ltv_state_est`](ltv_state_est.md) — RTS-smoother state estimation
- [`model_order`](model_order.md) — Model order selection

## Analysis

- [`compare`](compare.md) — Model comparison metrics
- [`detrend`](detrend.md) — Detrending utilities
- [`residual`](residual.md) — Residual diagnostics

## Plotting

- [`bode_plot`](bode_plot.md) — Bode plots from `FreqResult`
- [`spectrum_plot`](spectrum_plot.md) — Power spectrum plots
- [`map_plot`](map_plot.md) — Time-frequency map plots
- [`spectrogram_plot`](spectrogram_plot.md) — Spectrogram plots

## Result types

- [Result dataclasses](results.md) — `FreqResult`, `LTVResult`, etc.
