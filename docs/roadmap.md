# sid — function catalogue and roadmap

Language-neutral. This document defines the **public function catalogue** and the
**naming convention** every language port follows, and records what is planned,
reserved, and out of scope.

It deliberately says nothing about *behaviour*: [`spec/SPEC.md`](../spec/SPEC.md)
is the single source of truth for formulas, defaults, edge cases, output fields,
and error conditions, and [`spec/EXAMPLES.md`](../spec/EXAMPLES.md) for the
example-suite contract. Implementations conform to the spec, **not to each
other** ([ADR-0001](decisions/ADR-0001-spec-is-the-contract.md)).

Per-language coding standards, style, and testing conventions live in the
language guides: [`matlab/CONTRIBUTING.md`](../matlab/CONTRIBUTING.md) and
[`python/CONTRIBUTING.md`](../python/CONTRIBUTING.md).

## Naming convention

All public functions follow the pattern `sid` + `Domain` + `Method`:

```
sid  [Domain]  [Method/Variant]
 │      │          │
 │      │          └── BT, BTFDR, ETFE, ARX, N4SID, AR, ...
 │      └──────────── Freq, TF, SS, TS, LTV, ...
 └─────────────────── system identification (root)
```

**Per-language spelling.** MATLAB/Octave uses the camelCase form directly
(`sidFreqBT`). Python translates it to snake_case, dropping the `sid` prefix
(the package name supplies it) — and result fields translate the same way:

```
sidFreqBT(y, u, 'WindowSize', 30)   ->   sid.freq_bt(y, u, window_size=30)
result.Response                     ->   result.response
result.NoiseSpectrum                ->   result.noise_spectrum
```

Where a name collides with a language keyword, the port applies its own
escaping convention (e.g. MATLAB `Lambda` -> Python `lambda_`, PEP 8 trailing
underscore). Those rules belong to the language guide, not here.

**A new public function must fit this convention and appear in the catalogue
below.** A name that does neither is a review red flag — see
[`docs/REVIEW_CONTEXT.md`](REVIEW_CONTEXT.md) ("catalogue drift") and the
reviewer checklist in [`CONTRIBUTING.md`](../CONTRIBUTING.md). Adding or
renaming a public function is a major decision under `CLAUDE.md` §4.

## Function catalogue

Status is per language: ✅ implemented · ⬜ planned · — reserved (not started).
Julia is planned; no function is implemented there yet.

### Estimation

| Capability | MATLAB/Octave | Python | Julia |
|---|---|---|---|
| Frequency response via Blackman-Tukey | ✅ `sidFreqBT` | ✅ `freq_bt` | — |
| Blackman-Tukey, frequency-dependent resolution | ✅ `sidFreqBTFDR` | ✅ `freq_btfdr` | — |
| Empirical transfer function estimate | ✅ `sidFreqETFE` | ✅ `freq_etfe` | — |
| Time-varying frequency response map (BT or Welch) | ✅ `sidFreqMap` | ✅ `freq_map` | — |
| Short-time FFT spectrogram | ✅ `sidSpectrogram` | ✅ `spectrogram` | — |
| Discrete LTV state-space identification (COSMIC) | ✅ `sidLTVdisc` | ✅ `ltv_disc` | — |
| Regularisation tuning (validation, frequency) | ✅ `sidLTVdiscTune` | ✅ `ltv_disc_tune` | — |
| Frozen transfer function `G(ω,k)` | ✅ `sidLTVdiscFrozen` | ✅ `ltv_disc_frozen` | — |
| Partial-observation LTV identification (Output-COSMIC) | ✅ `sidLTVdiscIO` | ✅ `ltv_disc_io` | — |
| Batch LTV state estimation (RTS smoother) | ✅ `sidLTVStateEst` | ✅ `ltv_state_est` | — |
| LTI realisation from I/O frequency response (Ho-Kalman) | ✅ `sidLTIfreqIO` | ✅ `lti_freq_io` | — |
| Model order estimation (Hankel SVD) | ✅ `sidModelOrder` | ✅ `model_order` | — |

### Workflow utilities

| Capability | MATLAB/Octave | Python | Julia |
|---|---|---|---|
| Polynomial detrending | ✅ `sidDetrend` | ✅ `detrend` | — |
| Residual analysis (whiteness + independence) | ✅ `sidResidual` | ✅ `residual` | — |
| Model output comparison (NRMSE fit) | ✅ `sidCompare` | ✅ `compare` | — |

### Plotting

| Capability | MATLAB/Octave | Python | Julia |
|---|---|---|---|
| Bode diagram with confidence bands | ✅ `sidBodePlot` | ✅ `bode_plot` | — |
| Power spectrum with confidence bands | ✅ `sidSpectrumPlot` | ✅ `spectrum_plot` | — |
| Time-frequency colour map | ✅ `sidMapPlot` | ✅ `map_plot` | — |
| Spectrogram colour map | ✅ `sidSpectrogramPlot` | ✅ `spectrogram_plot` | — |
| Nyquist plot | — `sidNyquistPlot` | — `nyquist_plot` | — |
| Pole-zero map | — `sidPolePlot` | — `pole_plot` | — |

Private helpers are an implementation detail of each port, not part of this
catalogue: they are free to differ in decomposition, and pinning a
MATLAB↔Python helper mapping here would be an unverified drift surface. Each
port documents its own internals in its language guide.

## Planned

### Online/recursive COSMIC

**Theory:** [`spec/cosmic/online_recursion.md`](../spec/cosmic/online_recursion.md)

The COSMIC forward pass is a Kalman filter on parameter evolution and the
backward pass an RTS smoother, which admits a streaming formulation:

| Capability | MATLAB/Octave | Python | Julia |
|---|---|---|---|
| Initialise recursive/online estimator | ⬜ `sidLTVdiscInit` | ⬜ `ltv_disc_init` | — |
| Process one time step (filtered estimate) | ⬜ `sidLTVdiscUpdate` | ⬜ `ltv_disc_update` | — |
| Backward pass over window (smoothed estimates) | ⬜ `sidLTVdiscSmooth` | ⬜ `ltv_disc_smooth` | — |

Target properties: `O((p+q)³)` per step regardless of history length; filtered
estimates converging to smoothed as the window grows; filtered uncertainty ≥
smoothed uncertainty at every step.

### Reserved names

Not started, and reserved so the convention stays coherent if parametric
identification is taken on:

| Capability | MATLAB/Octave | Python | Replaces |
|---|---|---|---|
| Transfer function, ARX model | `sidTfARX` | `tf_arx` | `arx` |
| Transfer function, ARMAX model | `sidTfARMAX` | `tf_armax` | `armax` |
| State-space, N4SID subspace method | `sidSsN4SID` | `ss_n4sid` | `n4sid` |
| Time series, autoregressive | `sidTsAR` | `ts_ar` | `ar` |
| Time series, ARMA | `sidTsARMA` | `ts_arma` | `arma` |

## Out of scope for v1.0

Deferred, not rejected — these are excluded from the v1.0 surface, and none is
ruled out for a later version. The algorithmic subset of this list is normative
in [`spec/SPEC.md`](../spec/SPEC.md) §8.14 ("Deferred Extensions"); if the two
disagree, the spec wins. Items already earmarked for a specific future version
are marked **(v2)**.

- Online/recursive COSMIC — **(v2)**, see [Planned](#onlinerecursive-cosmic)
- Parametric identification: ARX, ARMAX, state-space subspace methods — **(v2)**
- LPV identification: structured parameter-varying models — **(v2)**
- Repository refactor for multi-language layout — **(v2)**
- Unknown observation matrix estimation (joint `H` + dynamics)
- Time-varying observation matrix `H(k)`
- Alternative regularisation norms (non-squared L2, L1 total variation)
- Alternative LTV algorithms (TVERA, TVOKID) — the `'Algorithm'` parameter is
  ready for them
- GCV lambda selection
- Frequency-domain input data
- Continuous-time models
- `maxSize` data segmentation
- Custom window functions (Hann only for `sidFreqBT`)
- idfrd-compatible class
- C reference implementation

## History

The v1.0 build-out phase logs are archived as historical records, not guidance:

- [`docs/plans/2026-07-29-matlab-v1.0-phase-log.md`](plans/2026-07-29-matlab-v1.0-phase-log.md)
- [`docs/plans/2026-07-29-python-port-phase-log.md`](plans/2026-07-29-python-port-phase-log.md)

The restructure that produced this document is recorded in
[ADR-0006](decisions/ADR-0006-roadmap-layer-split.md).
