# sid Examples

Runnable Jupyter notebooks demonstrating the full functionality of the sid
Python toolbox. Each notebook mirrors a MATLAB example script from
`matlab/examples/`, combining narrative, code, and inline plots.

## Why spring-mass-damper plants?

Every example in this directory is built around a **spring-mass-damper (SMD)
chain** — one, two, or three lumped masses connected by linear springs and
viscous dampers to a wall. The discrete-time state-space model is built with
`util_msd` (exact zero-order-hold via matrix exponential), its time-varying
variant `util_msd_ltv`, or the nonlinear Duffing simulator `util_msd_nl`, all
defined in the sibling module [`util_msd.py`](util_msd.py).

Using a single physical plant family across the tutorials keeps the narrative
concrete and cumulative. Instead of arbitrary `lfilter` coefficients or
hand-crafted state-space matrices, the reader works through a progression of
plants they can picture: a 1-DoF oscillator, a 2-mass chain with cross-coupling
through a shared spring, a 3-mass chain with three closely-spaced modes, a
structure whose stiffness drifts thermally, a structure that snaps, and
finally a Duffing-style nonlinear oscillator whose apparent resonance depends
on amplitude.

## Example Index

| Notebook | Plant | Functions demonstrated |
|---|---|---|
| [`example_siso`](example_siso.ipynb) | 1-DoF SMD (`m=1, k=100, c=2`, `ω_n=10 rad/s`, `ζ=0.1`) | `freq_bt`, `bode_plot`, `spectrum_plot`, `detrend`, `residual`, `compare` |
| [`example_etfe`](example_etfe.ipynb) | Same 1-DoF SMD as siso | `freq_etfe`, `bode_plot`, `spectrum_plot` |
| [`example_freq_dep_res`](example_freq_dep_res.ipynb) | 3-mass chain (`k=[300,200,100]`, modes at 6.4/15.1/25.1 rad/s) | `freq_bt`, `freq_btfdr`, `bode_plot` |
| [`example_coherence`](example_coherence.ipynb) | 2-mass chain with colored force disturbance at mass 2 | `freq_bt`, `bode_plot` (with coherence) |
| [`example_method_comparison`](example_method_comparison.ipynb) | Same 1-DoF SMD as siso | `freq_bt`, `freq_btfdr`, `freq_etfe`, `compare` |
| [`example_mimo`](example_mimo.ipynb) | 2-mass chain, 2 force inputs and 2 position outputs | `freq_bt` (MIMO mode) |
| [`example_freq_map`](example_freq_map.ipynb) | 2-mass LTI + continuous LTV (`k₁` ramp) + discrete LTV (step change) + Duffing hardening SDOF | `freq_map`, `map_plot` |
| [`example_spectrogram`](example_spectrogram.ipynb) | SDOF (`ω_n ≈ 32 Hz`, `Fs=1000 Hz`) driven by a 20→60 Hz chirp force | `spectrogram`, `spectrogram_plot` |
| [`example_ltv_disc`](example_ltv_disc.ipynb) | 1-DoF LTV (ramping stiffness) **plus** Duffing linearization | `ltv_disc`, `ltv_disc_tune`, `ltv_disc_frozen`, `compare`, `residual` |
| [`example_multi_trajectory`](example_multi_trajectory.ipynb) | 2-mass ensemble + step-change LTV + chirp-response spectrogram + LTV/freq_map consistency | `freq_bt`, `freq_map`, `spectrogram`, `ltv_disc` (all ensemble mode) |
| [`example_output_cosmic`](example_output_cosmic.ipynb) | 2-mass chain, position-only measurements (velocities hidden) | `freq_bt`, `model_order`, `ltv_disc_io`, `compare`, `residual` |

## Running

Launch any notebook locally (start Jupyter from the `python/examples/`
directory so the sibling `util_msd` import resolves):

```bash
cd python/examples
jupyter notebook example_siso.ipynb
```

Or run all notebooks non-interactively (used in CI):

```bash
pytest --nbmake python/examples/ -v
```

The notebooks import `util_msd` as a sibling module rather than from
`sid._internal`; the `pytest --nbmake` runner sets the working directory to
`python/examples/` automatically, and `jupyter notebook python/examples/...`
does the same. If you start Jupyter elsewhere, either `cd` into
`python/examples/` first or add that directory to `sys.path` before the
import cell.

## Contributing

See [`python/CONTRIBUTING.md`](../CONTRIBUTING.md) for notebook conventions:
naming, structure, auto-discovery, and CI validation.
