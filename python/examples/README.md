# sid Examples

Runnable Jupyter notebooks demonstrating the full functionality of the sid
Python toolbox. Each notebook mirrors a MATLAB example script from
`matlab/examples/`, combining narrative, code, and inline plots.

## Example Index

| Notebook | Description | Functions |
|---|---|---|
| [`example_siso`](example_siso.ipynb) | Basic SISO frequency response estimation | `freq_bt`, `bode_plot`, `spectrum_plot` |
| [`example_etfe`](example_etfe.ipynb) | Empirical transfer function estimate | `freq_etfe`, `bode_plot`, `spectrum_plot` |
| [`example_freq_dep_res`](example_freq_dep_res.ipynb) | Frequency-dependent resolution | `freq_btfdr`, `freq_bt`, `bode_plot` |
| [`example_coherence`](example_coherence.ipynb) | Coherence analysis and signal quality | `freq_bt`, `bode_plot` |
| [`example_method_comparison`](example_method_comparison.ipynb) | Comparing BT, BTFDR, and ETFE | `freq_bt`, `freq_btfdr`, `freq_etfe` |
| [`example_mimo`](example_mimo.ipynb) | Multi-input multi-output systems | `freq_bt` (MIMO mode) |
| [`example_freq_map`](example_freq_map.ipynb) | Time-varying frequency response maps | `freq_map`, `map_plot` |
| [`example_spectrogram`](example_spectrogram.ipynb) | Short-time FFT spectrogram | `spectrogram`, `spectrogram_plot` |
| [`example_ltv_disc`](example_ltv_disc.ipynb) | LTV state-space identification | `ltv_disc`, `ltv_disc_tune`, `ltv_disc_frozen` |
| [`example_multi_trajectory`](example_multi_trajectory.ipynb) | Multi-trajectory ensemble averaging | `freq_bt`, `freq_map`, `spectrogram`, `ltv_disc` |
| [`example_output_cosmic`](example_output_cosmic.ipynb) | LTV identification from partial observations | `freq_bt`, `model_order`, `ltv_disc_io` |

## Running

Launch any notebook locally:

```bash
jupyter notebook python/examples/example_siso.ipynb
```

Or run all notebooks non-interactively (used in CI):

```bash
pytest --nbmake python/examples/ -v
```

## Contributing

See [`python/CONTRIBUTING.md`](../CONTRIBUTING.md) for notebook conventions:
naming, structure, auto-discovery, and CI validation.
