# Examples

Every public estimator ships with a runnable example in both Python and
MATLAB. The two sets share a common cast of physical spring-mass-damper
plants defined in [`spec/EXAMPLES.md`](../spec/examples-spec.md) so the
same identification problem can be tackled from either language.

## Python notebooks

[Browse the notebook catalog](python/index.md) or pick one of the headline
examples:

- [SISO frequency response](python/example_siso.ipynb) — Blackman-Tukey on a 1-DoF SMD
- [LTV identification (COSMIC)](python/example_ltv_disc.ipynb) — time-varying state-space
- [Output-only COSMIC](python/example_output_cosmic.ipynb) — partial observations
- [Time-varying frequency map](python/example_freq_map.ipynb) — sliding-window analysis

Every notebook is executed at build time so the rendered outputs always
reflect the current codebase.

## MATLAB scripts

[Browse the MATLAB example catalog](matlab/index.md). The scripts are
runnable as-is from a MATLAB or Octave session after adding `matlab/sid/`
to the path.
