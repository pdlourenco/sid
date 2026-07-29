# Frequency-domain estimation

The frequency-domain estimators in `sid` are open-source replacements for
the System Identification Toolbox routines `spa`, `spafdr`, and `etfe`.

- [`freq_bt`](../api/python/freq_bt.md) / [`sidFreqBT`](../api/matlab/sidFreqBT.md) — Blackman-Tukey
- [`freq_btfdr`](../api/python/freq_btfdr.md) / [`sidFreqBTFDR`](../api/matlab/sidFreqBTFDR.md) — Frequency-dependent resolution
- [`freq_etfe`](../api/python/freq_etfe.md) / [`sidFreqETFE`](../api/matlab/sidFreqETFE.md) — Empirical transfer function estimate

The full mathematical derivation lives in the [Specification](../spec/index.md),
sections 2 (Blackman-Tukey), 3 (Uncertainty), and 4 (ETFE).
