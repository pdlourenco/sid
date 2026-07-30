# sid — Open-Source System Identification

[![MATLAB/Octave Tests](https://github.com/pdlourenco/sid/actions/workflows/tests.yml/badge.svg)](https://github.com/pdlourenco/sid/actions/workflows/tests.yml)
[![Python Tests](https://github.com/pdlourenco/sid/actions/workflows/python-tests.yml/badge.svg)](https://github.com/pdlourenco/sid/actions/workflows/python-tests.yml)
[![Cross-Language Validation](https://github.com/pdlourenco/sid/actions/workflows/cross-validate.yml/badge.svg)](https://github.com/pdlourenco/sid/actions/workflows/cross-validate.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://github.com/pdlourenco/sid/blob/main/LICENSE)
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/pdlourenco/sid/main?labpath=python%2Fexamples)

**`sid`** is a free, open-source toolbox for **system identification** —
covering both **non-parametric frequency response estimation** and
**time-varying state-space identification**. All implementations share a
single mathematical specification and cross-language reference test
vectors to ensure numerical consistency.

## Get started

<div class="grid cards" markdown>

-   :material-language-python: **Python**

    ```bash
    pip install -e "./python[plot]"
    ```

    [Install guide :material-arrow-right:](getting-started/install-python.md)

-   :simple-mathworks: **MATLAB / Octave**

    ```matlab
    run('/path/to/sid/matlab/sidInstall.m')
    ```

    [Install guide :material-arrow-right:](getting-started/install-matlab.md)

</div>

## Features

- **Blackman-Tukey spectral analysis** — frequency response and noise spectrum estimation with configurable window size.
- **Frequency-dependent resolution** — vary the smoothing bandwidth across the frequency axis.
- **Empirical transfer function estimate** — maximum resolution via FFT ratio, with optional smoothing.
- **Time-varying frequency maps and spectrograms** — sliding-window analysis for non-stationary signals.
- **LTV state-space identification (COSMIC)** — identify time-varying A(k), B(k) with O(N) complexity, automatic regularization, and Bayesian uncertainty.
- **Partial-observation identification (Output-COSMIC)** — identify dynamics from output-only measurements.
- **Multi-trajectory support** — ensemble averaging for frequency estimates; pooled least-squares for state-space.
- **Asymptotic uncertainty estimates** — confidence bands for all estimation functions.
- **SISO, MIMO, and time-series modes** — unified API across all estimation functions.

## Two paths, one specification

`sid` provides a **frequency-domain** path (Blackman-Tukey, ETFE) and a
**time-domain** state-space path (COSMIC). The two are derived from a
single [mathematical specification](spec/index.md) and validated against
shared test vectors so that both implementations agree to floating-point
tolerance.

[Read the overview :material-arrow-right:](concepts/overview.md)
