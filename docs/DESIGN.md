# sid — Design Document

> **Positioning.** This is the *rationale* half of the DESIGN / SPEC pair. It
> answers "why is sid shaped this way?" — architecture choices, the algorithmic
> approaches, the trade-offs accepted. The *contract* half lives in
> [`spec/SPEC.md`](../spec/SPEC.md): the formulas, defaults, output fields, and
> cross-language conventions every implementation must satisfy. When in doubt,
> link to the spec rather than duplicate it — two sources of truth are zero
> sources of truth. Tactical decisions beneath this architecture live in
> [`docs/decisions/`](decisions/).

## Overview

sid is an open-source system identification toolbox covering two families of
methods behind one interface:

- a **frequency-domain path** — non-parametric estimation of a system's transfer
  function and noise spectrum (`sidFreqBT`, `sidFreqETFE`, `sidFreqBTFDR`), with
  time-varying extensions (`sidFreqMap`, `sidSpectrogram`); and
- a **state-space path** — parametric identification of discrete linear
  time-varying dynamics `A(k), B(k)` via the COSMIC algorithm (`sidLTVdisc`),
  including partial-observation (`sidLTVdiscIO`) and uncertainty quantification.

It ships as independent per-language implementations (MATLAB/Octave stable,
Python stable, Julia planned) that share a single mathematical specification and
a set of cross-language reference vectors.

## The key architectural principle

**The interface between the language implementations is the written
specification in `spec/`, not any one implementation's code.**

Every port derives its behaviour independently from `spec/SPEC.md` and
`spec/EXAMPLES.md`; cross-language numerical equivalence is a *consequence* of
independent conformance, not a goal pursued by copying one port into another.
This is the load-bearing decision the whole project hangs on, recorded in
[ADR-0001](decisions/ADR-0001-spec-is-the-contract.md). It is what lets multiple
contributors (and agents) develop ports in parallel without silent drift, and
what makes the `testdata/` vectors meaningful: when two *independent*
derivations agree, that agreement is evidence; when one is copied from the
other, it is not. The cost is real — every normative behaviour must be written
into the spec before it is coded — and accepted deliberately.

## Why this shape

**Polyglot, contract-on-disk.** Users live in different ecosystems (MATLAB for
control engineering, Python for ML/data science, Julia emerging for scientific
computing). Rather than bind them to one runtime or ship thin wrappers around a
single core, sid implements each port natively in idiomatic style. The shared
artifact that keeps them honest is text (the spec) and data (JSON reference
vectors) — both language-agnostic by construction. A binary or pickle core would
privilege one language and reintroduce the single-implementation coupling the
spec exists to avoid.

**Auto-discovered tests and examples.** Runners discover files by naming
convention rather than a manifest, so the common failure mode — a test that
exists but is never run because nobody added it to a list — cannot happen. See
[`CONTRIBUTING.md`](../CONTRIBUTING.md) §"Test and Example Auto-Discovery".

## Why these algorithms

**Frequency path — Blackman-Tukey at the core.** The primary spectral estimators
use the correlogram (biased covariance → Hann lag window → Fourier transform)
rather than raw periodograms. The biased covariance estimator guarantees a
non-negative spectrum and has lower mean-squared error; the lag window gives a
single, interpretable resolution/variance knob (`M`). Asymptotic variance
formulas (Ljung 1999) provide per-frequency uncertainty in closed form. ETFE is
offered as the maximum-resolution / high-variance extreme (rectangular window,
`M = N-1`), and BTFDR as the frequency-varying-resolution variant. Welch is
available inside `sidFreqMap` for `tfestimate` compatibility. Rationale and the
bias/variance trade-offs are detailed per estimator in
[`spec/SPEC.md`](../spec/SPEC.md) §2–§7.

**State-space path — COSMIC's O(N) closed form.** Identifying time-varying
`A(k), B(k)` is posed as a regularized least-squares problem with a temporal
smoothness penalty. Setting the gradient to zero yields a *block-tridiagonal*
system solved by a forward-backward pass in `O(N(p+q)³)` — linear in horizon,
independent of the number of trajectories — instead of the `O(N³)` a dense solve
would cost. The same structure carries the Bayesian posterior covariance (via
left/right Schur complements, a Kalman-smoother analogue) and a naturally causal
forward pass that a future online variant can exploit. Multiple trajectories
(including variable-length) are pooled into the data matrices rather than
averaged after the fact, preserving the estimator structure. Full derivations
live in [`spec/SPEC.md`](../spec/SPEC.md) §8 and `spec/cosmic/`.

**Bridging the two paths.** The frozen transfer function `G(ω,k)` computed from
COSMIC's `A(k), B(k)` lets the parametric model be compared directly against the
non-parametric `sidFreqMap` estimate — including using the non-parametric
uncertainty bands to tune COSMIC's regularization (§8.11). The two paths are
deliberately designed to meet in the frequency domain rather than being
unrelated tools.

## Constraints

- **Permissive license only.** sid is MIT; dependencies must be license-compatible.
- **Runtime compatibility.** MATLAB R2016b+ and GNU Octave 8.0+ for the
  MATLAB/Octave port; Python 3.10+ with NumPy/SciPy for the Python port. Octave
  compatibility in particular constrains the MATLAB subset used (no
  newer-MATLAB-only syntax).
- **Dependency minimalism.** The numerical core leans on each ecosystem's
  standard array/linear-algebra facilities; heavy third-party dependencies are
  avoided so installation stays trivial and the ports stay comparable.
- **Numerical conventions are contract.** Normalization (System Identification
  Toolbox convention — no `Ts`, no `1/2π`), NaN/Inf substitutions, clamps, and
  warning identifiers are specified, not incidental — see `spec/SPEC.md` §2.8,
  §10, and the `Verified by:` annotations.

## Module map

The public function catalogue and the `sid + Domain + Method` naming convention
live in [`docs/roadmap.md`](roadmap.md); the binding behaviour of each function
lives in [`spec/SPEC.md`](../spec/SPEC.md). This document does not reproduce
either — see those for the authoritative lists. At a glance:

- **Frequency-domain estimators** — `spec/SPEC.md` §2–§7.
- **LTV state-space (COSMIC) and friends** — `spec/SPEC.md` §8 (+ `spec/cosmic/`).
- **Analysis / preprocessing** — detrending §13, residuals §14, comparison §15,
  model order §8.12.12.
- **Shared output struct** — §9 (one struct shape across the `sidFreq*` family).

## Future extensions

The deferred-with-conditions pattern, at the architectural level — what's out of
scope now, and what would bring it in. Contract-level deferrals are tracked in
`spec/SPEC.md` §8.14; this is the design-level view.

- **Parametric LTI identification** (`sidTfARX`, `sidTfARMAX`, `sidSsN4SID`,
  `sidTsAR`). The roadmap catalogue reserves the names; the `'Algorithm'`
  parameter on the LTV path is already shaped to admit alternatives. *Trigger:*
  demand for classical parametric methods alongside the LTV/non-parametric core.
- **Online / recursive COSMIC** (`sidLTVdiscInit/Update/Smooth`). The causal
  forward pass makes this natural; deferred to v2 (`spec/SPEC.md` §8.10).
  *Trigger:* a streaming/real-time use case.
- **Julia port.** The spec + reference vectors are designed to make a third port
  a matter of independent conformance, not redesign. *Trigger:* a maintainer or
  contributor to own it; implement the runner with auto-discovery from day one.
- **LPV identification.** Structured parameter-varying models via post-hoc
  regression on COSMIC output; design notes in `spec/lpv_extension_theory.md`.
