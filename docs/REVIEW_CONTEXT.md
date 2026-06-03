# sid — Reviewer Context

**Purpose.** This document seeds a reviewer (human or agent) with the context to
review changes to sid with judgment, not just surface-level linting. Point the
reviewer at this file alongside [`spec/SPEC.md`](../spec/SPEC.md),
[`spec/EXAMPLES.md`](../spec/EXAMPLES.md), and the relevant per-language
`CONTRIBUTING.md` before a review. It does double duty: a reviewer prompt **and**
a codified statement of what the project cares about — new contributors benefit
from reading it too.

Key documents:

- [`spec/SPEC.md`](../spec/SPEC.md) — the binding algorithm contract (with
  per-rule `Verified by:` mechanisms).
- [`spec/EXAMPLES.md`](../spec/EXAMPLES.md) — the example-suite contract.
- [`CONTRIBUTING.md`](../CONTRIBUTING.md) — spec discipline, workflow, CI,
  auto-discovery; per-language guides under `matlab/` and `python/`.
- [`docs/decisions/`](decisions/) — ADRs (the *why* of tactical choices).
- [`docs/roadmap.md`](roadmap.md) — function catalogue and naming convention.

## Verification vs validation

Two review modes; a reviewer can run either or both (default both):

- **Verification — *did we build it right?*** Does the diff match the binding
  contracts in `spec/SPEC.md` / `spec/EXAMPLES.md`, the `Verified by:` mechanisms,
  the function catalogue in `docs/roadmap.md`, the test/example naming
  conventions? Checks against **named artifacts**; findings are mechanical (a
  rule said X, the diff did Y).
- **Validation — *did we build the right thing?*** Does the diff match the
  principles below, the project goals, and the PR's stated scope? Checks against
  **intent**; findings are judgment calls.

## Project in one paragraph

sid is a multi-language system identification toolbox (MATLAB/Octave, Python,
Julia planned) built around a single mathematical specification and
cross-language reference vectors. The architecture: independent per-language
ports that each conform to the shared contract in `spec/`, developed largely by
coding agents working in parallel. Locked-in invariants: the spec is the sole
contract (ADR-0001); ports derive their behaviour independently; numerical
equivalence is a *consequence* of independent conformance, checked — not proven
— by the `testdata/` vectors. The goal that shapes everything: correctness you
can trust across languages because two independent derivations agree, not
because one was copied from the other.

## Core principles (review against these, cite by number)

1. **The spec is the contract.** `spec/SPEC.md` and `spec/EXAMPLES.md` define
   required behaviour — defaults, edge cases, error conditions, output fields,
   normalization, numerical diagnostics — not just the core formulas.
   Implementations conform to the spec, **not to each other** (ADR-0001).
2. **MATLAB is not ground truth.** "The MATLAB version does it this way" is not a
   justification. If MATLAB and the spec disagree, MATLAB is wrong.
3. **Spec first.** New or changed behaviour updates `spec/SPEC.md` *before* the
   code, in the same PR. Never silently encode behaviour the spec doesn't
   describe — the next port can't recover the decision.
4. **Shared helpers cause shared drift.** A bug in `sidValidateData` /
   `validate_data` (or any shared helper) makes *every* caller violate the spec
   identically, which cross-validation will **not** catch. Touching a shared
   helper requires auditing each caller against the relevant spec section.
5. **Reference vectors are a check, not a proof.** `testdata/` verifies the ports
   agree on sampled inputs; it does not prove either satisfies the spec. New
   behaviour needs a direct spec-to-implementation read-through, not only a
   vector comparison. A test of the form "assert it equals what MATLAB returned
   today" does not detect joint drift.
6. **Every normative rule names its verifier.** Per `spec/SPEC.md` §"Verification
   (right-side mechanisms)", each rule carries a `Verified by:` mechanism. A new
   or changed rule with `Verified by: none` is visible debt and should be flagged
   (or the gap tracked as an issue).
7. **Tests and examples are auto-discovered.** New tests/examples must match the
   naming convention (`test_*`, `example*`/`example_*`) so runners pick them up;
   never maintain a hardcoded manifest.
8. **Numerical diagnostics are contract, not incident.** NaN/Inf substitutions,
   clamps, and ill-conditioning warnings are specified behaviour with stable
   warning identifiers — changing them is a contract change.

## Terminology (enforce consistency)

- **LTI vs LTV** — time-invariant vs time-varying; `sidFreqBT`/`sidFreqETFE` are
  LTI, `sidFreqMap`/`sidLTVdisc` are time-varying. Don't conflate.
- **Trajectory vs segment** — a *trajectory* is an independent experiment
  (multi-trajectory = 3D arrays / cell arrays); a *segment* is a sliding window
  within one record (`sidFreqMap`, `sidSpectrogram`).
- **Lag window vs time-domain window** — the BT *lag window* (Hann over
  covariance lags, controls frequency resolution) is distinct from the
  spectrogram/Welch *time-domain window* (tapers each segment, reduces leakage).
- **Resolution vs window size** — `sidFreqBTFDR` takes a *resolution* `R(ω)`;
  window size `M_k = ceil(2π/R_k)` is derived.
- **COSMIC** — the closed-form LTV identification algorithm (§8); "Output-COSMIC"
  is the partial-observation variant (`sidLTVdiscIO`).
- **Frozen transfer function** — `G(ω,k)` from frozen `A(k),B(k)`, not a
  time-frequency map.

## Common red flags

- **Contract drift** — a port changed to match another port rather than the spec;
  a spec rule and an implementation disagreeing with no indication which is
  right (stop and ask, don't pick a side).
- **Copy-by-copy porting** — Python logic transcribed from MATLAB (or vice versa)
  instead of derived from the spec section.
- **Vector-only justification** — a new feature whose only evidence is "the
  cross-validation vectors match", with no spec read-through.
- **Unaudited shared-helper edit** — touching `validate_data` / `sidValidateData`
  / covariance / windowing helpers without checking every caller.
- **Spec-silent behaviour** — a new default, bound, threshold, or warning in code
  that isn't in `spec/SPEC.md`; or a new spec rule with no `Verified by:`.
- **Weakening the gate** — relaxing tolerances or skipping cases in
  `cross-validate` to make a PR green.
- **Catalogue drift** — a new public function whose name doesn't follow the
  `sid + Domain + Method` convention or isn't reflected in `docs/roadmap.md`.

## What to be lenient about

- Naming bikesheds, unless they break consistency with existing terminology.
- Imperfect prose in rationale docs (DESIGN, ADRs) — they aren't contracts.
- Minor stylistic inconsistencies across documents when the substance is right.
- Missing tests on clearly throwaway / exploratory scripts.

## What to be strict about

- Anything touching the cross-language contract (`spec/SPEC.md`,
  `spec/EXAMPLES.md`, `testdata/` formats).
- New normative behaviour landing without spec-first (principle 3).
- Terminology drift in the spec.
- Silently breaking or weakening cross-validation.
- A new magic number / threshold / fallback with no spec coverage and no ADR.

## Review output format

1. **Summary** — one paragraph: what the change does, is the direction right, is
   it ready to merge.
2. **What works well** — brief.
3. **Issues to address before merge** — numbered; each with `file:line`, severity
   (blocker / non-blocker), and a concrete suggested change. Cite the principle
   number or the affected spec rule / `Verified by:` mechanism.
4. **Follow-up suggestions** — non-blocking.
5. **Verdict** — approve / approve with changes / request changes (with the
   conditional spelled out if "approve with changes"). If no mechanical checks
   ran (CI unreachable, no toolchain), say so — a review that didn't verify is a
   validation-only review and should label itself as such.

## Tone

This is a small / solo-maintainer project. Reviews are a conversation, not an
approval gate — "I'd do this differently but yours works too" is a legitimate
thing to say. Be direct about what's a blocker versus a nit. When you disagree
with the author, propose concrete replacement wording and the reason ("this
should read X because Y"), don't just flag.
