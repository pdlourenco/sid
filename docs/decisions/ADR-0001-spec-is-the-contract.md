# ADR-0001: Implementations conform to the spec, not to each other

- **Status:** Accepted — _2026-06-02_
- **Deciders:** Project maintainer
- **Related:** [`spec/SPEC.md`](../../spec/SPEC.md), [`spec/EXAMPLES.md`](../../spec/EXAMPLES.md), [`CONTRIBUTING.md`](../../CONTRIBUTING.md) §"Specification as Source of Truth", [`CLAUDE.md`](../../CLAUDE.md) §3

> Retroactive ADR: this records a principle the project has operated under since
> its early multi-language structure. It is written now as the first entry in
> the decision trail so the rationale is recoverable, not because the decision
> is new.

## Context

sid is a multi-language toolbox (MATLAB/Octave, Python, Julia planned) where
independent language ports are developed largely in parallel by coding agents.
The central risk in that model is **silent drift**: two ports that disagree, or
— more insidiously — two ports that *agree* because one was copied from the
other and both inherited the same error.

If one port (historically MATLAB, the first written) is treated as the ground
truth that others are ported from, then a bug in that port propagates to every
other port, and the cross-language reference vectors in `testdata/` — which only
check that the ports agree — will pass green while every port is wrong. The
vectors cannot catch joint drift; only an independent derivation from a single
written specification can.

## Decision

[`spec/SPEC.md`](../../spec/SPEC.md) and [`spec/EXAMPLES.md`](../../spec/EXAMPLES.md)
are the **sole binding contract**. Every implementation conforms to the spec, not
to any other implementation. MATLAB is not ground truth — if MATLAB and the spec
disagree, MATLAB is wrong. Each port derives its behaviour independently from the
spec; cross-language numerical equivalence is a *consequence* of independent
conformance, never a goal pursued by copying one port into another. The
`testdata/` reference vectors are a check, not a proof.

## Consequences

- **Positive:** joint drift is structurally prevented — a bug must reproduce
  independently in two independent derivations to escape detection, which is far
  less likely than a single copied error. The spec stays the place where
  behaviour is decided, so the next port (Julia) has a complete contract to work
  from.
- **Positive:** changing behaviour has a clear, enforced workflow — update the
  spec first, then every port in the same PR.
- **Negative / cost:** every normative behaviour must be written into the spec
  *before* it is implemented, even when "just porting" would be faster. New
  features need a direct spec-to-implementation read-through, not only a vector
  comparison. Shared helpers (`sidValidateData` / `validate_data`) must be
  audited caller-by-caller against the spec when touched, because a helper bug
  drifts every caller identically.
- A standing follow-on: rules the spec does not yet name are tracked as
  verification gaps (see #113) rather than left implicit.

## Alternatives considered

- **MATLAB as the reference oracle (other ports ported from it).** Rejected:
  makes a bug in MATLAB undetectable by cross-validation, since the vectors are
  generated from MATLAB and every port is matched to them. Removes the
  independent second derivation that is the whole point of having a spec.
- **Reference vectors as the contract (no prose spec).** Rejected: vectors pin
  behaviour only on the specific inputs sampled; they say nothing about edge
  cases, error conditions, defaults, or output-field semantics, and they encode
  no rationale. A new port has no way to recover *why* a value is what it is.
- **No single source of truth (each port documents itself).** Rejected: this is
  the drift scenario by construction — N ports become N contracts that diverge
  the moment any one changes.
