# ADR-0002: Contract artifacts are drift-hardened, not just present

- **Status:** Accepted — _2026-07-24_
- **Deciders:** Project maintainer
- **Related:** [`CONTRIBUTING.md`](../../CONTRIBUTING.md) §"Contract artifacts and drift hardening", [ADR-0001](ADR-0001-spec-is-the-contract.md), #145, #147

## Context

The `testdata/` reference vectors are the mechanism that makes cross-language
equivalence checkable — the check that gives ADR-0001 its teeth (agreement
between two *independent* derivations is evidence). But a contract artifact's
failure mode is silent: the #145 repo-wide review found the cross-validation
machinery had blind spots that let it pass green while checking almost nothing.

- The reference generator double-scaled `S`, so the stored uncertainty `P`
  matched neither implementation — and the test never read `P` at all, so the
  whole uncertainty path had zero effective cross-language coverage.
- Four reference JSONs (`detrend`, `freq_domain_sim`, `ltv_io`, `uncertainty`)
  were committed but wired into no Python consumer — coverage on paper, none in
  fact.
- Per-field tolerances were hardcoded per language, so the two consumers could
  silently disagree on how tight a comparison actually was.

A contract artifact that exists but isn't rigorously produced and consumed is
worse than none: it reads as coverage while providing it. #145/#147 fixed the
specific instances; this ADR records the standing doctrine so the next artifact
doesn't reintroduce the class. (The seed's four-adopter study independently
converged on the same doctrine — evidence it's the right generalization, not a
local overfit.)

## Decision

Contract artifacts — today the `testdata/` vectors — are held to five standing
rules, stated in [`CONTRIBUTING.md`](../../CONTRIBUTING.md) §"Contract artifacts
and drift hardening": absolute tolerance floors; stored tolerances authoritative
(no per-language hardcoding); no orphan artifacts (every vector wired into every
port's validator, in the same PR); regeneration is the only sanctioned edit; and
structural gates over outcome tests. Recording generator provenance *inside* each
artifact is adopted as a **target**, tracked separately — not yet a gate.

## Consequences

- **Positive:** the vectors become trustworthy enough to *adjudicate* a
  numerical change, not merely accompany it. The silent-pass class #145 found is
  structurally excluded, not just patched at three sites.
- **Positive:** one place — CONTRIBUTING — states the rules a reviewer applies to
  any PR touching `testdata/`, and the pre-push reviewer prompt can cite them.
- **Negative / cost:** adding a reference vector is now strictly more work — its
  consumer, in every port, lands in the same PR, and a near-zero field needs an
  `atol` chosen deliberately rather than defaulted.
- **Cost / visible debt:** the generator-provenance target is stated but unmet,
  so it remains debt until a follow-up lands the metadata and a check for it.

## Alternatives considered

- **Leave it at the #145 fixes (no standing doctrine).** Rejected: the three
  fixes address instances, not the class; the next orphan vector or relative-only
  tolerance would have nothing to catch it. Standing policy is what turns a
  one-off remediation into a gate.
- **Make generator-provenance metadata a hard rule now.** Rejected: the metadata
  and its CI check don't exist yet; asserting a rule the repo currently fails
  would make the doctrine a lie. Stated as a tracked target instead — honest
  about the current state, which is itself one of the rules (no orphan claims).
- **Put the rules only in `testdata/README.md`.** Rejected: reviewers work from
  CONTRIBUTING and the reviewer prompt, so the binding rules belong where the
  review discipline already points; `testdata/README.md` stays the how-to.
