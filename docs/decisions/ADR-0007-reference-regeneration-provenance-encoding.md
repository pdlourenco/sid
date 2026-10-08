# ADR-0007: Reference vectors are regenerated on the PR, provenance-checked, and encode non-finite values

- **Status:** Accepted — _2026-10-08_. Supersedes [ADR-0002](ADR-0002-contract-artifact-hardening.md).
- **Deciders:** Project maintainer
- **Related:** [`docs/analyses/2026-10-04-repo-wide-review-v2.md`](../analyses/2026-10-04-repo-wide-review-v2.md) §4.1 (F1, F3, F5), [`docs/plans/2026-10-04-review-v2-remediation-plan.md`](../plans/2026-10-04-review-v2-remediation-plan.md) items A2 and A3 (decided 2026-10-08, §0.3), [`testdata/README.md`](../../testdata/README.md), [ADR-0005](ADR-0005-required-ci-checks.md)

## Context

[ADR-0002](ADR-0002-contract-artifact-hardening.md) made the `testdata/` vectors
contract artifacts held to five standing rules (absolute tolerance floors; stored
tolerances authoritative; no orphan artifacts; payloads change only by
regeneration; structural gates over outcome tests), with CI-pinned MATLAB R2025a as
the canonical generator and a `provenance` block in every vector. The 2026-10-04
review found three ways the machinery around those rules failed:

- **Regeneration lands unvalidated and churns (F1).** CI regenerates the vectors on
  every push to `main` and commits the result under `[skip ci]`, so no validator
  runs on the new bytes. The pinned runner is not byte-reproducible: one
  regeneration with no code change rewrote 21 of 32 payloads, `ltv_io.A` by
  4.5e-5. Three of the five regeneration commits since provenance landed were
  provenance-only re-stamps after docs-only merges. The mechanism produces the
  churn the rules forbid.
- **Provenance is weaker than stated (F5).** Both consumers accept
  `git_sha: "unknown"`, neither checks that the SHA belongs to the history being
  validated, and the PR-time regeneration is never compared with the committed
  vectors.
- **Non-finite values cannot be stored or compared (F3).** `jsonencode` writes `NaN`,
  `Inf` and `-Inf` all as `null`, so the §3.3 `Inf` sentinel is indistinguishable
  from `NaN` on disk. The Octave validator compares with `any(absDiff > thresh)`,
  which is false for any `NaN`, so a `NaN` on either side passes. No committed
  vector holds a non-finite value today, which is exactly why the degenerate-input
  regime — where the ports are known to diverge — has no cross-vector.

## Decision

sid keeps ADR-0002's five rules unchanged, and adds three:

1. **Regeneration lands on the PR, never on `main`.** No workflow commits vectors to
   `main`. The PR's MATLAB job (canonical R2025a) regenerates and diffs numerically
   against the committed vectors. It pushes the regenerated vectors to a same-repo
   PR branch when the diff reports a semantic change — a payload outside tolerance,
   or any change to a tolerance entry, key set or encoding — or when a changed
   vector's provenance names a non-canonical engine. That push's commit message
   never contains the skip-ci marker, and the job never pushes when the head commit
   is its own. For fork PRs the job uploads the vectors as an artifact and the
   maintainer commits them. On push to `main` the job still regenerates and diffs,
   without committing, and fails on a semantic difference, so two individually green
   PRs that merge into stale vectors are caught.
2. **Provenance is checked.** The generator stamps the full commit SHA and the
   engine version. A PR-time check requires each changed vector's SHA to be an
   ancestor of the PR head and its engine to be the canonical one. Both consumers
   reject a SHA that is `"unknown"` or not hexadecimal, and accept 7–40 hex digits.
3. **Non-finite values are stored as the JSON strings `"NaN"`, `"Inf"` and
   `"-Inf"`.** The generator writes them so; every consumer decodes them back to
   numbers before comparing; a `null` in a vector is an error. A comparison passes
   at a position only if both values are `NaN`, or both are infinite with the same
   sign, or both are finite and within tolerance.

## Consequences

- **Positive:** every vector that reaches `main` has been validated by every
  consumer on the PR that changed it, and engine-ULP churn no longer produces
  commits: a regeneration within tolerance is discarded.
- **Positive:** degenerate-input cross-vectors become possible (`NaN` and the §3.3
  `Inf` sentinel are distinct on disk and compared exactly), which removes the
  validator limitation that the `none` cross-vector tags in `SPEC.md` rested on.
- **Positive:** the strings `"NaN"`, `"Inf"` and `"-Inf"` are valid JSON, so any
  standard parser reads the vectors; only the sid consumers need to know the
  convention.
- **Negative / cost:** the three strings are reserved — a vector field can no longer
  hold one of them as text. Today no field does.
- **Negative / cost:** a PR that changes numerics gets a bot commit on its branch;
  the author must pull before pushing again and must not force-push over it.
- **Negative / cost:** the reference bot no longer needs to push to `main`, so its
  always-bypass in the `main` ruleset goes — but only once rule 1 is in place, since
  until then the bot still pushes regenerated vectors to `main`. That ruleset is
  ADR-0005's subject, and the README offers no way to supersede part of an ADR, so
  this ADR does not change ADR-0005's status: the ADR that supersedes ADR-0005 (the
  review-v2 plan's item A1) records the ruleset without the bypass.
- **Not closed by this ADR:** rule 3 is implemented with this ADR (plan item A3).
  Rules 1 and 2 are implemented by plan item A2; until that lands, CI still
  regenerates and commits on `main`, as `testdata/README.md` and `tests.yml`
  describe.

## Alternatives considered

- **Keep committing regenerations on `main`, but validate them** — a bot PR, or a
  `workflow_run` that validates the bot commit and reverts it on failure. Rejected
  because both keep the churn: every engine-ULP difference still becomes a commit,
  and every docs-only merge still re-stamps provenance.
- **Commit on `main` only when the numeric diff exceeds tolerance.** Rejected because
  the regeneration still happens after review, so the vector that lands has never
  been seen by the PR's validators or its reviewer.
- **Accept the churn and drop the "no churn" rule.** Rejected because a rule the
  machinery breaks on every push is not a rule, and churn hides real changes in
  noise.
- **Check that the provenance SHA is a real commit, or allow only merge commits.**
  A "real commit" check would accept commits from abandoned PRs and needs
  `refs/pull/*`, which a local clone lacks; allowing only merge commits is workable
  (every merge so far has been one) but constrains the repository to fit a check.
  An ancestor-of-the-PR-head check at PR time works under every merge method.
- **A parallel `<field>_mask` for non-finite positions.** Rejected because it doubles
  the fields a reviewer reads, and a mask and a payload can disagree.
- **Bare `NaN` / `Infinity` tokens** (`jsonencode` with `ConvertInfAndNaN` off).
  Rejected because they are not JSON (RFC 8259): JavaScript's `JSON.parse` and
  other strict parsers reject the file, and lenient ones read it lossily (`jq` turns
  `NaN` into `null` and `Infinity` into the largest double), so the vectors would
  stop being reliably readable by standard tools. sid's MATLAB and Octave consumers do use these tokens internally, as an
  intermediate step when decoding.
- **Keep `null` and treat it as `NaN`.** Rejected because it leaves the §3.3 `Inf`
  sentinel unrepresentable.
