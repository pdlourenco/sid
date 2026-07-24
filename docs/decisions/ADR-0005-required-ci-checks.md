# ADR-0005: Required CI checks are enforced as code

- **Status:** Accepted — _2026-07-24_
- **Deciders:** Project maintainer
- **Related:** [`.github/rulesets/main.json`](../../.github/rulesets/main.json), [`.github/workflows/branch-protection.yml`](../../.github/workflows/branch-protection.yml), [`CONTRIBUTING.md`](../../CONTRIBUTING.md) §CI, #174, #178, #179

## Context

Twice now an admin-merge has landed on `main` while a lint check was red or
incomplete: #174 and #178 both merged before the **MATLAB/Octave Lint** job
finished, and the breakage sat undetected on `main` until #155's local-ci runner
put every leg in one place on its first full outing (#179). The `main` ruleset
required **no** status checks, so nothing blocked the merges.

Two structural facts made this possible:

1. **No check was required.** The ruleset lived only in repo settings and
   required nothing, so a red gate was advisory.
2. **The lint checks couldn't safely be required as they stood.** Every PR-check
   workflow (`lint`, `python-lint`, `python-tests`, `cross-validate`) filtered
   its `pull_request` trigger by path. A path-filtered workflow does not run on a
   PR that misses the paths — and a *required* check that never runs stays
   `pending` forever, wedging the merge. So requiring them naively would have
   broken every docs-only PR.

This is exactly the seed's deferred "branch-protection-as-code + drift workflow"
feature; its epic trigger ("a drift incident happens") has now fired twice.

## Decision

The required-checks policy is enforced **as code**:

- The two **cheap** lint workflows (`lint.yml`, `python-lint.yml`, ~1 min each)
  drop their `pull_request` path filter and run on **every** PR to `main`, with
  distinct job names (`MATLAB/Octave lint`, `Python lint`) so their check
  contexts are unambiguous and safe to require.
- [`.github/rulesets/main.json`](../../.github/rulesets/main.json) declares the
  `main` ruleset — required status checks (`MATLAB/Octave lint`, `Python lint`),
  plus deletion / non-fast-forward protection — with **admin bypass retained**
  (`bypass_mode: always`), so a maintainer can still admin-merge, but as a
  *conscious* override rather than a silent default.
- [`.github/workflows/branch-protection.yml`](../../.github/workflows/branch-protection.yml)
  applies that JSON to the live ruleset on `workflow_dispatch` and drift-checks
  it on PRs that touch the config. It needs a repo-admin token
  (`RULESET_TOKEN`); without it the job no-ops rather than failing.

The expensive checks (`Python Tests` matrix, `Cross-Language Validation` with its
Octave install) are **not** required in this pass — they would need the same
always-run treatment and fork-secret handling first. Extending the required set
to them is a documented follow-up, not a blocker for closing the lint door.

## Consequences

- **Positive:** a red lint gate now blocks a normal merge; the two-admin-merge
  class is closed for non-admins and made a conscious choice for admins.
- **Positive:** the policy is version-controlled and drift-checked — the next
  time the required set changes, it changes in a reviewed PR, not silently in
  settings.
- **Negative / cost:** the two lint checks now run on every PR, including
  docs-only ones (~1 min each). Deemed worth it — these are the cheap checks and
  the ones that were bypassed.
- **Activation cost:** making it live is a one-time maintainer action — add a
  repo-admin `RULESET_TOKEN` secret and run the workflow once (or apply the JSON
  by hand). `main` must be green first (#179 fixes the outstanding lint break).

## Alternatives considered

- **Just flip the lint contexts to required in repo settings.** Rejected as the
  *sole* fix: it's invisible (not in the repo), undoes nothing when settings
  drift, and — without the always-run change — would wedge every PR that doesn't
  touch the filtered paths.
- **Require the checks but keep the path filters.** Rejected: the path-filter /
  required-check wedge is the well-known GitHub trap; a required check excluded
  by a path filter never reports and blocks the merge indefinitely.
- **A single always-running "CI gate" job aggregating the others' conclusions.**
  A valid pattern, rejected for now as heavier than needed: two cheap checks made
  always-run cover the door that actually broke, without a bespoke cross-workflow
  aggregator to maintain.
- **Remove admin bypass entirely.** Rejected: a single-maintainer scientific repo
  needs the escape hatch (e.g. to land a fix when a check is flaky); the goal is
  to make bypass *conscious*, not impossible.
