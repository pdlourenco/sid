# Contributing to sid

Contributions are welcome via issues and pull requests.

## Project Structure

sid is a multi-language system identification toolbox built around a shared
mathematical specification:

- [`spec/`](spec/) — Algorithm specification and mathematical derivations
  (single source of truth for all implementations)
- [`testdata/`](testdata/) — Cross-language reference test vectors (JSON; see [README](testdata/README.md) for format)
- [`matlab/`](matlab/) — MATLAB/Octave implementation (stable)
- [`python/`](python/) — Python implementation (stable)
- [`julia/`](julia/) — Julia implementation (planned)
- [`docs/`](docs/) — **internal** engineering documents: ADRs
  ([`docs/decisions/`](docs/decisions/)), analyses, plans, `DESIGN.md`,
  `REVIEW_CONTEXT.md`, and the function catalogue
  ([`docs/roadmap.md`](docs/roadmap.md)). Not published.
- `docsite/` — source for the **public** documentation site (MkDocs)

### `docs/` and `docsite/` are deliberately separate

`docs/` is internal, `docsite/` is what gets published. The site's `docs_dir` is
`docsite`, so a new file under `docs/` can neither leak onto the public site nor
break its build — which is exactly what happened when the two shared a directory
(the internal docs' relative links failed the strict build). Keep the split:

- Add internal engineering documents to `docs/`. No site changes needed.
- Publishing anything from `docs/` is a **deliberate decision**, made by adding
  it to the `nav` in `mkdocs.yml` — never a side effect of the build config.
- The site regenerates from source (API references, spec includes, executed
  notebooks), so most content changes need no `docsite/` edit at all.

Build it locally before changing anything under `docsite/`, `mkdocs.yml`, or
`scripts/build_matlab_*.py`:

```bash
pip install -r requirements-docs.txt && mkdocs build --strict
```

`--strict` is what CI runs; it turns link warnings into failures. The docs
workflow runs on every pull request, so a broken build blocks before merge.

## Specification as Source of Truth

**[`spec/SPEC.md`](spec/SPEC.md) is the binding contract for every
implementation.** This applies equally to human contributors and AI coding
agents. Read it before touching any algorithmic code.

### Core rules

1. **The spec defines required behaviour.** Defaults, edge cases, error
   conditions, output struct fields, normalization conventions, and
   numerical diagnostics are all part of the contract — not just the core
   formulas. If a requirement is in `SPEC.md`, every implementation must
   satisfy it.

2. **Implementations conform to the spec, not to each other.** MATLAB is
   *not* a ground truth. "The MATLAB version does it this way" is not a
   valid justification — if MATLAB and the spec disagree, MATLAB is wrong.
   Cross-language numerical equivalence is a *consequence* of every
   implementation independently satisfying the spec, not a goal pursued
   by copying one implementation into another.

3. **Fix the spec first when it is ambiguous or wrong.** If you find that
   the spec is silent, contradictory, or disagrees with a clearly correct
   algorithm, **update `SPEC.md` first** (in a dedicated commit, ideally
   with reviewer approval), *then* update the implementations to match.
   Never silently encode behaviour that the spec does not describe — the
   next language port will have no way to recover the same decision.

4. **Shared helpers can cause shared drift.** A bug in a helper like
   `sidValidateData` / `validate_data` can make *every* downstream caller
   silently violate the spec in the same way, which will not be caught by
   cross-validation tests (because the MATLAB and Python CI both exhibit
   the same drift). When you touch a shared helper, audit each caller
   against the relevant spec section.

5. **Cross-language reference vectors are a check, not a proof.** Tests
   in [`testdata/`](testdata/) verify that MATLAB and Python agree
   numerically on specific input data. They do *not* prove either
   implementation satisfies the spec. New features must come with a
   direct spec-to-implementation read-through, not only reference
   vector comparisons.

### Workflow for algorithmic changes

When adding or modifying an algorithm, a function default, an output
field, or any edge-case behaviour:

1. **Read the relevant `SPEC.md` section end to end** before touching
   code. Note every normative statement (defaults, bounds, NaN handling,
   regularization thresholds, output field names and shapes, warnings).
2. **If the spec does not cover what you need**, update `SPEC.md` in a
   dedicated commit before writing code. Include the rationale in the
   commit message.
3. **Implement in every maintained language** (MATLAB and Python today),
   referencing the same spec section. Use the `SPEC.md §X.Y` comment
   convention to mark each step.
4. **Cite the spec in the PR description.** State which sections the
   change implements or modifies, and call out any spec updates.
5. **Write tests against spec requirements**, not against the current
   output of a reference implementation. A test of the form "assert the
   result equals what MATLAB returned today" does not detect joint drift.

### Checklist for reviewers (and self-review)

- [ ] Is every new default/bound/threshold covered by `SPEC.md`?
- [ ] Do all touched functions cite the relevant `SPEC.md §` in comments?
- [ ] If behaviour was added or changed, was `SPEC.md` updated first?
- [ ] Are the MATLAB and Python behaviours derived independently from the
      spec, rather than ported copy-by-copy from one to the other?
- [ ] Are tests written against the spec's requirements, not against the
      current output of the other language?

## Contract artifacts and drift hardening

The [`testdata/`](testdata/) reference vectors are a **contract artifact**: a
generated, committed check that the language ports agree numerically. Their
failure mode is silent — a wrong vector that both ports happen to match keeps
cross-validation green while every port is wrong (core rule 5). The machinery
that produces and consumes them is therefore held to standing rules, hardened
after the #145 cross-validation remediation (see
[ADR-0002](docs/decisions/ADR-0002-contract-artifact-hardening.md)):

1. **Absolute tolerance floors.** Every reference field's tolerance carries an
   absolute floor (`<field>_atol`) alongside the relative one (`<field>_rel`),
   so a near-zero expected value can't make a relative-only comparison vacuously
   pass.
2. **Stored tolerances are authoritative.** Both consumers
   (`testdata/validate_reference.m`, `python/tests/test_cross_validation.py`)
   read the tolerance from the JSON, with one agreed default (`rtol 1e-6`,
   `atol 0`). Tests must not hardcode per-field overrides — tighten or loosen a
   comparison by editing the vector's `tolerance` block.
3. **No orphan artifacts.** Every committed `reference_*.json` is read by a
   validator in every language that has the port. A reference file no test reads
   is worse than none — it reads as coverage while providing zero. The consumer
   lands in the same PR as the vector.
4. **Payloads change only by regeneration.** A vector's *payload* (the `input` /
   `output` numbers) changes only by re-running `testdata/generate_reference.m`;
   never hand-edit a committed payload, and a payload change with no `matlab/**`
   diff to justify it is a red flag. The `tolerance` block *may* be edited
   directly (rule 2), but the generator's matching tolerance entries change in
   the **same PR**, so a fresh regeneration reproduces the committed file
   byte-for-byte — itself a checkable invariant (rule 5). Each JSON also carries a
   `provenance` block (generator + source commit SHA/date) stamped at generation
   and required by a structural gate in both validators, so a stale or
   hand-edited payload is caught (#172).
5. **Structural gates over outcome tests.** Prefer a gate that checks the
   *mechanism* — every stored field is read and compared, every tolerance
   honored — over one that asserts only a specific numeric outcome. Outcome
   tests pass right up until the generator drifts; structural gates catch it.

These make the vectors trustworthy enough to *adjudicate* a numerical change —
the role they play in every algorithmic PR — rather than merely accompany it.
Conformance is still to the spec, not to the vectors: agreement is evidence, not
proof (ADR-0001).

## Design decisions (ADRs)

Non-obvious tactical or engineering choices live in
[`docs/decisions/`](docs/decisions/) as Architecture Decision Records. See
[`docs/decisions/README.md`](docs/decisions/README.md) for the format,
lifecycle (ADR-first / issue-first), and the current index.

**Write an ADR when** you make a choice someone could reasonably challenge later
— a regularization threshold, a fallback ordering, NaN-handling policy, a
numerical-diagnostic cutoff — and future PRs will either follow it or explicitly
deviate. **Don't write one** when the choice is already fixed by
[`spec/SPEC.md`](spec/SPEC.md) (the spec owns the contract; ADRs capture the
*why* of tactical choices beneath it) or is purely mechanical (formatter
settings, import ordering, internal naming).

When an ADR motivates a spec change, both land in the same PR (spec first, per
"Workflow for algorithmic changes" above). Link ADRs from PR descriptions
(`Implements X per ADR-NNNN.`) and from code comments beside tactical values
(`% See ADR-NNNN`). Adding or revisiting a load-bearing decision is a "major
decision" — see [`CLAUDE.md`](CLAUDE.md) §4.

## Pre-push self-review (agent convention)

Before pushing on a PR branch, review the local diff against the project's
principles — not just lint — and act on findings before pushing. This catches
the "I'd have caught that if I'd thought harder" class of bug before it burns CI
minutes and reviewer attention. Note the outcome in the PR description
(`pre-push review: no findings` or `pre-push review flagged X, fixed in <sha>`).

Coding agents should launch a reviewer subagent on the diff, seeded with
[`docs/REVIEW_CONTEXT.md`](docs/REVIEW_CONTEXT.md) (project principles, red
flags, review modes) alongside [`spec/SPEC.md`](spec/SPEC.md) so it reviews
against what the project cares about. The prompt:

> Review the diff below against sid's principles in `docs/REVIEW_CONTEXT.md`,
> citing the principle number or the affected `spec/SPEC.md` rule:
>
> 1. **Contract conformance** — does any change make a port match another port
>    rather than the spec (principles 1, 2)? Is new/changed behaviour reflected
>    in `spec/SPEC.md` *first* (principle 3)? Does a touched shared helper have
>    every caller audited (principle 4)?
> 2. **Verification** — does a new/changed spec rule name a `Verified by:`
>    mechanism — or, until #113 lands, is it pinned by a test (principle 6)? Are
>    tests written against spec requirements, not against the other language's
>    current output (principle 5)? Are NaN/Inf substitutions, clamps, or warning
>    identifiers changed without a spec update (principle 8)?
> 3. **Scope & docs hygiene** — files or refactors outside the PR's stated
>    purpose; user-facing docs (README, API reference, examples) carrying
>    dev-tracking references (`ADR-NNNN` / `#issue`) or "recently changed"
>    narration (principle 9).
> 4. **Decisions deserving an ADR** — new thresholds, fallbacks, NaN policy,
>    magic numbers (see `docs/decisions/`); and missing ADR links in the PR.
> 5. **Catalogue / naming / auto-discovery** — new public functions follow the
>    `sid + Domain + Method` convention and appear in `docs/roadmap.md`; new
>    tests/examples match the discovery naming convention (principle 7).
>
> Report findings in under 200 words. Say "no findings" if the diff is clean.

**Exceptions** — one-line typo fixes, formatting-only changes, and pure reverts
don't warrant the ceremony.

## Known bugs and their lifecycle

sid ships known limitations and deferred fixes **visibly**, not silently. When a
change leaves a known bug, a rough edge, or a deliberately-deferred follow-up in
place, it carries a **visible-debt marker** at the site, pointing at an open
tracking issue — the style already used across both ports:

```matlab
% ... an unguarded search locks onto spurious tail gaps (AR(1) -> n=40).
% Data-aware-floor follow-up: issue #160.
```

The rules:

1. **Every known-limitation marker names an open issue.** The marker is a
   pointer; the issue is the register entry (context, why it's deferred, what
   would resolve it). A marker with no tracking issue is invisible debt — exactly
   what this convention exists to prevent. Use `issue #NNN` (optionally with the
   governing `SPEC §X.Y`) in a comment at the site, in **every port** where the
   limitation exists.
2. **A fix closes the loop in one PR.** The PR that fixes a known bug **closes
   its tracking issue and removes the visible-debt marker(s) in the same PR** —
   marker and open issue disappear together, so the code and the tracker never
   disagree about what's still broken. This is the known-bug row on the
   [pull request template](.github/pull_request_template.md).
3. **Partial fixes re-point, they don't orphan.** A PR that resolves part of a
   tracked limitation and defers the rest re-points the marker at the follow-up
   issue rather than deleting a marker while the debt remains.

See [ADR-0003](docs/decisions/ADR-0003-known-bug-lifecycle.md) for the rationale.

## General Guidelines

- The spec rules above apply to every implementation. The per-language
  guides below cover language-specific style, naming, and testing.
- Cross-language test vectors in [`testdata/`](testdata/) ensure numerical
  consistency across implementations. New algorithms should include
  reference vectors.
- The project is MIT-licensed. See [`LICENSE`](LICENSE).

## Code Style

The root [`.editorconfig`](.editorconfig) enforces basic formatting rules
(UTF-8, LF line endings, trailing whitespace). Language-specific linting
and style rules are documented in each language's contributing guide.

## Language-Specific Guidelines

Each language has its own contributing guide with conventions for naming,
documentation, code style, and testing:

- **MATLAB/Octave**: [`matlab/CONTRIBUTING.md`](matlab/CONTRIBUTING.md)
- **Python**: [`python/CONTRIBUTING.md`](python/CONTRIBUTING.md)

## Test and Example Auto-Discovery

Test and example runners in every language **discover files by naming
convention** — there is no hardcoded manifest to maintain. To add a test
or example, create a file matching the pattern for that language:

| Language | Tests | Examples |
|----------|-------|----------|
| MATLAB/Octave | `test_*.m` | `example*.m` |
| Python | `test_*.py` | `example_*.py` |
| Julia | `test_*.jl` | `example_*.jl` |

Runners sort discovered files alphabetically and execute them in order.
**Do not maintain hardcoded file lists** — auto-discovery prevents the
common failure mode where a new test exists but is never executed because
it was not added to a manifest.

### Templates

Each language provides template files for tests and examples. Copy the
template when creating a new file — it includes the runner instrumentation
variables that enable per-file progress tracking in CI output.

| Language | Test template | Example template |
|----------|--------------|-----------------|
| MATLAB/Octave | `matlab/tests/test_template.m` | `matlab/examples/example_template.m` |

Each language's contributing guide documents the discovery mechanism and
templates in detail. When starting a new language port, implement the runner
with auto-discovery from day one and provide starter templates.

## CI

CI workflows run per-language:

- **MATLAB Tests** — MATLAB and GNU Octave test suites
- **MATLAB Lint** — MISS_HIT style/lint checks and function header validation
- **Python Lint** — ruff style/lint checks and docstring validation
- **Python Tests** — pytest on Python 3.10–3.13
- **Cross-Language Validation** — reference test vector consistency

All checks must pass before merging. The **MATLAB/Octave lint** and **Python
lint** checks run on every PR and are **required** by the `main` ruleset — a red
lint gate blocks a normal merge. That policy is version-controlled in
[`.github/rulesets/main.json`](.github/rulesets/main.json) and applied /
drift-checked by [`.github/workflows/branch-protection.yml`](.github/workflows/branch-protection.yml)
(see [ADR-0005](docs/decisions/ADR-0005-required-ci-checks.md)); admins can still
bypass, but as a conscious override.

**Run them locally first.** [`scripts/local-ci`](scripts/local-ci) mirrors these
legs so the whole suite is one command before you push, instead of tribal
knowledge:

```
scripts/local-ci            # every leg whose toolchain is installed
scripts/local-ci python     # just the Python legs
scripts/local-ci lint       # just the lint legs
```

A leg whose toolchain isn't installed is **skipped** (with a hint), not failed,
so it's useful in a partial dev environment; the run exits non-zero only if a leg
that actually ran fails. Keep the script in sync with `.github/workflows/*.yml`.
