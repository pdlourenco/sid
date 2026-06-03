# CLAUDE.md — Agent operating rules

Binding rules for Claude (and any other coding agent) working on the **sid**
repository. Source-of-truth documents are linked inline; this file is a short
index, not a replacement for them.

sid is a multi-language system identification toolbox (MATLAB/Octave, Python,
Julia planned) built around a single mathematical specification. Most
development happens through coding agents editing independent language ports in
parallel against that shared contract — these rules are the guardrails that keep
the ports from silently drifting apart.

## 1. Follow `CONTRIBUTING.md`

[`CONTRIBUTING.md`](CONTRIBUTING.md) is authoritative for the spec-as-source-of-truth
discipline, the workflow for algorithmic changes, the reviewer checklist, CI, and
test/example auto-discovery. The per-language guides cover language-specific
style, naming, and testing:

- MATLAB/Octave: [`matlab/CONTRIBUTING.md`](matlab/CONTRIBUTING.md)
- Python: [`python/CONTRIBUTING.md`](python/CONTRIBUTING.md)

Read the relevant guides before opening a PR. Deviations require an explicit note
in the PR description.

## 2. Pre-push self-review is expected

Before pushing on a PR branch, review the local diff against the project's
principles — not just surface-level lint — and act on what you find before
pushing. This catches the "I'd have caught that if I'd thought harder" class of
bug before it burns CI minutes and reviewer attention.

Until the dedicated reviewer-context document lands (tracked in #115), use the
**"Checklist for reviewers (and self-review)"** in [`CONTRIBUTING.md`](CONTRIBUTING.md)
as the review basis, paying particular attention to:

- Is every new default / bound / threshold covered by `spec/SPEC.md`?
- Do the MATLAB and Python sides derive their behaviour **independently from the
  spec**, rather than one being ported copy-by-copy from the other?
- Are tests written against the spec's requirements, not against the current
  output of the other language?

The narrow exceptions — one-line typo, formatting-only, pure revert — don't
warrant the ceremony. Note the outcome briefly in the PR description.

## 3. Implementation is bound to the spec as a contract

[`spec/SPEC.md`](spec/SPEC.md) (algorithm specification) and
[`spec/EXAMPLES.md`](spec/EXAMPLES.md) (example-suite contract: plant catalog,
helper API, per-example structure every port conforms to) define the binding
contracts every implementation must satisfy. Implementations conform to the
spec, **not to each other**. Concretely:

- **The spec defines required behaviour** — defaults, edge cases, error
  conditions, output fields, normalization conventions, and numerical
  diagnostics are all part of the contract, not just the core formulas.
- **MATLAB is not ground truth.** "The MATLAB version does it this way" is not a
  justification. If MATLAB and the spec disagree, MATLAB is wrong. Cross-language
  numerical equivalence is a *consequence* of each port independently satisfying
  the spec, never a goal pursued by copying one port into another.
- **If you need to change a contract, update the spec first** (ideally a
  dedicated commit with the rationale), then update every implementation side in
  the same PR. Follow the "Workflow for algorithmic changes" in
  [`CONTRIBUTING.md`](CONTRIBUTING.md). Never silently encode behaviour the spec
  does not describe — the next language port has no way to recover the decision.
- **Shared helpers cause shared drift.** A bug in a helper like `sidValidateData`
  / `validate_data` makes *every* caller violate the spec the same way, which
  cross-validation will **not** catch (both languages exhibit the same drift).
  Audit each caller against the relevant spec section when you touch a shared
  helper.
- **Reference vectors are a check, not a proof.** The [`testdata/`](testdata/)
  vectors verify MATLAB and Python agree numerically on specific inputs; they do
  not prove either satisfies the spec. New features need a direct
  spec-to-implementation read-through, not only vector comparison.
- The cross-validation gate in CI is not optional. Do not weaken it to make a PR
  green. If the spec and an implementation drift and you cannot tell which is
  correct, **stop and ask** — do not pick a side unilaterally.

## 4. Discuss major decisions before deciding; record one if it sticks

A "major decision" is anything that:

- changes a contract in `spec/SPEC.md` or `spec/EXAMPLES.md`, or the public
  function catalogue / naming convention in [`docs/roadmap.md`](docs/roadmap.md);
- changes a load-bearing governance document the whole workflow leans on — this
  file, [`CONTRIBUTING.md`](CONTRIBUTING.md), or the reviewer-context document
  (`docs/REVIEW_CONTEXT.md`, landing in #115);
- introduces a new external dependency, a new public function, or a new on-disk
  artifact (e.g. a new reference-vector format);
- locks in a trade-off a future PR could reasonably want to revisit
  (regularization thresholds, fallback ordering, NaN handling, error-handling
  policy);
- materially changes the scope or shape of the work in flight.

For any of the above:

1. **Pause and surface the decision — recommend, don't decide.** Lay out the
   choice, the alternatives, and the trade-off in the conversation, and mark the
   option you'd choose with a one-line *why*; then let the maintainer choose.
   Wait for an explicit go-ahead before implementing. This
   recommend-don't-decide posture is the default for *every* maintainer-facing
   question, not only the major decisions enumerated above: surface the options
   and your recommendation rather than picking silently.
2. **If the decision is accepted and non-obvious, record it as an ADR** in
   [`docs/decisions/`](docs/decisions/) following
   [`docs/decisions/README.md`](docs/decisions/README.md), and link it from the
   PR description. When an ADR motivates a spec change, both land in the same PR
   (spec first).
3. **Tactical and mechanical choices do not need this** — formatter settings,
   import ordering, internal naming, obvious refactors. When in doubt, ask; the
   cost of a question is lower than the cost of an unwanted commit.

## 5. Opening PRs, and the two-session authoring/review split

You may **open** a PR for work that is already planned or pre-approved — a filed
issue, an epic sub-task, an agreed refactor — without asking first. Opening a PR
is how work is proposed for review, not how it is committed to `main`.
**Merging always requires explicit maintainer approval:** never merge your own
PR, and never weaken or skip a required check to get one green (§3).

Authoring and review belong in **different sessions**. The session that writes a
change runs the pre-push self-review (§2) on its own diff, but a reviewer that
shares the author's context also shares its blind spots — so a *separate*
reviewer session (or a subagent with fresh context) reviews the pushed PR
against the same principles §2 lists, seeded by `docs/REVIEW_CONTEXT.md` once it
lands (#115). Keep the two roles apart: don't approve, in substance, a diff your
own session authored.

Session-level direction from the maintainer overrides this file. If, in a
session, the maintainer tells you to do something this document defers or
forbids, the in-session instruction wins for that task — this file is the
default, not a veto over the person you're working with.
