# SID — disciplined-project-seed adoption plan

**Status:** Proposed — 2026-06-02. Awaiting maintainer go-ahead per item.

This plan sequences the adoption of the higher-value pieces of the
[disciplined-project-seed](https://github.com/pdlourenco/disciplined-project-seed)
into SID. It is deliberately scoped: SID already realizes the seed's *core* — a
single binding specification (`spec/SPEC.md`) with cross-language reference
vectors (`testdata/`) — and is more mature there than the empty seed template.
The gap is in the layers the seed grew *after* SID forked the concept: agent
governance, a decision trail, verification annotations, and infra-as-code.

None of the work below touches algorithmic code. Every change is process, docs,
or CI scaffolding — low blast radius, easy to revert, no spec-behaviour impact.

## Guiding principles for this adoption

1. **Adapt to SID's layout, don't impose the seed's.** SID keeps contracts in
   `spec/` (not `docs/SPEC.md`) and `CONTRIBUTING.md` at the repo root. The seed
   is "opinionated about shape, agnostic about contents" — we keep SID's paths
   and fold seed conventions into them.
2. **Highest leverage first.** Order by value-to-effort for an agent-driven,
   spec-bound, largely single-maintainer scientific toolbox.
3. **Each item is its own PR.** Small, reviewable, independently revertible. One
   topic per PR keeps CI failures diagnosable.
4. **Dogfood the deferred-with-conditions discipline.** Lower-value pieces are
   not dropped — they are parked in §Deferred with explicit triggers.

---

## Phase 1 — Agent governance + verification (highest leverage)

The single biggest gap. SID is developed almost entirely by coding agents (the
git history is overwhelmingly `claude/*` branches), yet agents get no binding
operating rules and the spec's rules name no verifier. This phase closes both.

### PR 1.1 — `CLAUDE.md` (the keystone)

**Add** a root `CLAUDE.md` adapted from the seed, re-pointed at SID's layout.

- Section 1: follow `CONTRIBUTING.md` (authoritative for CI, spec discipline,
  commit/branch conventions).
- Section 2: pre-push self-review is mandatory (forward-reference PR 2.2; until
  that lands, point at the existing reviewer checklist in `CONTRIBUTING.md`).
- Section 3: implementation is bound to `spec/SPEC.md` **and** `spec/EXAMPLES.md`
  as contracts — change the spec first, then every implementation in the same
  PR; never silently encode behaviour the spec doesn't describe.
- Section 4: discuss major decisions before deciding; write an ADR if it sticks
  (forward-reference PR 2.1).

**Files:** `CLAUDE.md` (new).
**Effort:** ~1–2 hours. Mostly re-pointing `docs/` → `spec/` and folding in
SID's existing "MATLAB is not ground truth" / "shared helpers cause shared
drift" rules.
**Success criteria:** every binding rule references a real SID file; no dangling
links to docs that don't exist yet (use forward-references with a note, or stub).
**Dependencies:** none. Do this first.

### PR 1.2 — `.github/pull_request_template.md`

**Add** a PR template adapted from the seed: spec sections touched, pre-push
review outcome row, local-CI outcome row. Drop the seed's label/ADR rows until
those conventions land (PR 2.1, Phase 3).

**Files:** `.github/pull_request_template.md` (new).
**Effort:** ~30 min.
**Success criteria:** opening a PR pre-fills the template; rows map to SID's
actual workflow (spec citation is already informal practice).
**Dependencies:** none.

### PR 1.3 — `Verified by:` annotations + a Verification section in `spec/SPEC.md`

**Add** the seed's right-side verification convention to SID's *already-real*
spec: a short "Verification (right-side mechanisms)" section near the top, and a
`Verified by:` line on each normative rule naming the mechanism that gates it —
a specific test, a `testdata/` reference vector, a CI job, or manual inspection.

**Files:** `spec/SPEC.md` (and optionally the `spec/cosmic/*.md` sub-specs).
**Effort:** Medium — the framing is easy; the honest accounting is the work.
SID's spec is large (v1.0.0, many sections), so retrofitting will surface real
test gaps. Recommend doing it section-by-section across a couple of PRs rather
than one monster diff.
**Success criteria:** every normative rule either names a verifier or is visibly
marked as uncovered debt; the uncovered set is enumerated so it can be burned
down.
**Dependencies:** none, but highest payoff *after* CLAUDE.md makes the spec the
agent-governed contract.

**Why this phase first:** `CLAUDE.md` + `Verified by:` together turn SID's
strongest asset (a real, enforced spec) into an agent-governed, auditable
contract — which is the seed's entire thesis. Lowest effort, highest leverage.

---

## Phase 2 — Decision trail + reviewer discipline

Codifies *why* decisions were made and makes the existing reviewer checklist a
reusable, agent-runnable convention.

### PR 2.1 — ADR infrastructure (`docs/decisions/`)

**Add** `docs/decisions/README.md` (conventions + ADR-first/issue-first
lifecycle + index), `docs/decisions/ADR-TEMPLATE.md`, and a short
"Design decisions (ADRs)" section in `CONTRIBUTING.md`. Optionally seed ADR-0001
retroactively for one already-made SID decision (e.g. "MATLAB is not ground
truth; implementations conform to the spec independently") to make the trail
concrete from day one.

**Files:** `docs/decisions/README.md`, `docs/decisions/ADR-TEMPLATE.md` (new);
`CONTRIBUTING.md` (edit); optional `docs/decisions/ADR-0001-*.md`.
**Effort:** Low to install; ongoing cost is cultural.
**Success criteria:** the next non-obvious tactical decision (a threshold, a
fallback order, NaN handling, λ-tuning policy) lands with an ADR linked from its
PR.
**Dependencies:** referenced by `CLAUDE.md §4` — close that forward-reference
once this lands.

### PR 2.2 — `docs/REVIEW_CONTEXT.md` + pre-push self-review convention

**Add** `docs/REVIEW_CONTEXT.md` capturing SID's actual principles as **numbered
assertions** a reviewer agent can cite by number — e.g. "(1) the spec is the
contract; implementations conform to it, not to each other; (2) shared helpers
cause shared drift — audit every caller when touching `validate_data`;
(3) reference vectors are a check, not a proof; (4) every normative behaviour
must be in `spec/SPEC.md` first." Then add a "Pre-push self-review (agent
convention)" section to `CONTRIBUTING.md` and wire `CLAUDE.md §2` to it.

**Files:** `docs/REVIEW_CONTEXT.md` (new); `CONTRIBUTING.md`, `CLAUDE.md` (edit).
**Effort:** Medium — value depends on writing SID-specific principles well, not
pasting the template. SID already has the raw material in `CONTRIBUTING.md`'s
"Core rules" and reviewer checklist; this promotes them into a reviewer prompt.
**Success criteria:** a reviewer subagent run against a diff produces findings
that cite SID principles by number; the convention is noted in PR descriptions.
**Dependencies:** strongest after PR 1.3 (so verification-mode findings can cite
`Verified by:` mechanisms).

### PR 2.3 — `docs/DESIGN.md` (rationale half of the DESIGN/SPEC pair)

**Add** a `docs/DESIGN.md` that gives SID's architectural rationale a home,
separating it from the README's "How It Works" and from the spec's contract
text: why Blackman-Tukey, why the block-tridiagonal O(N) COSMIC formulation, why
pooled least-squares for multi-trajectory, why polyglot-by-shared-spec. Link to
`spec/SPEC.md` for contract details rather than duplicating field tables.

**Files:** `docs/DESIGN.md` (new); cross-links from `README.md`, `spec/SPEC.md`.
**Effort:** Medium — this is writing, not scaffolding.
**Success criteria:** the README's "How It Works" can shrink to a pointer; the
"why" questions a new contributor asks have a single home.
**Dependencies:** none, but most useful once REVIEW_CONTEXT references it as a
context doc.

---

## Deferred (lower value at SID's current scale — adopt on trigger)

These are real seed features, deliberately parked with explicit triggers rather
than dropped. Same discipline the seed applies to its own deferred items.

### Label-as-code + issue templates

- **What:** `docs/LABELS.md` + `.github/labels.yml` + sync-labels workflow;
  `.github/ISSUE_TEMPLATE/{bug,decision-proposal}.yml`.
- **Why deferred:** the lifecycle taxonomy (`discussion`/`decided`/`ready`/
  `deferred`) is overhead without a real issue backlog; SID is near-single-
  maintainer today.
- **Trigger to wire in:** issue volume grows past ~one screen of open issues, or
  a second regular contributor joins. Pairs naturally with the reviewer's
  label/lifecycle-currency checks from the seed.

### Branch-protection-as-code + drift workflow

- **What:** `.github/branch-protection.yml` + `scripts/setup-branch-protection.sh`
  + weekly `check-branch-protection.yml`.
- **Why deferred:** SID's required checks are few and stable; the
  YAML+script+jq+drift-workflow machinery (and the admin `gh auth` / `issues:
  write` it needs) is heavyweight for that. Branch protection can stay
  configured in the GitHub UI for now.
- **Trigger to wire in:** required-check set starts changing often, or a
  drift/misconfiguration incident actually happens.

### "Reviewing an open PR" subscribe-and-follow-through convention

- **What:** the parameterised reviewer-invocation section (mode / output /
  CI-handling / subscription) from the seed's `CONTRIBUTING.md`.
- **Why deferred:** it's the most elaborate convention in the seed and only
  earns its keep once pre-push review (PR 2.2) and REVIEW_CONTEXT exist.
- **Trigger to wire in:** immediately after Phase 2 lands and the team wants
  agent review on *remote* PRs, not just local diffs.

### `docs/plans/` phase plans

- **What:** the ROADMAP + `plans/` split (terse roadmap entry → detailed
  per-phase plan).
- **Why deferred:** SID's roadmap is a function-catalogue table, not a phase
  sequence — the phase-plan shape doesn't fit the current planning style.
- **Trigger to wire in:** a multi-PR phased effort (e.g. the Julia port, or the
  deferred online/recursive COSMIC `§8.10`) is scheduled and needs a PR-sequence
  execution view. This very document is a one-off plan; the trigger is *recurring*
  planning need.

### Project-level `CHANGELOG.md`

- **What:** a root Keep-a-Changelog consolidating the per-language release notes.
- **Why deferred:** SID already has per-language `RELEASE_NOTES.md`; a single
  changelog duplicates them.
- **Trigger to wire in:** the per-language notes start diverging or a unified
  release cadence emerges.

---

## Suggested sequencing summary

| Order | PR | Value | Effort | Blocks |
|------|------|-------|--------|--------|
| 1 | 1.1 `CLAUDE.md` | ★★★ | Low | — |
| 2 | 1.2 PR template | ★★ | Trivial | — |
| 3 | 1.3 `Verified by:` in spec | ★★★ | Med | — |
| 4 | 2.1 ADR infra | ★★ | Low | closes CLAUDE §4 ref |
| 5 | 2.2 REVIEW_CONTEXT + pre-push review | ★★ | Med | closes CLAUDE §2 ref |
| 6 | 2.3 `DESIGN.md` | ★★ | Med | — |

Phase 1 is the high-leverage core and can land in days. Phase 2 follows as the
decision/review discipline takes hold. Everything in §Deferred waits for its
named trigger.
