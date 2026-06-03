# Architecture Decision Records

This directory holds sid's Architecture Decision Records. Each ADR documents a
non-obvious tactical or engineering decision with enough context that a future
contributor (human or agent) can understand why the decision was made and when
to revisit it.

## What an ADR is — and isn't

**An ADR records a decision that will stick.** If a future PR might reasonably
want to revisit the decision, or if contributors following or diverging from the
pattern will want a reference, it's ADR-worthy.

**An ADR is NOT:**

- A place for the binding contract. [`spec/SPEC.md`](../../spec/SPEC.md) and
  [`spec/EXAMPLES.md`](../../spec/EXAMPLES.md) own the cross-language contract.
  An ADR can *motivate* a spec change, but the authoritative text lives in the
  spec. When an ADR motivates a spec change, both land in the same PR (spec
  first, per [`CONTRIBUTING.md`](../../CONTRIBUTING.md) §"Workflow for
  algorithmic changes").
- A place for obvious or purely mechanical choices (formatter settings, import
  ordering, internal naming). These churn too much to earn a permanent record.
- A narrative of what you did. It records what you *decided* and why the
  rejected alternatives were rejected.

When to write one vs. not is also covered in
[`CONTRIBUTING.md`](../../CONTRIBUTING.md) §"Design decisions (ADRs)" and
[`CLAUDE.md`](../../CLAUDE.md) §4 (major decisions).

## ADR lifecycle

ADRs can be written in either of two shapes:

- **ADR-first:** the ADR is opened as a PR; discussion happens on the PR; the
  ADR captures the decision when the PR merges.
- **Issue-first:** an issue is opened for the decision (Context + Alternatives +
  where it lands); discussion/resolution happens in the issue; once it closes, a
  follow-up PR writes the ADR from the closed-issue thread.

Issue-first is useful when the decision is part of a batch where opening N empty
ADR-PRs would be heavier than opening N issues. Either shape produces the same
ADR here; the difference is where the discussion lives.

## Format

Use [`ADR-TEMPLATE.md`](ADR-TEMPLATE.md) as the starting point. Required
sections:

1. **Status** — with date. Proposed / Accepted / Rejected / Superseded.
2. **Context** — what forced the decision. Prefer concrete evidence (bugs that
   happened, drift that was observed) over abstract argument.
3. **Decision** — what was chosen. Terse; 1–3 sentences.
4. **Consequences** — what this means going forward. Mix of positive and
   negative.
5. **Alternatives considered** — **not optional.** Each alternative named,
   described, and rejected with enough specificity that a reader can evaluate
   the reasoning.

The "Alternatives considered" section is what makes an ADR useful six months
later. An ADR without it is a decision-log entry, not an ADR.

## Numbering and filenames

- Filenames: `ADR-NNNN-kebab-case-title.md`, `NNNN` zero-padded.
- Numbers are assigned sequentially as ADRs merge, never reused, never
  reordered.
- Titles should be short and specific.
- If an ADR supersedes another, the new one references the old one in its Status
  section, and the old one's Status flips to "Superseded by ADR-NNNN".

## Linking to ADRs

- From PR descriptions: `Implements X per ADR-NNNN.`
- From code comments beside tactical values: `% See ADR-NNNN` (MATLAB) or
  `# See ADR-NNNN` (Python) — useful for magic numbers, thresholds, and fallback
  choices whose rationale lives in an ADR.
- From other docs: markdown links to the ADR file.

## Index

<!-- Append-only, numeric order. ADR number — short title — status — one line. -->

- [ADR-0001](ADR-0001-spec-is-the-contract.md) — Implementations conform to the
  spec, not to each other — **Accepted** — `spec/SPEC.md` + `spec/EXAMPLES.md`
  are the sole contract; each language port derives independently; cross-language
  vectors are a check, not a proof.
