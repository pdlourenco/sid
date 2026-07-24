# Analyses

This directory holds **analyses**: point-in-time investigations of the codebase
— reviews, audits, numerical cross-checks, root-cause write-ups. An analysis
records what was found *at a specific commit*, and the findings feed the binding
artifacts (issues, ADRs, spec changes); the analysis itself is a dated record,
not a living document.

The convention (four rules):

1. **Dated.** The filename is `YYYY-MM-DD-slug.md` (e.g.
   `2026-07-12-repo-wide-review.md`). The date is when the analysis was
   performed, and files sort chronologically.
2. **Commit-anchored.** The header states the commit the analysis was performed
   at, and any `file:line` references are relative to that commit — so a reader
   can check a finding against the exact tree it was found in, even after the
   code moves on.
3. **Immutable.** An analysis is a record of a moment, not a page to keep
   current. Don't rewrite it as the code changes: fix a typo, but a *finding*
   that was later resolved, refined, or refuted is followed up in the issue it
   spawned or in a new dated analysis — never by editing the original into
   disagreement with what it actually said. This is what makes it citable
   evidence (and pairs with review principle 9 / ADR-0001's "vectors cannot
   catch joint drift" being provable *because* the finding is fixed in time).
4. **Not a contract.** Analyses *motivate* decisions; they don't bind. The
   authoritative artifacts are [`spec/SPEC.md`](../../spec/SPEC.md) (the
   contract), [`docs/decisions/`](../decisions/) (ADRs), and the issue tracker
   (known-bug register). An analysis links *out* to those; nothing conforms *to*
   an analysis.

Convention adopted from the disciplined-project-seed (ADR-0010); recorded here
rather than as a sid ADR because the README *is* the convention statement.
`docs/analyses/2026-07-12-repo-wide-review.md` already complies (dated,
`Reviewed at: <commit>` anchor, immutable record).
