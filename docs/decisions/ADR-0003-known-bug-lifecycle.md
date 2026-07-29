# ADR-0003: Known bugs are visible debt with a tracked lifecycle

- **Status:** Accepted — _2026-07-24_
- **Deciders:** Project maintainer
- **Related:** [`CONTRIBUTING.md`](../../CONTRIBUTING.md) §"Known bugs and their lifecycle", [`.github/pull_request_template.md`](../../.github/pull_request_template.md), [ADR-0001](ADR-0001-spec-is-the-contract.md)

## Context

sid is developed largely by coding agents in parallel, and not every rough edge
can be fixed the moment it's found — a fix may need a spec decision, or belong to
a later remediation phase. The failure mode to avoid is *invisible* debt: a
workaround or known-wrong path that lives in the code with nothing pointing at
it, so the next contributor (or agent) can't tell a deliberate deferral from an
oversight, and the tracker and the code disagree about what is still broken.

sid already does this de-facto — the #134–#160 remediation left inline
`issue #NNN` markers at deferred sites in both ports (e.g. `sidModelOrder.m` and
`model_order.py` both carry `Data-aware-floor follow-up: issue #160`). This ADR
records the convention so it is applied consistently rather than ad hoc, and
pairs it with the known-bug row already on the PR template.

## Decision

A known limitation or deferred fix is **visible debt**: it carries an inline
marker at the site naming an open tracking issue (`issue #NNN`, optionally with
`SPEC §X.Y`), in every port where it exists. The tracking issue is the register
entry. The lifecycle closes in one PR — the PR that fixes the bug closes its
issue and removes the marker(s) together (the known-bug row on the PR template).
A partial fix re-points the marker at the follow-up issue rather than orphaning
it. Stated in [`CONTRIBUTING.md`](../../CONTRIBUTING.md) §"Known bugs and their
lifecycle".

## Consequences

- **Positive:** the code and the issue tracker cannot silently disagree about
  what is broken — every marker has a live tracker, and every closed tracker
  leaves no stale marker behind.
- **Positive:** a reviewer can tell deliberate deferral from oversight at the
  site, and the known-bug PR-template row makes closing the loop a checkbox
  rather than a thing to remember.
- **Negative / cost:** fixing a bug now includes tracker + marker bookkeeping in
  the same PR, and a marker that spans both ports must be removed from both.

## Alternatives considered

- **Track known bugs only in the issue tracker (no in-code marker).** Rejected:
  the tracker records *that* something is deferred but not *where* — a reader at
  the workaround site has no breadcrumb, and the debt is invisible in the place
  it actually bites.
- **In-code markers only (no required tracking issue).** Rejected: a bare
  `FIXME` / `TODO` with no issue rots — no context, no owner, no resolution
  criteria, and nothing that ever closes it. The issue is the register; the
  marker is the pointer to it.
- **A separate `KNOWN_BUGS.md` register file.** Rejected: a hand-maintained list
  duplicates the issue tracker and drifts from it. GitHub issues already are the
  register, with the open/closed lifecycle built in.
