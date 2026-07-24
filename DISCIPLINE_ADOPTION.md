# Disciplined-project-seed adoption marker

sid adopts engineering-discipline conventions from the
[disciplined-project-seed](https://github.com/pdlourenco/disciplined-project-seed).
This file is the **marker**: it pins the seed version sid is synced to, records
per-artifact what was adopted / adapted / dropped, and keeps an append-only sync
log so each future seed release is a one-pass catch-up instead of archaeology.

**Seed version pinned:** `v0.4.1` (`49860ad`, 2026-07-22)

sid realizes the seed's *core* — a single binding specification (`spec/SPEC.md`)
with cross-language reference vectors (`testdata/`) — natively and more maturely
than the template; sid is recorded as the seed's fourth adopter ("governance
retrofit" profile). The table below covers the layers the seed grew after sid
forked the concept: agent governance, a decision trail, verification discipline,
and infra-as-code. **Adapted** means the seed convention was folded into sid's
layout (contracts in `spec/`, `CONTRIBUTING.md` at root), not copied verbatim.

## Adopted / adapted

| Seed convention | Seed ref | sid artifact | Status | Issue (PR) |
|---|---|---|---|---|
| Agent operating rules + PR lifecycle / two-session split | ADR-0009 | [`CLAUDE.md`](CLAUDE.md) | Adapted | #111 (#125) |
| PR template + known-bug row + ADR row | ADR-0012 | [`.github/pull_request_template.md`](.github/pull_request_template.md) | Adapted | #112 (#126, #128) |
| ADR infrastructure + parallel-track numbering | — | [`docs/decisions/`](docs/decisions/) | Adapted | #114 (#128) |
| Reviewer context (verification/validation modes) + pre-push self-review | ADR-0007 | [`docs/REVIEW_CONTEXT.md`](docs/REVIEW_CONTEXT.md) + `CONTRIBUTING.md` | Adapted | #115 (#129) |
| Architectural rationale doc | — | [`docs/DESIGN.md`](docs/DESIGN.md) | Adapted | #116 (#130) |
| Least-privilege workflow permissions + `GH_REPO` | seed PRs #21/#23 | `.github/workflows/*.yml` | Adopted | #131 (#132) |
| Contract-gate catalogue / drift-hardening doctrine | ADR-0011 | `CONTRIBUTING.md` + [ADR-0002](docs/decisions/ADR-0002-contract-artifact-hardening.md) | Adapted | #152 (#167) |
| Known-bug lifecycle (visible-debt markers + register) | ADR-0012 | `CONTRIBUTING.md` + [ADR-0003](docs/decisions/ADR-0003-known-bug-lifecycle.md) | Adapted | #153 (#168) |
| `analyses/` doc-type convention | ADR-0010 | [`docs/analyses/README.md`](docs/analyses/README.md) | Adopted | #154 (#169) |
| Local-CI runner | ADR-0004 | [`scripts/local-ci`](scripts/local-ci) | Adapted | #155 (#170) |
| Adoption marker | ADR-0013 | this file | Adopted | #151 (#171) |
| Branch-protection-as-code + drift workflow | seed BP-as-code | [`.github/rulesets/main.json`](.github/rulesets/main.json) + [`branch-protection.yml`](.github/workflows/branch-protection.yml) + [ADR-0005](docs/decisions/ADR-0005-required-ci-checks.md) | Adopted | — (#180) |
| Randomized-exploration / PBT (seeded MC: fixed-seed gate smoke + `SID_MC_CAMPAIGN`-gated sweep) | ADR-0014 | [`python/tests/test_ltv_uncertainty_calibration.py`](python/tests/test_ltv_uncertainty_calibration.py) | Adopted | #137 (#175) |

## Held

| Seed convention | Seed ref | sid artifact | Why held |
|---|---|---|---|
| `Verified by:` per-rule spec annotations | seed SPEC option | `spec/SPEC.md` | Re-derived after remediation Phase 3 completes (remaining: #137), landing with #113 — the June annotations predate the repo-wide review and over-claim. #127 held; un-defer worklist tracked on that PR. |

## Deferred (with triggers)

Real seed features parked with explicit triggers rather than dropped, same
discipline the seed applies to its own deferred items:

- **Label-as-code + issue templates** — *trigger:* issue volume past ~one screen, or a second regular contributor.
- **Release-as-code CHANGELOG gate** (ADR-0015) — *trigger:* next release-process change.
- **Project-level `CHANGELOG.md`** — *trigger:* per-language `RELEASE_NOTES.md` diverge, or a unified cadence emerges.
- **Traceability matrix** (seed SPEC option) — *trigger:* a formal V&V / ECSS process requirement.

## Sync log

Append-only. One entry per adoption wave or seed-release catch-up; on each seed
release, read the seed's `meta/CHANGELOG.md` between the pinned and new version,
take/skip each entry with a one-line reason, and bump the pin above.

- **2026-06-04 — initial adoption authored** against seed `18b10d1`. Stacked
  series #125–#132 (CLAUDE.md, PR template, ADR infra, REVIEW_CONTEXT, DESIGN,
  workflow perms; `Verified by:` as #127).
- **2026-07-22 — re-baseline planned** to seed `v0.4.1` (`49860ad`). The June
  series was seven weeks stale (seed released v0.4.0 + v0.4.1, five new ADRs);
  refresh plan `docs/plans/2026-07-22-seed-v0.4.1-adoption-refresh.md` (#150).
- **2026-07-24 — refreshed and merged** to `v0.4.1`. #125/#126/#128/#129/#130/#132
  refreshed with the v0.4.1 deltas (two-session split, known-bug row, parallel-track
  numbering, principle 9, runtime-compat) and merged. New increments adopted:
  contract-gate doctrine (#152), known-bug lifecycle (#153), analyses convention
  (#154), local-ci (#155), this marker (#151). `Verified by:` (#127) held for
  post-Phase-3 re-derivation. Pin set to `v0.4.1` (`49860ad`).
- **2026-07-24 — PBT convention adopted** (ADR-0014), moved deferred → adopted.
  Its named trigger (Phase 3 landing; #137's calibration-test tightening as the
  vehicle) fired and was swept in the same PR: `test_ltv_uncertainty_calibration.py`
  uses the two-tier shape — deterministic fixed-seed gate smoke + a
  `SID_MC_CAMPAIGN`-gated 500-trial λ-sweep out of the gate (#137, #175). With
  Phase 3 complete, `Verified by:` (#127) is now unblocked.
- **2026-07-24 — branch-protection-as-code adopted** (ADR-0005), moved deferred →
  adopted. The deferred trigger ("a drift incident") fired twice — #174 and #178
  both admin-merged while the MATLAB/Octave Lint check was red, undetected until
  #155's local-ci runner surfaced it (#179). The two cheap lint workflows now run
  on every PR with distinct contexts, and `.github/rulesets/main.json` +
  `branch-protection.yml` require them as code (admin bypass retained). Activation
  is a one-time maintainer step (add `RULESET_TOKEN`, run the workflow).
