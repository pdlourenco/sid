# ADR-0006: The roadmap is split along the contract/implementation layers

- **Status:** Accepted — _2026-07-29_
- **Deciders:** Pedro Lourenço (maintainer), Claude (analysis + implementation)
- **Related:** issue #195; [ADR-0001](ADR-0001-spec-is-the-contract.md);
  `CLAUDE.md` §4; [`docs/REVIEW_CONTEXT.md`](../REVIEW_CONTEXT.md)

## Context

sid has two development layers — a math-and-contracts layer (the spec, the
naming convention, the public function catalogue, project scope) and a
per-language implementation layer. The two roadmap files cut *across* those
layers rather than along them, and both had decayed into something other than
roadmaps:

- **Neither was a roadmap.** The bulk of `docs/roadmap.md` (265 of its 445 lines)
  was a checked-off execution log of the MATLAB v1.0 build-out (Phases 1–11 with
  day estimates); `docs/roadmap_python.md` was the same for the Python port
  (P1–P11, all ✅). The genuinely forward-looking content — reserved names,
  online/recursive COSMIC, the out-of-scope lists — was a small fraction of each.
- **Layer mixing.** `docs/roadmap.md` held the cross-language contract-adjacent
  material that `CLAUDE.md` §4, the `CONTRIBUTING.md` reviewer checklist, and
  `docs/REVIEW_CONTEXT.md`'s "catalogue drift" red flag all point at — the
  `sid + Domain + Method` convention and the function catalogue — while being
  titled *"for MATLAB/Octave"* and mixing in Octave compatibility rules and a
  MATLAB result struct. A Python or Julia contributor had no reason to read the
  document that defines the naming rule binding on their port.
- **A governance contradiction in a living document.**
  `docs/roadmap_python.md` stated *"The MATLAB code in `matlab/sid/` is the
  ground truth for numerical behaviour"* — directly contradicting ADR-0001 and
  `CLAUDE.md` §3 — and its porting workflow ("read MATLAB source, write
  Python") is the copy-by-copy-porting pattern `docs/REVIEW_CONTEXT.md` flags as
  a red flag. A contributor following the roadmap and a contributor following
  the governance docs would have been told opposite things.
- **Spec duplication.** The roadmap's "Result Struct" section duplicated the
  output-field contract owned by `spec/SPEC.md` — a drift surface with no
  verifier.

## Decision

`docs/roadmap.md` is the **language-neutral function catalogue and roadmap**: the
naming convention, one merged per-capability catalogue across all ports, and
forward-looking material only. Implementation-layer material lives in the
per-language guides, behaviour lives in `spec/SPEC.md`, and the executed phase
logs are archived under `docs/plans/` as historical records.

## Consequences

- The document that governance references for catalogue and naming is now
  actually language-neutral, so the rule binding on every port is stated once,
  in a file every port's contributors have reason to read.
- No living document claims MATLAB is ground truth. The claim survives only in
  an archived phase log, explicitly disclaimed in its header.
- The catalogue is a single merged table keyed by *capability*, so adding a port
  means adding a column, not forking a document. Adding a language no longer
  multiplies roadmap files.
- The roadmap no longer restates output fields; a reader after behaviour is sent
  to the spec. One fewer unverified duplication of the contract.
- **Cost:** per-capability rows must be updated in one more place when a port
  implements a function — the catalogue is not generated from code, and nothing
  verifies it. Catalogue drift remains a review-time check
  (`docs/REVIEW_CONTEXT.md`), not a mechanized gate.
- **Cost:** the MATLAB↔Python private-helper mapping table is dropped rather than
  relocated. It was already drifting, and helper decomposition is explicitly an
  implementation detail each port may choose differently.
- **Not closed by this ADR:** only the *ground-truth claim* is gone. The
  copy-by-copy porting workflow cited in Context above still stands as house
  policy in [`python/CONTRIBUTING.md`](../../python/CONTRIBUTING.md) §"Porting
  Workflow", which instructs contributors to read the MATLAB source and
  translate it, without mentioning `spec/SPEC.md`. Rewriting it spec-first is a
  governance change beyond this restructure's scope and is tracked separately
  (#196). Until then, a Python contributor is still pointed at MATLAB rather
  than the spec — this ADR narrows the contradiction, it does not eliminate it.

## Alternatives considered

- **Leave both files, fix only the ground-truth sentence.** Delete the offending
  line in `roadmap_python.md` and move on. Rejected because it treats the symptom:
  the layer mixing that put a cross-language contract inside a MATLAB-titled
  document, and the copy-porting workflow, would both remain — and the next port
  (Julia) would still have no obvious home for its column.
- **Keep per-language roadmaps, add a third neutral catalogue.** Three files: a
  shared catalogue plus one roadmap per language. Rejected because the two
  per-language files had almost no live content left once the phase logs were
  archived — it would institutionalise two near-empty documents and reintroduce
  the drift risk of parallel catalogues.
- **Generate the catalogue from source.** Derive the table by scanning
  `matlab/sid/*.m` and `python/sid/*.py`. Rejected as premature: it solves the
  update-in-one-more-place cost but not the layer split this ADR is about, and it
  needs a verifier and a CI leg of its own. A reasonable follow-up once the
  catalogue's shape has settled — the split has to happen first either way.
- **Delete the phase logs outright.** They are fully checked off and superseded.
  Rejected because they are the decision trail for how v1.0 was sequenced, and
  the repo's convention (`docs/plans/`) is to archive executed plans rather than
  drop them.
