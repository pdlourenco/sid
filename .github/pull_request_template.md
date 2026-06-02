<!-- Keep this short. Delete sections that don't apply. -->

## Summary

<!-- What this PR does and why, in a few sentences. -->

## Spec sections touched

<!-- Cite the spec/SPEC.md §X.Y (and spec/EXAMPLES.md) sections this PR
     implements or modifies. Write "none — no algorithmic change" if it
     touches no contract. If a contract changed, state that the spec was
     updated first and every implementation side updated in the same PR. -->

## Type of change

- [ ] Algorithmic change (touches a spec contract — spec updated first)
- [ ] New feature / function (non-contract)
- [ ] Bug fix
- [ ] Docs / tooling / CI only

## Checks

- [ ] Implementations derive behaviour **independently from the spec**, not ported copy-by-copy between languages
- [ ] Tests written against spec requirements, not against the other language's current output
- [ ] Shared helpers touched? Audited every caller against the relevant spec section
- [ ] New defaults / bounds / thresholds are covered by `spec/SPEC.md`
- **Pre-push review:** <!-- "no findings" or "flagged X, fixed in <sha>" — see CLAUDE.md §2 -->
- **Local CI:** <!-- "green" or "<job> failed, fixed in <sha>" -->
