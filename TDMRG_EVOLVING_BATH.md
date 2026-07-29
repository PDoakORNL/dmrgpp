# tDMRG evolving-bath: chain of assumptions and test map

This document exists because the evolving-bath work in
`cincuenta/src/ImpuritySolverNeqTdmrg.h` was built as a chain of small,
individually-verified links, each resting on an assumption established by
the previous one. After enough links, the chain itself becomes hard to hold
in your head. This is the map: for each link, what is assumed, why we
believe it, and which test (if any) actually protects it. When a link's
test is missing or a link is known-broken, that is stated explicitly rather
than glossed over.

See `/Users/epd/.claude/plans/fancy-painting-moon.md` for the full narrative
history (how each link was discovered/debugged). This doc is the *current
state* summary, not the story.

## How to read this

Each link has:
- **Assumption** — the specific claim this link depends on.
- **Why we believe it** — the empirical check that established it.
- **Protected by** — the Catch2 test case(s) that would fail if the
  assumption stopped holding. "NONE" means: nothing currently catches a
  regression here.
- **Status** — GOOD (currently holds, tested) / BROKEN (confirmed false,
  blocks downstream work) / UNTESTED (plausible but never checked).

All tests below live in `cincuenta/src/tests/test_ImpuritySolverNeqTdmrg.cpp`
unless noted. Run via `ctest -R ImpuritySolverNeqTdmrg` — **never** by
invoking the compiled binary directly with no filter; `catch_discover_tests`
registers each `TEST_CASE` as its own ctest process, and running all cases
in one manual process has produced file-collision false failures that don't
occur under ctest (see Link 11).

---

### Link 1 — Chained restart-based `TimeEvolve` segments are safe across a Connectors change

**Assumption:** nothing in DMRG++'s restart/checkpoint validation pins a new
segment's `Connectors=`/`potentialV=` to match the previous segment's.

**Why we believe it:** read `Checkpoint.h` directly — it only validates
`IsComplex` and `operatorsPerSite` across a restart.

**Protected by:** indirectly, every test below that chains more than one
segment with different Connectors (Links 5-8). No standalone test for this
specific claim.

**Status:** GOOD, but see Link 9 — this claim is about the restart being
*permitted*, not about the new Hamiltonian actually being *used*. Those
turned out to be different claims.

---

### Link 2 — `RestartSourceTvForPsi` selects which TV becomes the new `|gs>`

**Assumption:** `RestartMappingTvs=[i0,i1,...]` is positional-by-TV-index
remap; `RestartSourceTvForPsi=k` separately selects which of the *old* run's
TVs seeds the *new* run's `|gs>` reference.

**Why we believe it:** `RestartStruct.h:33-38` (mapping semantics) +
`TargetingCommon.h:304` (consumption, itself flagged
`"FIXME TODO WARNING: Need better spec for TvForPsi"` in upstream code —
noted as a real upstream rough edge, not our bug).

**Protected by:** `test_ImpuritySolverNeqTdmrg.cpp`'s column-0 chaining test
(`solveChainedColumn0` vs monolithic `solve()`, n=1 and n=2, exact match).

**Status:** GOOD for what it does (selecting the seed state). Does NOT
establish that the new segment's Hamiltonian is respected post-restart —
see Link 9.

---

### Link 3 — A single `FiniteLoops` row per segment is insufficient; need two (there-and-back)

**Assumption:** one `[@auto,m,2]` row does not reliably sweep back to site 0
(the measurement site); two rows do.

**Why we believe it:** empirically — one-row segments produced logs with a
sweep from site 4→0 (worked) then a later segment sweeping 1..5 only (missed
site 0 entirely, gave silently-wrong all-zero output).

**Protected by:** every multi-column/multi-step test in this file
(`buildStepInput` always emits two rows now). No test specifically asserts
"if this regressed to one row, X would fail" — it's baked into
`buildStepInput` itself, so a regression would need someone to edit that
function.

**Status:** GOOD, structural (encoded in code, not a runtime check).

---

### Link 4 — Each segment's log has two site-0 measurement occurrences; which one is "the" answer depends on whether this is the column's first advance

**Assumption:** a segment restarting from a **plain single-TV** checkpoint
(column's first advance since birth) reports the wanted post-advance value
FIRST, throwaway SECOND. A segment restarting via `RestartSourceTvForPsi`
from a **prior chained segment** reports the wanted value SECOND (the
`takeLast` case).

**Why we believe it:** diffed raw logs directly at U=0 and U=0.5 (ruling out
U=0 degeneracy as an alternative explanation) — pattern was identical at
both, and switching `parseSingleMeasurement`'s `takeLast` flag per this rule
made n=1 and n=2 match `ImpuritySolverNeqExactDiag` exactly.

**Protected by:** `solveChainedColumn0` test (24/24 assertions, both U
values) + `computeFullGrid` test (48/48, both U values).

**Status:** GOOD.

---

### Link 5 — A column's first-advance-since-birth FIRST measurement occurrence IS the equal-time diagonal `G(born,born)`

**Assumption:** this diagonal value is a free byproduct of the existing
off-diagonal chaining, not something requiring a separate measurement
scheme.

**Why we believe it:** matches `ImpuritySolverNeqExactDiag`'s `G(0,0)`
exactly for column 0; generalizes because every column's first advance is
structurally identical to column 0's.

**Protected by:** `computeFullGrid` test's diagonal assertions (part of the
72/72 total).

**Status:** GOOD. Gauge-invariant by construction (P1/P2 both descend from
the same birth state), so no phase correction is needed for this value
specifically.

---

### Link 6 — The `<P2.last|P2>` phase-gauge correction is unnecessary and was actively wrong

**Assumption:** `<P2|c|P1>`'s bra/ket share the same loaded reference phase
and cancel it already; dividing by `<P2.last|P2>` on top injects a spurious
rotation whenever that overlap isn't trivially 1 (which it stopped being
once an eps-split extended geometry made the reduced density matrix more
degenerate).

**Why we believe it:** (a) `<P2.last|P2>` was confirmed constant across
sweep positions AND across `t=0`/`t=dt` rows within a segment — a genuine
dynamical phase would differ between those rows, it didn't; (b) extended-
vs-unextended GS sector energies differ by exactly `-2*eps*L`, confirming
the impurity block's physics is undisturbed and the phase is pure
DMRG-basis arbitrariness.

**Protected by:** `computeFullGridWithInertSecondBath` test (24/24,
`Vplus=0` inert-spectator check) — this is the first test that exercises a
genuinely non-trivial `<P2.last|P2>`, so it's the one that would catch a
regression here.

**Status:** GOOD (correction removed from `advanceColumn`; `Column::gaugePRaw`/
`gaugeHRaw` fields gone).

**Known untested sibling:** `solveChainedColumn0`'s own separate
`gaugeP`/`gaugeH` correction was never tested against a non-trivial-gauge
case and likely has the same latent bug. Also, the monolithic `solve()`
path's `applySignFlip`/`applyGlobalPhase` have never been tightly tested
against a genuine seed≠evolve-Hamiltonian quench (production input
`inputNeqTsujiTdmrg.ain` runs exactly that, but only passes via a
self-consistent stored-reference comparison, not an analytic one).
**Status: UNTESTED**, not currently blocking anything.

---

### Link 7 — Complex-valued `Connectors=` entries are correctly parsed and Hermitian-conjugated by the live engine

**Assumption:** Ainur's `<real>i<imag>` literal syntax, fed through
`usecomplex`, produces a correct complex hopping bond (not silently
truncated to its real part, not mis-conjugated).

**Why we believe it:** standalone `dmrg` binary runs — GS energy of a single
hopping bond depends only on `|V|`; confirmed exactly for `V` real, complex,
and pure-imaginary (`-0.5831`/`-0.583095`/`-0.7`, all matching `-|V|`).

**Protected by:** NONE as a Catch2 regression test — this was a one-off
standalone smoke test in `build/tmp/complex_connector_smoke/`, not committed
as a test case. `buildConnectorsStrWithSecondBath`'s complex-literal
formatting path IS exercised indirectly by every `NeqBathRank=1` test with
nonzero `Vplus`, but there is no test that would specifically catch, e.g., a
future refactor of `formatComplexLiteral` reintroducing a sign or
conjugation error.

**Status:** GOOD but **UNPROTECTED** — worth promoting the standalone smoke
test into a real Catch2 case if this code is touched again.

---

### Link 8 — eps-split GS potential correctly seeds the second bath's L-occupied/L-empty configuration, provided eps is scaled to the largest other coupling present

**Assumption:** `eps = 5 * max(|first-bath hoppings|, |bathEps| range, U)`
reliably seeds the intended occupation; a fixed small `eps` (e.g. 0.01) does
not, once real couplings are present.

**Why we believe it:** standalone `dmrg` runs directly measuring `⟨n_p⟩`:
`eps=0.01` failed once a `V=0.5` first-bath hopping was present (electron
went to the coupled impurity/bath bonding orbital instead); `eps=1.0`
(comparable to V) worked.

**Protected by:** `measureSecondBathOccupations` test, tag
`[Phase2][SecondBathSeeding]` — asserts both the working case and the
documented-failure case through the real C++ code path.

**Status:** GOOD, protected.

---

### Link 9 — A column's advance segments beyond the first respect that segment's own `Connectors=` value — **BROKEN, confirmed 2026-07-29**

**Assumption (as originally held):** since each segment's Ainur input text
is freshly generated by `buildStepInput` on every call (confirmed by
reading the code — there is no cached/reused input text), a segment
restarting via `RestartSourceTvForPsi` from a prior chained segment
(`takeLast=false`) evolves under **that segment's own**, freshly-declared
Connectors, exactly like the first-advance case does.

**What we actually found:** it does not. `diagnosticSecondAdvanceConnectors`
births a column, advances it once with `connectorsStep1`, then advances it
again with `connectorsStep2`:
- Holding `connectorsStep1` fixed and varying `connectorsStep2` (even a
  large change, `0.26` → `0.52-0.008i`) — **result unchanged**
  (`(0.493548,-0.063098)` both times).
- Positive control: holding `connectorsStep2` fixed and varying
  `connectorsStep1` instead — **result changes**
  (`(0.493548,-0.063098)` vs `(0.489516,-0.062906)`).

Conclusion: a column's second-and-later advances are frozen at whatever
Hamiltonian the FIRST advance used; a newly-declared Connectors value on any
later segment is silently ignored by the engine. This is a genuine,
previously-undiscovered gap in DMRG++'s restart machinery (likely in how
`RestartSourceTvForPsi`/checkpoint restore interacts with Hamiltonian/block
operator construction — not yet root-caused at that level), not a plumbing
mistake in this project's Ainur-generation code.

**Why Phase 1 never caught this:** Phase 1 always used a **static** bath —
the same Connectors value at every step of a run. A "frozen at the first
advance" bug is invisible when every advance would have used the same value
anyway. It only becomes visible once Connectors genuinely vary step-to-step,
which is exactly what Phase 2's self-consistent `Vplus(n,p)` requires.

**Protected by:** `diagnosticSecondAdvanceConnectors` +
`"...positive control"` tests, tag `[Diagnostic]` — these are TEMPORARY,
investigation-only test cases (not permanent regression coverage; should be
replaced by a proper permanent regression test once this is fixed, asserting
non-frozen behavior with `CHECK`/`Approx` rather than just printing).

**Status:** BROKEN. This is the reason Task 15's full gate
(`NeqBathRank=1` tDMRG vs GBEK) currently fails at `n=2` — column 0's second
advance can't pick up `Vplus(2,·)`. **Blocks all of Phase 2's
self-consistent correctness for `nT ≥ 2`.** Everything downstream of this
link (Links 10-11) is built on top of infrastructure that is real but not
yet correct for more than one time step.

---

### Link 10 — `decomp_`/`prepareTimeStep` wiring and the "truncated batch recompute" design

**Assumption:** fully rebuilding the column fan-out from t=0 on every
`computeGimp(gimp,n)` call, reading `Vplus(k,p)` at its current state, is a
correct (if expensive) way to keep a corrector loop simple — no
propagated-through/rewind bookkeeping needed, since nothing is cached across
calls.

**Why we believe it:** this part of the design is *logically* sound
independent of Link 9 — the recompute-from-scratch strategy is not itself
wrong. But it inherits Link 9's bug: since the recompute chains segments the
same way Phase 1 does, its `n≥2` columns are equally frozen at whatever
Connectors were used for their own first advance, regardless of `decomp_`'s
current state.

**Protected by:** the wiring/dispatch itself (`solve()`/`computeGimp()`/
`gimp()` routing on `neqBathRank_`) is protected by the
`"...NeqBathRank=1 self-consistent wiring runs end-to-end without crashing"`
and `"...can be driven through NeqDmftSolver's own solve()"` tests — both
only check `isfinite`, not correctness. **The actual correctness claim (does
`Vplus` reaching the engine change the physical answer) is what Task 15's
still-failing full gate is meant to protect, and it currently does not
pass.**

**Status:** wiring GOOD (crashes nothing), correctness BROKEN pending
Link 9's fix. Do not trust `NeqBathRank>0` tDMRG results for `nT≥2` until
Link 9 is resolved and Task 15 passes.

---

### Link 11 — Test-isolation: distinct file/checkpoint roots per call within one process

**Assumption:** reusing the same `RootOutputname`/internal chain-root prefix
across multiple `DmrgRunner` calls in one process (whether across
`TEST_CASE`s or across repeated calls to the same method) can silently
corrupt results; each call needs a distinct prefix.

**Why we believe it:** found via two concrete collisions — two back-to-back
`computeFullGridWithInertSecondBath` calls sharing a root, and two
`NeqBathRank=1` tests sharing a `RootOutputname` — both fixed by giving each
call a distinct root/suffix and confirmed clean afterward.

**Protected by:** the fix is structural (explicit `rootSuffix` parameters,
distinct `RootOutputname`s, `scCallCounter_`-based per-call chain roots in
`fillSelfConsistentRow`) rather than test-asserted. **Not a CI risk**:
`catch_discover_tests` isolates every `TEST_CASE` into its own ctest
process, so this class of bug cannot occur under `ctest -R
ImpuritySolverNeqTdmrg` regardless — it only ever showed up when manually
running the whole compiled binary in one process, which is not how CI or
`ctest` invoke it.

**Status:** GOOD for CI. A convention to keep in mind for future manual
debugging sessions: don't take "manually running the whole binary" failures
as ground truth for this file.

---

## Summary table

| # | Link | Status |
|---|------|--------|
| 1 | Chained restarts across a Connectors change are permitted | GOOD |
| 2 | `RestartSourceTvForPsi` selects the new `\|gs>` seed | GOOD |
| 3 | Two `FiniteLoops` rows needed per segment | GOOD (structural) |
| 4 | First-vs-second measurement occurrence per segment type | GOOD |
| 5 | First-advance-since-birth's first occurrence = diagonal | GOOD |
| 6 | `<P2.last\|P2>` gauge correction removed | GOOD (siblings UNTESTED) |
| 7 | Complex `Connectors=` parsed/conjugated correctly | GOOD, UNPROTECTED |
| 8 | eps-split seeding, scaled to largest coupling | GOOD |
| 9 | Later advances respect their own segment's Connectors | **BROKEN** |
| 10 | Truncated-batch recompute + `decomp_` wiring | wiring GOOD, correctness **BROKEN** (depends on #9) |
| 11 | Per-call file-root isolation | GOOD |

**Bottom line:** everything through Link 8 is solid and tested. Link 9 is a
confirmed, real bug in the underlying chained-restart architecture — it
blocks Task 15 (the Phase 2 full correctness gate) and means **no
`NeqBathRank>0` tDMRG result for `nT≥2` should be trusted** until it's fixed.
Links 10 is built correctly on top of a foundation that isn't correct yet;
it will not need rework once Link 9 is fixed, just re-verification.
