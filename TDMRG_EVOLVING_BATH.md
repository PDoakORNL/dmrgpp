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

### Link 9 — Each segment applies exactly one dt-advance under its own `Connectors=` value — **FIXED, 2026-07-29/30, engine-level `maxAdvances` cap**

**Assumption (as originally held, WRONG):** a segment restarting via
`RestartSourceTvForPsi` from a prior chained segment (`takeLast=false`)
evolves under **that segment's own**, freshly-declared Connectors, applying
exactly one dt-advance, exactly like the first-advance case does.

**First (refuted) hypothesis: "the engine ignores Connectors on continuation
segments".** `diagnosticSecondAdvanceConnectors` (birth a column, advance
once with `connectorsStep1`, again with `connectorsStep2`) showed the
harvested value insensitive to `connectorsStep2` and sensitive only to
`connectorsStep1`. This looked like a frozen Hamiltonian. **Refuted** by a
third-advance discriminating test (`diagnosticThirdAdvanceConnectors`,
Connectors `C1,C2,C3`, harvest at step 3): varying `C2` (the *middle*
segment's Connectors, held fixed at both ends) DID change the step-3
harvest. A frozen-at-first-advance engine could not produce that — the
value clearly responds to *some* later segment's Connectors, just not the
one intuition expected.

**Second (refuted) hypothesis: one-segment indexing lag.** The natural
explanation for "step-3 harvest tracks C2, not C3" is an off-by-one: a
continuation segment's occurrence #1 (the measurement kept for
`takeLast=false`, per Link 4) is a readout of whatever the *previous*
segment produced, so `harvest(m)` tracks `connectors` passed at call
`m-1`. Fixing this looked like a pure index shift (`vMidConnectors(k+1)`
instead of `vMidConnectors(k)` for continuation calls). Attempting the
shift immediately hit a real conflict: a single call's Connectors argument
would need to simultaneously satisfy the direct requirement for its OWN
row (first-advance calls, where harvest = argument directly, no lag) and
the fed-forward requirement for the NEXT row — two different required
values whenever `Vplus` genuinely varies. That conflict was the tell that
the model was still wrong, not just the fix.

**Root cause, confirmed via the engine's own instrumentation (not
guessed):** grepping the `Steps without advance` trace that
`ApplyOperatorExpression::evolve` already prints
(`dmrg/Engine/ApplyOperatorExpression.h:795-798`) in the generated logs
showed fires happening at `site=0` and `site=7` (both sweep BORDERS) within
a single two-row segment, each resetting the counter. Reading
`ApplyOperatorExpression.h:769-773` explains why: unless `SolverOptions`
includes `advanceUnrestricted`, a TimeEvolve advance can only fire when the
sweep is AT A BORDER, regardless of `advanceEach`'s count. With
`advanceEach = nsitesExt-2` (this project's convention), each of the two
`FiniteLoops` rows independently accumulates enough count to fire once, at
its own border. **Net effect: every segment applies TWO real dt-advances
under its single, segment-level `dir0:Connectors=` value, not one.**

Direct confirmation from the measurement logs' own printed elapsed-time
field across a 3-call static-Connectors chain (`testTdmrgChaindiagcalib6_*`):
imaginary parts `0.0, -0.031695, -0.063098, -0.093923, -0.123895,
-0.152755` — linear at ≈`0.0308`/advance, six consecutive advances across
three calls. So: occurrence #1 of call `m` sits at `2(m-1)` advances,
occurrence #2 at `2m-1` advances.

**Consequences, now precisely characterized (not hand-waved):**
- The `n=2` gate failure is a pure **bath** error, not a time error: call
  2's occurrence #1 lands at exactly 2 advances (the physically correct
  elapsed time for row 2), but BOTH of call 1's advances used
  `vMid(1)` — the interval `[t1,t2]` silently got `[t0,t1]`'s bath.
- **`nT≥3` is wrong in TIME, not just bath**: call `m`'s harvest sits at
  `2(m-1)` advances, mislabeled as step `m-1`. No prior test caught this
  because no passing test ever drove a column past its *second* advance
  (Phase 1 at `nT=2`: column 0 gets calls at `k=1,2` only; column 1 gets
  one advance; column 2 only a throwaway).
- **Why Phase 1 never caught any of this:** Phase 1 always used a
  **static** bath (same Connectors at every step). Two advances under an
  unchanging Hamiltonian just means "twice as much of the same physics" —
  invisible unless you're checking absolute step count, which Phase 1's
  own `NtNeq=2` tests never had reason to.

**Two candidate fixes tried and empirically refuted, both via direct
calibration (not inspection):**
- Scaling `advanceEach` upward under the SAME border-gating: plateaued —
  every trial 7 through 12 gave the identical "2 advances" result (7-12 all
  fail to trigger within one row alone, but the counter carries over
  across the segment's two rows and still fires exactly twice by the
  segment's end); no value produced "1 advance".
- Adding `advanceUnrestricted` (removing the border gate) and scaling
  `advanceEach` as a pure count: also plateaued — trials 11-13 gave "2
  steps total across 3 calls", trials ≥14 gave "0 advances at all", never
  converging on "1 advance/segment" for any tested value. `firstSeeLeftCorner_`
  (a precondition for the counter to start incrementing at all, per
  `ApplyOperatorExpression.h:783-790`) appears to cost the FIRST segment
  since restart its ability to fire, shifting where subsequent segments'
  accumulated count lands — `advanceEach` cannot express "exactly one
  advance per two-row segment" under either gating mode with this engine.
- A one-`FiniteLoops`-row-per-segment alternative (bypassing the two-row
  requirement entirely) was also tried: site 0 was never reached at all in
  the extended geometry for the second call (0 measurement occurrences) —
  confirms the ORIGINAL reason two rows exist (a single row's sweep isn't
  guaranteed to reach the measurement site), so this is dead, not viable.

**Fix chosen and implemented: Option 1 (engine instrumentation), per explicit
user direction ("let's go with scope choice 1").** Added an opt-in
`maxAdvances` cap, fully backward compatible (0 = unlimited, the default,
preserving every existing behavior byte-for-byte):
- `GroupOfOneTimeEvolutions.h`'s `OneTimeEvolution` gained an
  `advancesSoFar_` counter (lifetime of the `OneTimeEvolution`, i.e. spans
  every `FiniteLoops` row of the segment that owns it — not reset per row),
  incremented each time `advanceTime` actually fires.
- `NonLocalForTargetingExpression.h`'s `TimeParams` gained a `maxAdvances`
  field (parsed from `TimeEvolve{...,maxAdvances=N}`); `advanceInTimeOrNot`
  ANDs `!advanceCapped` (`advancesSoFar() >= maxAdvances`) into its firing
  condition, on top of the existing `advanceOnlyAtBorder` gate — so a
  segment can still only fire AT a border (unchanged), but is refused a
  SECOND fire once it's used its quota, regardless of how many more borders
  the segment's two `FiniteLoops` rows cross.
- `cincuenta/src/ImpuritySolverNeqTdmrg.h`'s `SecondBathExt`/`buildStepInput`
  thread `maxAdvances=1` through every self-consistent-path segment (the
  static/Phase-1-only path is unaffected, `maxAdvances` defaults to 0 there).
  `advanceColumn`'s `takeLast` logic was widened: under the cap, a
  continuation segment now ALSO takes the last occurrence (there is only
  ever one real advance now, so the old "continuation ⇒ first occurrence"
  rule — which existed specifically to skip the second, over-advanced
  occurrence — no longer applies).

**Verified directly against the engine, not just by re-running the gate:**
- Engine trace (`NonLocalForTargetingExpression`'s own `"Steps without
  advance"` print) for a `maxAdvances=1` segment shows the counter reach a
  border, fire once (elapsed time jumps `0.0→0.1`), and then continue
  incrementing all the way to the segment's second border WITHOUT a second
  jump — the cap visibly holds.
- `diagnosticThirdAdvanceConnectors` re-run with `maxAdvances=1`: step-3
  harvest now differs when C3 (the segment's OWN Connectors) differs, with
  C1/C2 held fixed (`resultA=(0.485552,-0.093923)` vs
  `resultB=(0.483213,-0.093723)`) — confirms the one-segment lag described
  above is gone; a segment's harvest now tracks its own Connectors, not the
  previous segment's.
- `vMidConnectors(k)` (the midpoint-averaging of `Vplus(k-1,·)` and
  `Vplus(k,·)` for the segment covering `[t_{k-1},t_k]`) was independently
  re-derived from raw debug traces during this verification and confirmed
  arithmetically correct — this rules out midpoint-averaging as a
  contributor to any remaining discrepancy (see Link 12).

**Protected by:** `"...maxAdvances=1 gives exactly one advance per segment"`
and `"...maxAdvances=1 -- step-3 harvest now tracks its OWN Connectors
(C3)"` tests, tag `[Diagnostic]` — investigation-only, not yet promoted to
permanent regression coverage. `diagnosticSecondAdvanceConnectors` and
`"...positive control"` remain from the root-causing phase, now stale
(they document the OLD 2-advances-per-segment behavior) — candidates for
removal or repurposing once Link 12 resolves and Task 15 passes.

**Status:** FIXED and verified in isolation. This closes the root cause
Link 9 originally identified. Task 15's gate, however, still does not pass
— see Link 12, a narrower, still-open discrepancy uncovered only after this
fix removed the coarser 2-advances-per-segment error.

---

### Link 10 — `decomp_`/`prepareTimeStep` wiring — originally a "truncated batch recompute", **replaced 2026-07-30 by an incremental design (Task 17)**

**Original assumption (superseded, kept for record):** fully rebuilding the
column fan-out from t=0 on every `computeGimp(gimp,n)` call, reading
`Vplus(k,p)` at its current state, was a correct (if expensive,
`O(nT^3*neqDmftIter)`) way to keep a corrector loop simple — no
propagated-through/rewind bookkeeping needed, since nothing was cached
across calls. This part of the design was *logically* sound independent of
Link 9 (the recompute-from-scratch strategy was not itself wrong), but it
inherited Link 9's bug while that was still open.

**Replaced by an incremental design (Task 17, `fillSelfConsistentRow`):**
`scColumns_` now persists across calls for the life of a solve.
Correctness rests on the same fact Link 12's investigation established
directly against source: `NeqBathDecomposition::update(n,.)` only ever
mutates row `n` of its `V_` — rows `<n` are frozen forever once step `n`'s
own corrector loop finishes. So a column already advanced through step `k`
never needs redoing for a later `n>k`; `prepareTimeStep`'s new rollback
(a one-level cursor-field snapshot/restore per column, NOT a bare
watermark decrement — `advanceColumn` destructively overwrites a column's
on-disk cursor, so simply relabeling an integer the way GBEK's own
`propagatedThrough` does would leave the cursor pointing at the wrong
checkpoint) undoes exactly the current row's own last advance whenever a
corrector refines `Vplus(n,.)`. New complexity:
`O(nT^2*(1+neqDmftIter))` — a full factor of `nT` better. See
`fancy-painting-moon.md`'s "Task #17 scope" section for the full design
derivation, and `project_tdmrg_evolving_bath` memory for the implementation
account (one real bug found and fixed along the way: a dangling reference
to `scColumns_[0]` held across an `emplace_back` call that could reallocate
the vector).

**Protected by:** the wiring/dispatch itself (`solve()`/`computeGimp()`/
`gimp()` routing on `neqBathRank_`) by the `"...NeqBathRank=1
self-consistent wiring runs end-to-end without crashing"` and `"...can be
driven through NeqDmftSolver's own solve()"` tests (both exercise the real
corrector loop, not just `isfinite` — they are what caught the dangling-
reference bug during implementation). The actual correctness claim (does
`Vplus` reaching the engine change the physical answer, and does the
incremental rewrite reproduce the original batch recompute bit-for-bit) is
Task 15's full gate, re-verified bit-identical after every one of Task 17's
five implementation steps.

**Status:** wiring GOOD, incremental design GOOD (bit-identical to the
prior batch recompute, confirmed after each implementation step). Task 15's
gate still passes 23/24 assertions, with the same one small characterized
residual as before this refactor (see Link 12) — unaffected by it, as
expected for a pure performance change. `NeqBathRank>0` tDMRG results for
`nT≥2` remain trustworthy to approximately the `1e-4` level demonstrated by
Task 15's gate, now at a fraction of the previous computational cost.

---

### Link 12 — Task 15 gate still fails at `(n=2,j=1)` after Link 9's fix — **RESOLVED, 2026-07-30: sticky `std::cout` precision leak in `NeqDmftSolver.h`**

**Observation:** with `maxAdvances=1` correctly wired (Link 9 fixed and
independently verified), Task 15's gate improved — `(n=2,j=0)` now passes —
but `(n=2,j=1)` still failed: `retarded.imag` expected `-0.98364`, got
`-1.0`; `lesser.real` expected `-0.063`, got `-0.1`; `lesser.imag` expected
`0.4918`, got `0.5`. The errors were small (a few percent), not the
order-of-magnitude symptom Link 9 produced — a materially different,
narrower discrepancy.

**Ruled out, with direct evidence, before the actual cause was found:**
- **Not a missing bath update, midpoint-averaging bug, or missing/extra
  advance.** The segment's generated Ainur input, hand-reconstructed
  `vMidConnectors(2)` arithmetic, and the engine's own `"Steps without
  advance"` trace were all independently confirmed correct for the exact
  failing segment.
- **Not column 1's birth checkpoint killing Connectors sensitivity.** A new
  diagnostic (`diagnosticColumn1FirstAdvanceConnectors`) built the same
  birth-then-first-advance segment in isolation and confirmed it DOES
  respond to its own Connectors (`(0.498382,-0.031689)` vs
  `(0.497370,-0.031684)` for two different post-birth Connectors values) —
  and critically, this isolated reproduction gave a value close to
  `0.4974`, NOT the real run's `(0.5,-0.0)`. Same method, same geometry,
  same `maxAdvances=1` — different answer. That mismatch was the actual
  discriminator.
- **Not eps-split seeding for this specific (5-site-bath) config.**
  Task 8's own validation used a different (1-bath-site) config;
  re-measuring `⟨n_p⟩` with the gate's own `bathParams`/`eps=3.0` (a new
  permanent-style check) confirmed correct seeding (`occ[0]=2.0`,
  `occ[1]=0.0`) — not the cause.
- **`applySignFlip` confirmed clean.** Column 1's `ggtRaw` has exactly one
  key (its single off-diagonal entry), so the function's `size()<2`
  early-return makes it a structural no-op there; it cannot be responsible
  for `(n=2,j=1)`'s error.

**Root cause, found by re-examining the raw log text itself:** the suspect
log's measurement lines showed only ONE decimal digit —
`0 (0.5,0.0) 0.0 <P2|c|P1> ...`, `0 (0.5,-0.0) 0.1 <P2|c|P1> ...` — not
truncated/stale (a fresh, correctly-rewritten file), but genuinely
LOW-PRECISION TEXT. `TargetingCommon::test()` (`dmrg/Engine/
TargetingCommon.h:864-865`), the function that prints every in-situ
measurement line, writes to the global `std::cout` with whatever format
flags/precision `std::cout` currently has — it does not set its own.
`cincuenta/src/NeqDmftSolver.h`'s `solve()` prints a per-step wall-clock
timing line using `std::cout << std::fixed << std::setprecision(1) << ...`
directly on `std::cout`, and never restores the prior state. `std::fixed`/
`setprecision` are STICKY on the stream object they're applied to — they
persist across `DmrgRunner`'s per-call log-file redirects, which redirect
the underlying buffer, not reconstruct the stream. So once step 1 of
`NeqDmftSolver::solve()`'s outer loop finished and printed its timing line,
EVERY subsequent in-situ measurement logged anywhere in the process for
the rest of that `solve()` call — including all of step 2's segments —
was silently rendered at `std::fixed`, 1 decimal digit. `(0.4974...,
-0.0317...)` round-trips through that formatting as exactly `(0.5,-0.0)`.
This explains every earlier observation: it's why isolated diagnostics
(which never call `NeqDmftSolver::solve()`) always showed full precision
while the real gate run didn't, and why it looked like a physics bug for
two sessions.

**Fix:** format the timing strings into local `std::ostringstream`s instead
of setting precision on `std::cout` directly (`NeqDmftSolver.h`). Confirmed:
`(n=2,j=1)` now matches GBEK to five decimal places
(`tDMRG (-0.031684,0.497352)` vs `GBEK (-0.031743,0.497350)`), a ~50x
improvement over the prior ~2-4% discrepancy. 23/24 assertions in Task 15's
gate now pass.

**Protected by:** the fix itself is structural (format into a local stream,
never touch the shared stream's state). `diagnosticColumn1FirstAdvanceConnectors`
and the gate-config eps-split occupation check, both added while narrowing
this down, are kept as permanent `[Diagnostic]`-tagged regression coverage.

**Status:** RESOLVED. One small, characterized residual remains (see the
note below) — NOT believed related to this bug, and not blocking use of
the evolving-bath machinery.

**Residual note — `(n=2,j=0)` `retarded.real` = `0.000118` vs GBEK's `~0`,
margin `1e-4`, still red.** `(n=2,j=0)` is the only entry in this gate
assembled from column 0's SECOND chained advance (two restarts via
`RestartSourceTvForPsi`, vs one birth-restart for columns 1/2) — the only
entry where `applySignFlip` is even structurally live (`ggtRaw` has 2 keys
there, confirmed a no-op by inspection: `sum_norm` for two nearly-equal
values is far above `0.1*diff_norm`, so it never flips). The same row shows
a consistent, opposite-sign `5.9e-5` offset in `lesser.real` at both `j=0`
and `j=1` — comparable magnitude on both, reads as a small accumulated
offset across the row rather than a localized defect. Checked whether it's
ordinary bond-dimension truncation by re-running the gate with
`InfiniteLoopKeptStates` raised 100→300: the residual was bit-identical
(`0.000118` both times) — but this check is weaker than it looks for this
system (`~4900`-dimensional Hilbert space, `LanczosCore` reports
`mat.rank=551`; 100 kept states is very plausibly already exact here, so
"no change at 300" doesn't distinguish truncation from something else).
`TridiagEps` is not binding either (engine trace shows `actual eps=3e-13`).
Ran out of cheap discriminating checks; both sides of the comparison
(tDMRG's DMRG+Krylov chain, GBEK's own ED/Lanczos) are themselves
approximate, and there is no evidence of a defect, only an unexplained
small residual on the exact-zero reference with no relative-error slack to
absorb it. **Deliberately NOT widening the gate's margin** (currently
`1e-4`) to paper over this without justification — left red, characterized,
for a future session to either explain or accept with a justified number.

---

### Link 13 — Hazard: sticky `std::cout`/`std::cerr` format state leaks across `DmrgRunner` calls in one process

**Not a link in the assumption chain — a standing hazard for anyone
debugging this codebase's in-process, multi-`DmrgRunner`-call solvers
(tDMRG's chained columns, GBEK's per-step calls, or any future one).**
`DmrgRunner`'s per-call log redirection (`dealWithConsoleOutput`) redirects
`std::cout`'s underlying buffer to a new file each call — it does NOT
reconstruct the stream object, so any format state (`std::fixed`,
`setprecision`, `std::hex`, etc.) ever applied directly to `std::cout` or
`std::cerr` PERSISTS across every subsequent call in the same process,
silently changing how later, unrelated log files render numbers.

Link 12 was exactly this: `NeqDmftSolver.h`'s own progress-printing line
set `std::cout` to `std::fixed`+`setprecision(1)` and never restored it,
silently truncating every in-situ measurement logged for the rest of that
`solve()` call to 1 decimal digit — which was misread as a ~2-4% physics
discrepancy for two full sessions before the log text itself was checked
character-by-character.

**First-thing-to-check, before any physics investigation:** if a
log-parsed numeric value looks implausibly round (`0.5`, `-0.0`, `1.0`,
etc.) in a context where it plausibly shouldn't be, open the RAW log file
and look at the literal text — not a value already parsed into memory —
and check how many digits are actually printed. If it's suspiciously few,
grep the calling code (and everything upstream of it in the same process)
for direct `std::cout`/`std::cerr` manipulator use (`std::fixed`,
`std::setprecision`, `std::cout.precision(...)`) that isn't scoped to a
local stream. The established, safe idiom already used almost everywhere
else in `dmrg/Engine/` is `PsimagLite::OstringStream msgg(std::cout.precision())`
— construct a LOCAL stream, format into it, print its `.str()` — never set
format flags on the shared global stream directly.

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
| 9 | Each segment applies exactly one dt-advance under its own Connectors | **FIXED** (engine-level `maxAdvances` cap, verified in isolation) |
| 10 | `decomp_` wiring, now an INCREMENTAL design (Task 17) | wiring GOOD, correctness GOOD, bit-identical to prior batch recompute (one small characterized residual, see #12) |
| 11 | Per-call file-root isolation | GOOD |
| 12 | Task 15 gate at `(n=2,j=1)` after Link 9's fix | **RESOLVED** — sticky `std::cout` precision leak in `NeqDmftSolver.h`, fixed |
| 13 | Sticky iostream format state hazard (standing, not a chain link) | Documented — check first when log values look implausibly round |

**Bottom line:** everything through Link 8 is solid and tested. Link 9's
originally-root-caused bug (border-gated advance firing causing TWO real
dt-advances per two-row segment) is FIXED via an opt-in engine-level
`maxAdvances` cap (`GroupOfOneTimeEvolutions.h`/
`NonLocalForTargetingExpression.h`), verified directly against the engine's
own trace and via `diagnosticThirdAdvanceConnectors`'s restored
C3-sensitivity. Wiring this into `fillSelfConsistentRow` improved Task 15's
gate (`(n=2,j=0)` passed) but `(n=2,j=1)` still failed — root-caused to a
SEPARATE, unrelated bug (Link 12): `NeqDmftSolver.h`'s own progress-printing
code left `std::cout` in `std::fixed`+`setprecision(1)` state, silently
truncating every subsequently-logged in-situ measurement to 1 decimal
digit. Fixed by formatting into a local stream instead. Task 15's gate now
passes 23/24 assertions (`(n=2,j=1)` matches GBEK to five decimal places).
One small, characterized, deliberately-unresolved residual remains at
`(n=2,j=0)`'s `retarded.real` (`0.000118` vs exact-zero reference, `1e-4`
margin) — see Link 12 for what was ruled out. **`NeqBathRank>0` tDMRG
results for `nT≥2` are now trustworthy to approximately the `1e-4` level
demonstrated by Task 15's gate.**

**Task 17 (same day, following session): `fillSelfConsistentRow`'s
"truncated batch recompute" (Link 10) replaced with an incremental
design**, cutting the self-consistent path's cost from
`O(nT^3*neqDmftIter)` to `O(nT^2*(1+neqDmftIter))`. Verified bit-identical
to the prior batch recompute after each of 5 implementation steps
(determinism check, hoist, promote-to-member, wire rollback, real
algorithm, cleanup); Task 15's gate gives the exact same 23/24 assertions
and the exact same `(n=2,j=0)` residual throughout. One real bug found and
fixed during implementation: a `const Column&` reference to
`scColumns_[0]` held across a `scColumns_.emplace_back()` call, which can
reallocate the vector's backing storage and dangle the reference —
manifested as an empty `RestartFilename=` Ainur parse error, caught by the
two tests that drive `fillSelfConsistentRow` through `NeqDmftSolver`'s real
corrector loop. Full design derivation in `fancy-painting-moon.md`'s
"Task #17 scope" section.
