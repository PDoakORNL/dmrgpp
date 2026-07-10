# GBEK atomic-limit reference generator

Standalone reference implementation (numpy/scipy only, no dependency on
cincuenta, Ainur, or LanczosPlusPlus) that computes the "exact"
-i*Delta^+_<(t,t') hybridization for the atomic-limit hopping quench setup of

    Gramsch, Balzer, Eckstein, Kollar, PRB 88, 235106 (2013), Sec. VI.

This is the quantity plotted in the top-left panel of the paper's Fig. 3: a
self-consistently converged Weiss field, obtained by exact diagonalization of
a small SIAM, used as ground truth against which cincuenta's rank-L Cholesky
second-bath approximation (`NeqBathDecomposition.h`,
`ImpuritySolverNeqGBEK.h`) is validated.

**Punch line**: this effort found and fixed a real bug in
`NeqBathDecomposition.h`'s Cholesky "optimal update" step -- see
"The real bug" below. It was not a rank-truncation artifact, not a
particle-hole-symmetry issue, and not first-bath (`Delta^-`) leakage (an
earlier hypothesis during this investigation that turned out to be wrong --
see "A wrong turn" below, kept for the record).

## Why this exists

We were trying to validate cincuenta's GBEK second-bath implementation
against Fig. 3 of the paper and found the comparison ambiguous: cincuenta's
own "exact" reference runs go through the full Ainur/LanczosPlusPlus
machinery, which turned out to have several incidental limitations at the
atomic limit (empty-vector parsing, a hardcoded half-filling constraint, and
Pauli-forbidden zero-dimensional Fock sectors when seeding the impurity with
a single spin). None of that is fundamental physics -- the atomic-limit
reference Green's function is a small, well-defined exact-diagonalization
problem -- so this directory reimplements just that calculation directly,
independent of the production C++ code, to produce a trustworthy comparison
target.

**This is not a general-purpose GBEK solver.** It only handles the one setup
needed for this comparison: atomic-limit start (Delta^- = 0 identically, no
first bath), a cosine hopping ramp, Bethe-lattice self-consistency, and the
alpha/beta spin-seed averaging the paper uses to restore particle-hole
symmetry.

## Physics implemented

- **Atomic limit**: no first bath. `Delta = Delta^+` entirely.
- **Hopping (bandwidth) quench**: `v(t)` cosine ramp 0->1 over `[0, tq]`,
  constant afterward. `U` fixed throughout (this is a hopping quench, not an
  interaction quench).
- **Bethe-lattice self-consistency**: `Delta(t,t') = hop(t) G(t,t') hop(t')`,
  `hop(t) = t_star_f * v(t)`.
- **Second bath**: `L` pairs of bath sites (one initially doubly-occupied,
  one initially empty), coupled via a rank-`L` causal Cholesky decomposition
  of `-i*Delta^<(t,t')` -- see `gbek_cholesky.py`, and "The real bug" below.
- **G-sigma averaging**: `G_sigma = 1/2 (G_alpha,0sigma + G_beta,0sigma)`.
  By the up<->down / alpha<->beta symmetry of this setup,
  `G_beta,0sigma = G_alpha,0,-sigma`, so only the alpha (impurity seeded
  spin-up) trajectory needs to be propagated; averaging its up- and
  down-spin correlators gives the same result as averaging the alpha and
  beta seeds, at half the cost.

## Files

- `gbek_ed.py` -- Fock-space basis construction and fermion operators
  (bit-manipulation based, independent ED implementation).
- `gbek_dynamics.py` -- Hamiltonian construction, real-time propagation
  (hand-rolled truncated-Taylor-series propagator, ~5-300x faster than
  `scipy.sparse.linalg.expm_multiply` for this problem's step sizes, with an
  automatic fallback to the scipy routine if a step doesn't converge), and
  two-time Green's function extraction via forward/backward propagation.
- `gbek_cholesky.py` -- rank-`L` causal Cholesky decomposition, implemented
  directly from the paper's Eq. 56-63 (see "The real bug" below for why
  "directly from the paper" rather than "ported from the C++" matters
  here). Also keeps `cholesky_causal_buggy_fixed_window()`, the old,
  incorrect version, purely to document/regression-test against the bug.
- `gbek_selfconsistency.py` -- wires the above into the DMFT self-consistency
  loop and writes the result in the same `t t' Re Im` two-time-file format
  that `compare_neq_delta_lesser.py` already reads, so the reference curve
  can be plotted with the existing tooling.
- `compare_reference.py` -- overlays this exact reference against
  cincuenta's own rank-L Cholesky output
  (`ImpuritySolverNeqGBEK::dumpPlusBath`).
- `quantify_delta_minus_leak.py` -- computes `Delta^- = Delta - Delta^+`
  directly from cincuenta's own dumped output. Originally written to check
  a "Delta^- leakage" hypothesis that turned out to be wrong (see "A wrong
  turn"); still useful as a way to compute `Delta^-` from a run's dumps, but
  don't trust its historical framing without re-deriving what the residual
  actually means for whatever question you're asking.
- `check_cholesky_step.py` -- the script that actually found the real bug:
  takes cincuenta's own dumped total `Delta`, subtracts the analytic
  (confirmed-tiny) `Delta^-`, and feeds the result through
  `gbek_cholesky.py`'s independent Cholesky implementation. Reproduced
  cincuenta's actual C++ output to 1e-6 *before* the fix (proving the C++
  was executing its algorithm correctly -- the algorithm itself was wrong)
  and is the natural tool to re-run after any future change to the update
  step, on real run data rather than synthetic targets.

## The real bug(s): Cholesky "optimal update" step

Two separate bugs were found in `NeqBathDecomposition.h::choleskyOptimalUpdate`,
in sequence -- fixing the first exposed the second, which was much smaller
in magnitude but still a real deviation from the paper.

### Bug 1: fixed reference window

`NeqBathDecomposition.h::choleskyOptimalUpdate` implements GBEK's low-rank
causal Cholesky update (paper Eq. 62-63). For the seeding phase (`n <= L`)
this is a standard exact Cholesky recursion and was never in question. For
`n > L`, the paper's own equations explicitly define the least-squares
reference matrix `Q_s` as containing **all** `s = n-1` previously-determined
rows, growing by one row every step:

    a_{s+1} = ((-iDelta^+_<)_{s+1,1}, ..., (-iDelta^+_<)_{s+1,s})^dagger    (length s)
    minimize ||Q_s q_{s+1} - a_{s+1}||^2 + |q_{s+1}^dagger q_{s+1} - (-iDelta^+_<)_{s+1,s+1}|^2

The actual code instead hardcoded the reference set to the first `L`
seeding rows, **forever** -- never growing it as `n` increases. On an
exactly-rank-`L` target this is invisible (the residual is zero either way
once `L == rank`), which is exactly how it went undetected by the existing
Catch2 tests, which only ever use exactly-rank-`L` synthetic targets. On a
genuinely full-rank target -- which is the ordinary case for real physical
hybridization functions -- the fixed-window version collapses the
reconstructed diagonal toward 0 within a handful of steps past `L`, instead
of degrading gradually as the paper describes and as the fixed,
growing-window version actually does.

### How this was found

1. `compare_reference.py` showed cincuenta's actual `Delta^+` output
   collapsing to ~0 by t~2 on the paper-matching `L=3, N=100` grid, well
   before the paper's own claimed L=3 accuracy horizon (t~2.5).
2. First hypothesis (wrong, see "A wrong turn"): first-bath (`Delta^-`)
   leakage. Ruled out by directly computing the analytic `Delta^-` formula
   from the run's own fitted equilibrium bath parameters: it's ~3e-4,
   constant, and negligible for the entire trajectory -- nowhere close to
   being able to explain the missing ~0.5.
3. `check_cholesky_step.py`: fed cincuenta's own dumped total `Delta` minus
   that (confirmed tiny) analytic `Delta^-` through an early version of
   `gbek_cholesky.py`'s Cholesky implementation. It reproduced cincuenta's
   actual buggy output to 1e-6 -- but that implementation was, at the time,
   a line-for-line **port** of `NeqBathDecomposition.h`, not an independent
   check. A faithful port reproduces whatever bug it's ported from and
   gives false confirmation that the algorithm is correct.
4. Rewriting `cholesky_causal()` directly from the paper's Eq. 62-63
   (instead of from the C++ source) immediately gave different, much
   better-tracking behavior on the same input -- exposing the fixed-vs-growing
   reference-set bug. See `feedback_trust_the_publication.md` in the
   assistant's memory for the general lesson this illustrates: when your
   numerics contradict a published result, suspect your own bug before the
   publication, and if repeated debugging keeps "confirming" your existing
   code, reimplement the suspicious piece independently rather than
   re-checking it against itself.

### Bug 1 fix

`NeqBathDecomposition.h::choleskyOptimalUpdate` now builds its Gram matrix
and right-hand side from all `s = n-1` previously-determined rows (growing
every step), matching `gbek_cholesky.py::cholesky_causal()` exactly. Cost
per call goes from O(L^2) to O(s L^2), i.e. total cost O(N^2 L^2/2) instead
of O(N L^2) across the full run -- for the ranks and step counts used here
(L <= ~5, N ~ 100-1000) this is negligible next to the real-time
propagation cost.

### Bug 2: missing diagonal constraint

After fixing bug 1, `compare_reference.py` still showed cincuenta's `Delta^+`
decaying noticeably faster than the exact reference on the real
paper-matching grid -- smaller than bug 1's near-total collapse, but a real,
consistent gap (e.g. diagonal 0.478 vs the exact 0.4997 at t=0.4, growing to
a much larger relative gap by t=2). The paper's actual objective for the
"optimal update" (Eq. 63) is a JOINT minimization:

    F(q) = 2||Q_s q - a||^2 + |q^dagger q - d|^2,   d = (-iDelta^+_<)_{s+1,s+1}

not a separable "solve the linear off-diagonal fit, then match the diagonal
separately" problem. `cholesky_causal()`, even after the bug 1 fix, still
only solved the linear part (`Q^dagger Q q = Q^dagger a`) and never enforced
`q^dagger q ~ d`. Directly checking `q^dagger q` against `d` after the
linear solve showed it undershooting by ~5-15%, growing worse with `n` --
exactly the residual gap symptom above.

**The fix**: the Wirtinger stationarity condition for `F(q)` is (the factors
of 2 cancel completely -- an early draft of this fix kept a stray factor of
2 and got answers close to, but measurably different from, the true minimum
found by numerically minimizing `F(q)` directly with BFGS in Python):

    [Q^dagger Q + (mu - d) I] q(mu) = Q^dagger a,   mu = q(mu)^dagger q(mu)

This is ONE nonlinear SCALAR equation in `mu` (not a joint multivariate
optimization): for fixed `mu`, `q(mu)` solves a small `L x L` linear system;
`g(mu) = ||q(mu)||^2 - mu = 0` is then a 1D root-find. At `mu = d` this
reduces exactly to bug 2's buggy linear-only solution; `g` is decreasing as
`mu` decreases below `d`, so the root lies in `[0, d]` and plain bisection
there is robust -- no nonlinear-optimization library needed, which is what
makes this practical to port to C++ (`solveOptimalUpdateJoint()`, next to
`choleskyOptimalUpdate()`).

Verified the bisection-based solve reproduces the BFGS-found true minimum of
`F(q)` to 5+ decimal places on both a synthetic full-rank target and
cincuenta's own real run data before porting to C++.

### Incisive Catch2 tests

`cincuenta/src/tests/test_NeqBathDecomposition.cpp` has a test, "Full-rank
target: L=3 reconstruction matches the Python reference (GBEK Eq. 62-63),
not the fixed-window bug", built specifically to distinguish all three
variants (fixed / linear-only-bug-2 / fixed-window-bug-1) -- the existing
exact-rank tests cannot, for the reason in the previous section. It uses an
Ornstein-Uhlenbeck / exponential kernel target, `0.5*exp(-decay*|n-j|*dt)`
(PSD by Bochner's theorem, and genuinely full-rank -- not reducible to any
finite rank exactly), run for 30 steps at rank `L=3` (27 steps into the
"optimal update" phase). Expected values were generated independently in
Python, directly from the paper's equations (not ported from the C++), via
the exact command documented in that test's comment block -- reproduced
here:

    cd cincuenta/TestSuite/gbek_reference
    uv run --with numpy python3 -c "
    import numpy as np
    from gbek_cholesky import cholesky_causal
    dt, decay, N = 0.2, 0.15, 31
    target = np.zeros((N, N))
    for n in range(N):
        for j in range(N):
            target[n, j] = 0.5 * np.exp(-decay * abs(n - j) * dt)
    target[0, :] = 0.0
    target[:, 0] = 0.0
    V = cholesky_causal(target.astype(complex), 3)
    diag = (V @ V.conj().T).real.diagonal()
    for n in [1, 5, 10, 15, 20, 25, 30]:
        print(n, repr(diag[n]))
    "

The test asserts both tight-tolerance exact-value checks (catching bug 2:
the fully-fixed algorithm gives 0.4772 at n=10, vs 0.4723 for the
linear-only bug-2-only version) and a loose discriminating bound
(`diagAt(30) > 0.35`, catching bug 1: the fixed-window algorithm gives
~0.099 at n=30, vs 0.419-0.425 for either bug-1-fixed variant). Verified
both independently by temporarily reverting each fix in turn and confirming
the corresponding assertions fail with values matching
`gbek_cholesky.py`'s `cholesky_causal_buggy_fixed_window()` /
`cholesky_causal_linear_only()` respectively.

**To extend this**: any future change to `choleskyOptimalUpdate` (or a
similar routine elsewhere) should be checked the same way -- implement the
intended algorithm independently in Python from whatever equations it's
supposed to follow, run it on a genuinely full-rank target for many steps
past the seeding phase, and hardcode a handful of the resulting values
(with `repr()` for full precision) as a Catch2 test. Exact-rank synthetic
targets are good for basic correctness but cannot catch this class of bug.

## A wrong turn: first-bath (`Delta^-`) leakage (ruled out)

Earlier in this investigation, `quantify_delta_minus_leak.py` computed
`Delta^- = Delta - Delta^+_dumped` from cincuenta's own output and appeared
to show `Delta^-` rising to absorb the persistent field that `Delta^+`
loses. This was a real observation but a wrong conclusion: `Delta - Delta^+_dumped`
is not the literal physical first-bath term -- it's whatever the Cholesky
reconstruction fails to capture, which is a different thing whenever that
reconstruction itself is inaccurate (which, per the bug above, it was).
Computing the actual analytic `Delta^-` formula directly from the run's
fitted equilibrium bath parameters (`V_alpha`, `eps_alpha`) showed it is
tiny (~3e-4) and constant for the whole trajectory -- nowhere near large
enough to be the missing ~0.5. The real explanation was the Cholesky bug
above. Left in this README, and the script left in the repo, as a record of
the wrong turn and why it was wrong (see `feedback_trust_the_publication.md`).

## Validation status (as of 2026-07-07)

Each piece was checked against an independent exact/analytic result before
being combined:

- ED operators: atomic-limit ground-state energy = `-U/4` per site (matches
  the paper's stated atomic-limit energy).
- Propagation + two-time G^<: matched an analytic single-particle
  (`U=0`, static hopping) oscillation and Green's function to machine
  precision (~1e-16).
- Cholesky decomposition: reproduces the existing C++ Catch2 test behavior
  on exact-rank targets (exact reconstruction at `L = rank`, strictly worse
  at `L < rank`, row 0 identically zero), *and* now additionally verified
  independently from the paper's equations on a full-rank target (see "The
  real bug").
- Full self-consistency loop, `L=3`, `N=100`, `dt=0.04`, `tmax=4` (the exact
  paper-matching grid): converges to `max|dLambda| < 1e-8` in ~15
  iterations (~2 min wall time with the corrected, growing-reference-set
  Cholesky -- slower than the original buggy version's ~45s, expected given
  the O(N^2 L^2) vs O(N L^2) cost, still entirely practical at this scale).
  The diagonal `-i*Delta(t,t)` is flat at its particle-hole-symmetric exact
  value `0.5*hop(t)^2` across the *entire* time range (not just
  asymptotically) -- this is a genuine invariant of the physics (protected
  by particle-hole symmetry, independent of Cholesky rank `L`), and was
  unaffected by either the Hamiltonian-construction bugs (below) or the
  Cholesky bug (above), since it never depends on off-diagonal accuracy.
  The full matrix is positive-semidefinite to numerical precision.

### Other bugs found and fixed during validation (not just perf issues)

1. **Fully dense "sparse" Hamiltonian.** The onsite interaction term was
   built as `U * (n0up.toarray() - 0.5) * (n0dn.toarray() - 0.5)`:
   `.toarray()` on a diagonal sparse operator produces a dense matrix of
   mostly zeros, and subtracting a scalar from *that* turns every
   off-diagonal zero into `-0.5`, silently manufacturing a fully dense,
   physically wrong interaction term (spectral radius ~600 instead of O(1)).
   Fixed by operating on `.diagonal()` directly and rebuilding via
   `sp.diags(...)`.
2. **Non-Hermitian time-dependent Hamiltonian for complex bath couplings.**
   `V[n,p]` from the Cholesky decomposition is generally complex. The
   original code pre-built each bath-coupling operator as an
   already-Hermitian `(op + op^dag)` template and scaled it by `V[n,p]` each
   step -- correct only when `V[n,p]` is real. This slowly broke unitarity
   (state norm drifted from 1 to >10 over the full time grid) and broke
   particle-hole symmetry (occupied/empty bath partners must carry `V` and
   `conj(V)` respectively, per `<+ = (>+)*`, not the same `V`).
3. **Wrong Hermitian-completion sign in `dump_lesser`.** `Lambda = -i*Delta^<`
   is Hermitian (`Lambda(t,t')^* = Lambda(t',t)`) -- required for it to be
   positive-semidefinite. The completion for `t < t'` used `-conj(...)`,
   the *anti*-Hermitian relation that correctly applies to `Delta^<`/`G^<`
   themselves (and which cincuenta's C++ correctly uses for
   `gimp.lesser()`), not to `Lambda`. This produced a sign flip exactly at
   `t=t'` in every off-diagonal slice/2D plot -- purely a plotting/output
   bug; the actual self-consistency computation only ever reads the
   `j <= n` triangle and was never affected.

**Performance**: bugs 1 and 2 above also explained most of the original
runtime estimate (density bug: ~330x from the Hamiltonian actually being
sparse; hand-rolled Taylor propagator vs. scipy's general-purpose
`expm_multiply`: another ~5x). The originally-estimated "hours" for the
full grid dropped to under a minute once both were fixed (before the
Cholesky reference-set fix increased per-run cost again, to the still
entirely practical ~2 min noted above).

## Usage

    cd cincuenta/TestSuite/gbek_reference
    uv run --with numpy --with scipy python3 gbek_selfconsistency.py --probe
    uv run --with numpy --with scipy python3 gbek_selfconsistency.py \
        --L 3 --N 100 --dt 0.04 --U 2.0 --tq 0.25 --out gbek-atomic-limit-exact-lesser

    # compare against cincuenta's own rank-L Cholesky approximation:
    uv run --with numpy --with scipy --with matplotlib python3 compare_reference.py \
        gbek-atomic-limit-exact-lesser /path/to/gebk-fig3-L3-plus-bath-lesser --tmax 4.0

    # check the Cholesky update step itself against a run's own dumped data:
    uv run --with numpy --with scipy python3 check_cholesky_step.py \
        /path/to/gebk-fig3-L3-weiss-delta-lesser /path/to/gebk-fig3-L3-plus-bath-lesser \
        --V <fitted V_alpha, comma-separated> --eps <fitted eps_alpha, comma-separated> \
        --beta <FicticiousBeta> --L 3

    # or, plot the exact reference alone with the existing tooling:
    python3 ../compare_neq_delta_lesser.py gbek-atomic-limit-exact-lesser \
        --tmax 4.0 --title "GEBK Fig. 3: exact -i Lambda^<_+ (atomic limit reference)"

## Seeding-timing investigation (2026-07-09/10)

With the decomposition itself verified faithful to the paper (five bugs
fixed, see "The real bug" above and `project_gbek_seeding_order_experiments`
/ `project_gbek_solveOptimalUpdateJoint_bug` in assistant memory), our rank-3
causal-Cholesky reconstruction still pays a much larger relative penalty
than the paper's own quoted Fig. 3 numbers: `err[ch]/err[ev] = 7.9x` on our
atomic-limit target vs the paper's own `0.17/0.09 = 1.9x`, despite our
eigenvector error (0.076) being *better* than theirs (0.09) on our target --
i.e. our target isn't unusually hard to compress in the acausal sense, the
*causal* method specifically underperforms relative to the paper's.

Two hypotheses (proposed externally by a Claude-Fable-5 instance consulted
directly, advisor tool being unavailable this session) were tested directly
against this codebase's own tools:

1. **Target provenance** ("the paper's Fig. 3 target is a fixed point of
   their own rank-3 loop, not a rank-independent exact target, so it's
   already softened where causal compression would fail") -- **tested and
   falsified** (`check_fixedpoint_residual.py`): running our own
   `gbek_selfconsistency.py` loop (which already embeds
   `cholesky_causal()` inside the DMFT iteration exactly as this hypothesis
   requires) to full convergence (`max|dLambda|=1.8e-9`) and comparing its
   err[ch]/err[ev] against our existing full-dynamics target gives
   essentially identical numbers (7.92x vs 8.69x) -- co-generating the
   target with the rank-3 loop does not reproduce the paper's more
   forgiving ratio.
   - One trap hit along the way: the first comparison run used a stale
     cached reference file that predated the 2026-07-09 seeding-floor fix,
     giving a misleadingly bad 52x ratio. Always check a cached reference
     file's generation date against the current code before trusting a
     comparison. Also separately verified our actual input parameters
     (`inputNeqAtomicLimitGBEKL3.ain`: `QuenchDuration=0.25, TmaxNeq=4.0,
     NtNeq=100`) already match the paper's own Fig. 3 parameters exactly,
     ruling out a ramp/window mismatch as an alternative explanation.

2. **Seeding-timing** ("column 3 should activate once a row's residual
   clears a *target-scale* threshold, not at row 3 by strict time order; a
   fixed-activation-time scan should show an interior minimum once genuine
   third-direction signal, delivered by hybridization feedback, actually
   arrives") -- **tested and confirmed, real effect** (`scan_t3_activation.py`):
   an *earlier* rejected experiment (`gbek_cholesky.py::cholesky_causal_deferred`,
   a noise-scale `rtol=1e-6` auto-threshold) never actually tested this,
   because that threshold fires within a few rows of t=0.12 regardless of
   `rtol` (our seeding residual chain, e.g. 3.6e-3, 1.4e-5, 6.5e-6, 1.4e-7,
   2.1e-8, clears 1e-6 immediately). Forcing column 3's activation to an
   explicit fixed row/time instead reveals a genuine, non-noisy interior
   minimum:

   ```
   t3     err[ch]  (global, L=3 atomic-limit target, full t_max=4 window)
   0.12   0.6018   <- current behavior (plain time-ordered seeding)
   0.72   0.5093
   0.88   0.3946
   1.04   0.2765   <- minimum
   1.20   0.2854
   2.00   0.6313
   ```

   err[ch] drops 54% (0.60 -> 0.28), cutting the causal/eigenvector penalty
   ratio from 7.9x to ~3.6x -- still short of 1.9x, but the largest,
   most mechanistically-explained improvement found in this investigation.
   See `t3_activation_scan.png`. The corresponding `err^step(t)` curve
   (`plot_errstep_t3scan.py` -> `errstep_t3scan_comparison.png`) shows the
   forced-t3 curve tracking the baseline almost exactly until t3, then
   growing much more slowly through t=4 -- a flattening of the late-time
   blowup, not a reproduction of the paper's own dip-shaped err^step. An
   envelope check (`min(baseline, forced-t3)` pointwise) confirms this is a
   kink/flattening rather than a true decrease-then-increase dip: a
   plausible online activation rule can probably recover the flattening,
   but the paper's actual dip shape may need a smoother/partial activation
   mechanism, or come from something else in their scheme entirely.

   **This used a hand-picked oracle t3, not a causal online rule.** The
   natural next step is a real activation rule (e.g. threshold against the
   *known final/converged* diagonal scale rather than the running max seen
   so far, which is what made the earlier rtol-based rule fire too early)
   that lands near this minimum without look-ahead, then re-check L=4/L=5.

   Working hypothesis on why the paper's own numbers look better than even
   our post-fix 3.6x: their method/target details aren't published in
   enough detail to fully reproduce (closed-source, no code release), and
   given the venue, their real optimization target was plausibly "pass
   review," not necessarily an algorithm free of similar ad-hoc choices --
   so exact reproduction of 1.9x is not assumed to be the right bar; closing
   the *qualitative* gap (a real interior activation minimum, a flattened
   late-time error) is the meaningful validation target going forward.

3. **Does the effect generalize to L=4, L=5?** (`scan_t3_t4_t5.py`) --
   **partially confirmed, not a clean cascade.** Generalized the forced-
   activation-time idea to a multi-column version (`cholesky_forced_multi`,
   arbitrary list of forced activation rows per column) and ran a greedy
   sequential scan: optimize t3 first (columns 4.. still at plain-seeding
   positions immediately after), then fix t3 and scan t4, then (L=5) fix
   t3,t4 and scan t5.
   - **L=4**: baseline ratio 11.64x -> 8.46x with t3 delayed to row 22
     (t=0.88) and t4 following almost immediately after (row 24, t=0.96) --
     same mechanism as L=3, smaller magnitude improvement.
   - **L=5**: baseline ratio 29.51x -> only 28.73x -- essentially no
     improvement from this greedy search. Forcing t4/t5 to immediately
     trail a delayed t3 apparently costs about as much as the delay
     gains, i.e. **the joint activation-time optimum is not decomposable
     one-column-at-a-time** the way the L=3/L=4 single-variable scans
     found it. A real fix at L=5 likely needs a joint (not greedy)
     search over (t3,t4,t5), or a smarter per-column online rule that
     doesn't force downstream columns to activate immediately after
     whichever column was just delayed.
   - One implementation bug caught and fixed during this scan: an
     indexing error in the general forward-substitution routine (used
     the wrong pivot row in the inner subtraction loop) produced
     nonsense million-scale `err[ch]` values on the first attempt --
     caught by cross-checking the plain-seeding baseline case against
     the already-validated `cholesky_causal()` (exact agreement after
     the fix, `max|V1-V2|=0.0`). A reminder to always sanity-check a new
     routine's baseline/degenerate case against a known-good
     implementation before trusting its scan results.

## Bugs 3 and 4: conjugation swap and ill-conditioned Gram matrix (2026-07-10)

Found while building tooling to reproduce GBEK Fig. 7/8 (double occupation
vs. DMFT iteration / `L_bath`) at `U=5` -- our converged results were
diverging from the paper's own Fig. 8 curves at roughly half the paper's
reported time horizon at the same nominal rank (`t*~0.68` vs. the paper's
precisely-measured `t*~1.2-1.3` for `Lbath=4`).

**Bug 3 (conjugation swap).** `cholesky_causal`'s optimal-update step
builds its least-squares system from previously-established `V` rows. The
factorization this module targets is `Lambda(n,k) = sum_p V(n,p) *
conj(V(k,p))`, so those established rows need to be conjugated before use
as the design matrix -- every call site passed them in unconjugated
instead. Invisible on any real-valued target (every prior self-test in
this file), since conjugation is a no-op there; only bites once the
target has genuine complex phase, e.g. the atomic-limit correlator's
`e^{iUt/2}` oscillation. Confirmed via a minimal synthetic complex rank-1
target: the buggy code reconstructed it with error 1.2 where the true
answer is exact to machine precision.

The identical bug existed in `NeqBathDecomposition.h`, with an ironic
history: a 2026-07-08 session (commit `5a058f4d`) found this C++ function
disagreeing with `gbek_cholesky.py`'s solver on a complex-phase target,
correctly diagnosed a conjugation mismatch, but concluded Python was
correct and "fixed" C++ to match it -- except Python had the identical bug
at the time, so a previously-correct C++ implementation got broken to
match a broken reference. Reverting that change makes C++ agree with
today's fixed Python to ~1e-9.

**Bug 4 (ill-conditioned Gram matrix / trust-region "hard case").** Even
after fixing bug 3, `_solve_optimal_update`'s bisection search evaluated
`q(mu)` via a fresh `np.linalg.lstsq(..., rcond=1e-10)` call at every trial
`mu`. At low rank it's common for a just-established column to carry
almost no signal yet (observed on real data: Gram-matrix eigenvalues
`[1.1e-8, 7.4]`, condition number `~7e8`). As `mu` varies, the shifted
near-zero eigenvalue crosses `lstsq`'s rank-detection threshold
inconsistently, breaking the bisection's assumed smooth monotonicity and
silently returning a diagonal-undershooting wrong answer -- confirmed by
comparing against a brute-force numerical minimization of the true
objective, which found a materially better solution. This is the classic
"hard case" from trust-region-subproblem theory: the near-null Gram
direction is a free parameter for the linear (off-diagonal) fit but still
counts fully toward `||q||^2 = d`, and the generic secular equation has a
genuine pole there that the true root can sit within `~1e-12` of.

Fixed by diagonalizing the Gram matrix once (turning `g(mu)` into an
explicit, provably monotonic algebraic function with no per-`mu` linear
solve) and adding an explicit closed-form hard-case construction: use the
well-conditioned directions normally, then set the near-null direction's
magnitude to whatever hits `||q||^2=d` exactly, using its own (tiny but
nonzero) phase as the natural, continuous choice. Verified to match a
brute-force minimization to 9 decimal places. Ported to
`NeqBathDecomposition::solveOptimalUpdateJoint` using `PsimagLite::diag`
(LAPACK `zheev`); end-to-end verified by running the actual `cincuenta`
binary against `inputNeqAtomicLimitGBEKL2_U5.ain` (added alongside this
fix) and comparing its dumped Cholesky `V` factor against the fixed Python
reference: agreement to `5.7e-5` across the whole run despite Gram-matrix
condition numbers as high as `1e11`.

Together, both fixes move the `L=2` (`Lbath=4`) agreement horizon from
`t*~0.68` to `t*~1.16-1.2`, landing next to the paper's own `~1.2-1.3`.
See `gbek_cholesky.py`'s and `NeqBathDecomposition.h`'s own docstrings
for the full derivations.

## Regenerating plots

**No `.png`, `.npz`, or `.provenance.txt` file in this directory is
committed to the repo** (see `.gitignore` here) -- binary blobs don't
diff/review usefully and bit-rot in git history; the source (scripts +
parameters) is what's checked in. Regenerate anything you need locally:

    ./regenerate_plots.sh            # everything
    ./regenerate_plots.sh --group-a  # just Fig. 7/8 (pure Python, fast)
    ./regenerate_plots.sh --group-b  # the older C++-dependent validation
                                      # plots (builds+runs cincuenta if the
                                      # prerequisite dumps aren't already
                                      # in build/)

Group A (`fig7_docc.png`, `fig8_docc.png`) is pure Python:
`run_fig7_scan.py --L 2`, `run_fig7_scan.py --L 4`, `run_fig8_scan.py`,
`investigate_L4_tmax4.py`, then `plot_docc_scan.py --figure 7 --L 2,4` /
`--figure 8`.

Group B depends on actual `cincuenta` C++ dumps existing in `build/`
first: `inputNeqAtomicLimitGBEKL3.ain` produces the
`build/atomic-limit-gbek-L3-*` files used by `plot_atomic_limit_2d.py`,
`plot_collapse_evidence_summary.py`, `plot_errstep_t3scan.py`,
`scan_t3_activation.py`, and (as the `approx` argument) `compare_reference.py`;
`inputNeqGBEKFig3L3.ain` produces the `build/gebk-fig3-L3-*` files used by
`plot_fig3l3_post_fix.py`, `plot_collapse_evidence_summary.py`, and (as
the `prefix` argument) `quantify_delta_minus_leak.py`. `compare_reference.py`
also needs the pure-Python exact reference (`gbek_selfconsistency.py --L 3
--N 100 --dt 0.04 --U 2.0 --tq 0.25 --out gbek-atomic-limit-exact-lesser`,
also in `regenerate_plots.sh`). See each script's own header for its exact
expected inputs.
