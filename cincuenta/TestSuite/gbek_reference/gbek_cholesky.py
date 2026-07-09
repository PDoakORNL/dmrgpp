"""
Rank-L causal Cholesky decomposition of -i*Lambda^<(t_n,t_j), independently
implemented directly from Gramsch, Balzer, Eckstein, Kollar, PRB 88, 235106
(2013), Eqs. 56-63 -- the reference this file's implementation is checked
against is the paper's equations, not cincuenta/src/NeqBathDecomposition.h.

V[n, p] for n=0..nT, p=0..L-1.  V[0, :] = 0 exactly (n=0 skipped, matching
the "n==0 is skipped" fix in NeqBathDecomposition.h -- that part of the
algorithm matches the paper and was never in question).

## Two bugs this file caught, in sequence

An earlier version of this file was a line-for-line port of
NeqBathDecomposition.h::choleskyOptimalUpdate, and "confirmed" that C++ code
was correct because it reproduced the same output -- but a faithful port
shares any bug in what it's ported from; it does not independently verify
the algorithm. See feedback_trust_the_publication.md in the assistant's
memory for the general lesson. Both bugs below were found by instead
implementing directly from the paper's own equations and comparing.

**Bug 1 -- fixed reference window.** The paper's "optimal update" step
(n > L, Eq. 62) explicitly builds its least-squares reference matrix Q_s
from ALL s = n-1 previously-determined rows, growing by one row every step.
The original C++ (and the first port) instead hardcoded the reference set
to the first L seeding rows forever. See
cholesky_causal_buggy_fixed_window()'s docstring.

**Bug 2 -- missing diagonal constraint.** The paper's objective (Eq. 63) is

    F(q) = 2||Q q - a||^2 + |q^dagger q - d|^2

a JOINT minimization over the off-diagonal fit (Qq ~ a) AND the diagonal
target d = (-iDelta^+_<)_{s+1,s+1} -- not a separable linear-least-squares
problem. After fixing bug 1, cholesky_causal() still only solved the linear
part (Qq ~ a) and never enforced q^dagger q ~ d, which systematically
undershoots the diagonal (verified: q^dagger q ends up ~5-15% below the
true target d on realistic inputs, growing worse with n). See
cholesky_causal_linear_only()'s docstring and _solve_optimal_update() below
for the fix.

NeqBathDecomposition.h::choleskyOptimalUpdate has been fixed to match this
file's cholesky_causal() (both bugs).
"""
import numpy as np


def _solve_optimal_update(Q, a, d, L):
    """
    Solve the paper's Eq. 63 joint minimization exactly:

        minimize_q  F(q) = 2||Q q - a||^2 + (q^dagger q - d)^2

    Setting the Wirtinger gradient d F / d q* to zero gives the stationarity
    condition (the factors of 2 cancel completely -- easy to get wrong; an
    earlier draft of this fix kept a stray factor of 2 and got answers that
    were close to, but measurably different from, the true minimum found by
    numerically minimizing F(q) directly with BFGS):

        [Q^dagger Q + (mu - d) I] q(mu) = Q^dagger a,   mu = q(mu)^dagger q(mu)

    This is ONE nonlinear scalar equation in mu (not a joint multivariate
    optimization): for any fixed mu, q(mu) solves a small L x L linear
    system; substituting back into mu = q(mu)^dagger q(mu) gives a scalar
    equation g(mu) = ||q(mu)||^2 - mu = 0, solved here by bisection (after
    bracketing) rather than a general-purpose multivariate optimizer -- this
    is what makes the fix straightforward to port to C++ (see
    NeqBathDecomposition.h::choleskyOptimalUpdate): no nonlinear-optimization
    library needed, just a linear solve inside a 1D root-find.

    At mu = d, this reduces exactly to the plain linear-least-squares
    solution (bug 2's buggy behavior): [Q^dagger Q] q = Q^dagger a.
    ||q(mu)||^2 is monotonically decreasing in mu (A(mu) = Q^dagger Q +
    (mu-d)I only gains positive-definiteness as mu grows, shrinking its
    inverse), so g(mu) = ||q(mu)||^2 - mu is strictly decreasing too. Since
    g(0) >= 0 generically and g(mu) -> -infinity as mu -> infinity, exactly
    one root exists in (0, infinity) -- but NOT always in [0, d]. When
    g(d) >= 0 (the diagonal-unconstrained fit's norm already meets/exceeds
    d), monotonicity forces the true root ABOVE d: this function used to
    just return q(d) in that case (bracketing only ever searched [0, d]),
    silently accepting an overshoot -- q(d) is not even a stationary point
    of Eq. 63 when its norm exceeds d. Found 2026-07-09 on the atomic-limit
    self-consistency target, whose off-diagonal magnitude does not decay
    (it's a pure oscillating phase), unlike the decaying-ridge target this
    bisection was originally validated against; that's why the gap went
    unnoticed. Fixed by bracketing upward in that regime too, PROVIDED
    d >= 0 and Q^dagger Q is well-conditioned, both true for a converged,
    physically PSD target. During self-consistency ITERATION
    (gbek_selfconsistency.py), intermediate/not-yet-converged targets need
    not satisfy either: d can go slightly negative from numerical noise, or
    -- the actual failure hit in practice -- Q^dagger Q can be NEAR-singular
    this early past the seeding phase (observed: two of three columns of Q
    built entirely from ~1e-9 seed values, at the very first "optimal
    update" step). np.linalg.solve has NO protection against a near-singular
    (as opposed to exactly singular) system: it doesn't raise, it silently
    returns enormous garbage (observed: |q| ~ 1e6 on real self-consistency
    data). That then makes every downstream propagation step catastrophically
    slow (their cost scales with the operator norm) rather than crashing, so
    it looks like a hang rather than an error. Fixed by using lstsq (SVD
    based) as the primary solve throughout instead of solve() -- lstsq
    naturally regularizes a near-singular system by dropping tiny singular
    values instead of dividing by them, which is exactly the failure mode
    here. d is also clamped to be non-negative for the same
    intermediate-estimate reason.
    """
    QtQ = Q.conj().T @ Q
    Qta = Q.conj().T @ a
    I = np.eye(L)

    def q_of(mu):
        A = QtQ + (mu - d) * I
        return np.linalg.lstsq(A, Qta, rcond=1e-10)[0]

    d = max(d, 0.0)
    q_linear = q_of(d)
    n2_linear = np.vdot(q_linear, q_linear).real
    if not np.isfinite(n2_linear):
        return np.zeros(L, dtype=complex)  # QtQ itself is degenerate; give up safely

    def g(mu):
        q  = q_of(mu)
        n2 = np.vdot(q, q).real
        return n2

    n2_at_d = g(d)
    if not np.isfinite(n2_at_d):
        return q_linear

    if n2_at_d >= d:
        # True root lies ABOVE d (see docstring): bracket upward instead of
        # silently accepting q(d), whose norm n2_at_d overshoots d.
        mu_lo, mu_hi = d, d
        step  = max(d, 1.0) * 0.5
        found = False
        for _ in range(60):
            mu_hi += step
            n2 = g(mu_hi)
            if np.isfinite(n2) and n2 - mu_hi <= 0:
                found = True
                break
            step *= 1.5
        if not found:
            return q_linear
    else:
        mu_lo, mu_hi = 0.0, d
        step  = max(d, 1.0) * 0.5
        cur   = d
        found = False
        for _ in range(60):
            cur = max(cur - step, 0.0)
            n2  = g(cur)
            if np.isfinite(n2) and (n2 - cur >= 0 or cur <= 0.0):
                mu_lo, found = cur, True
                break
            step *= 1.5
        if not found:
            return q_linear

    for _ in range(100):
        mu_mid = 0.5 * (mu_lo + mu_hi)
        n2 = g(mu_mid)
        if np.isfinite(n2) and n2 - mu_mid >= 0:
            mu_lo = mu_mid
        else:
            mu_hi = mu_mid
        if mu_hi - mu_lo < 1e-14 * max(mu_hi, 1.0):
            break

    q_final = q_of(mu_lo)
    n2_final = np.vdot(q_final, q_final).real
    # Final sanity check rather than guarding every intermediate step (which
    # is fragile near a near-singular shift): a genuine root satisfies
    # n2_final ~ mu_lo; a numerically-blown-up spurious solve does not, and
    # in any case should never wildly exceed the natural scale of the data.
    bound = max(d, n2_linear, n2_at_d, 1.0) * 10
    if (not np.isfinite(n2_final) or n2_final > bound
            or abs(n2_final - mu_lo) > 0.05 * max(mu_lo, 1.0)):
        return q_linear
    return q_final


def cholesky_causal(lam, L, floor_rtol=1e-10):
    """
    lam: (nT+1, nT+1) complex array, lam[n, j] = -i*Lambda^<(t_n, t_j)
      (only need j <= n filled; must be Hermitian PSD conceptually).
    L: rank.
    floor_rtol: seeding-phase floor (see below), as a fraction of the
      largest diagonal magnitude seen so far. Default 1e-10.
    Returns V of shape (nT+1, L).

    Seeding-phase zero-clamp floor (added 2026-07-09, see
    project_gbek_solveOptimalUpdateJoint_bug and the seeding-fragility
    follow-on): the exact-Cholesky seeding residual `d` for column L-1 (the
    last column, seeded from the most information-starved early row) is
    often the tiny difference of O(1)-scale terms -- confirmed via an
    independent 50-digit mpmath replay that this is NOT a linear-algebra
    precision problem, it's that the true residual, given the available
    precision of Lambda itself, genuinely sits at the edge of what's
    resolvable. Previously this clamped hard to exactly 0 whenever it
    rounded negative (`sqrt(max(d,0))`), which is numerically catastrophic:
    a column seeded at EXACTLY zero makes that entire column's row/column
    in every future Gram matrix (QtQ) identically zero, which then
    structurally FORCES every subsequent optimal-update solve to keep that
    component at exactly zero forever too (`(mu-d)*q[p] = 0` with mu != d
    has no other solution) -- i.e. one unlucky sign flip on a
    borderline-noise residual permanently kills an entire bath direction
    for the rest of the run (confirmed: this is exactly what happened in
    cincuenta's own C++ output at L=5, while an independent Python
    recomputation of the identical algorithm on a very slightly different
    (but equally valid) floating-point path did not collapse and got a
    good rank-5 reconstruction instead).
    A small positive floor (rather than 0) fixes this: the seeded column
    stays numerically ALIVE (nonzero Gram row/column), so later
    optimal-update rows -- which see genuine, non-noise data at those
    rows -- remain free to grow that direction into something meaningful
    instead of being locked at zero by construction. The floor's absolute
    scale barely matters (only this one row's value is fixed; everything
    after is re-solved from real data) -- it just needs to be nonzero and
    small enough not to distort anything, hence tying it to a tiny
    fraction of the run's own diagonal scale rather than a fixed constant.
    A rank-revealing/deferred-establishment alternative (skip seeding a
    column at low-residual rows entirely) was also tried and found NOT to
    help (see git history around 2026-07-09): deferring forces every
    skipped row to assert EXACTLY zero coupling in the not-yet-established
    column, which can be just as wrong as a noisy nonzero value, and
    measured no better overall. This floor is the narrower, verified fix.
    """
    nT = lam.shape[0] - 1
    V = np.zeros((nT + 1, L), dtype=complex)
    max_diag_seen = 0.0

    def lam_at(n, j):
        # lam = -i*Delta^< is Hermitian (lam(t,t')^* = lam(t',t)); the causal
        # algorithm only ever calls this with j <= n in practice, so this
        # branch is never exercised during real self-consistency use, but
        # must still use the correct (Hermitian, no extra minus sign)
        # completion for standalone use (e.g. this file's own smoke test).
        if j <= n:
            return lam[n, j]
        return np.conj(lam[j, n])

    for n in range(1, nT + 1):
        max_diag_seen = max(max_diag_seen, abs(lam_at(n, n).real))
        if n <= L:
            # Standard (exact) Cholesky seeding phase, Eq. 56-58: column p is
            # seeded at pivot row p+1 (row 0 is degenerate, per the header
            # docstring). Never in question -- identical across all variants
            # in this file.
            col = n - 1
            for p in range(col):
                val = lam_at(n, p + 1)
                for k in range(p):
                    val -= V[n, k] * np.conj(V[p + 1, k])
                pivot = V[p + 1, p]
                V[n, p] = 0 if abs(pivot) < 1e-14 else val / np.conj(pivot)
            d = lam_at(n, n).real
            for k in range(col):
                d -= abs(V[n, k]) ** 2
            floor = floor_rtol * max(max_diag_seen, 1e-300)
            V[n, col] = np.sqrt(max(d, floor))
        else:
            # Optimal update, Eq. 62-63: Q_s is the s x L matrix of ALL
            # previously-determined rows 1..s (s = n-1), growing every step
            # (bug 1's fix). Solve the JOINT off-diagonal + diagonal
            # objective exactly via _solve_optimal_update (bug 2's fix).
            a = np.array([lam_at(n, k) for k in range(1, n)])  # length s, ALL prior rows
            Q = V[1:n, :]                                       # (s, L), ALL prior rows
            d = lam_at(n, n).real
            V[n, :] = _solve_optimal_update(Q, a, d, L)
    return V


def cholesky_causal_pivoted_seed(lam, L, verbose=False):
    """
    Same algorithm as cholesky_causal() for rows n>L (Eq. 62-63 optimal
    update, unchanged), but replaces the SEEDING phase (rows 1..L) with a
    pivoted (rank-revealing) Cholesky decomposition of the fully-known L x L
    Hermitian block, instead of forcing strict time-order seeding (column p
    always from row p+1).

    Motivation (2026-07-09): the paper's own text describing this step
    ("During the first L time steps we have enough free parameters V^ch_np
    for AN exact matrix decomposition of -iLambda^<_+ and can thus rely on
    the Cholesky decomposition to fill up the matrix Q^L") is looser than
    "use the same Eq. 57 recursion in strict time order": an L x L Hermitian
    PSD matrix admits a whole family of exact square-root factorizations
    related by unitary rotation, of which strict-time-order Cholesky is
    only one member. Once L full timesteps have been collected, the ENTIRE
    L x L block is already known (not a one-row-at-a-time stream), so
    nothing about causality prevents choosing a BETTER-CONDITIONED member
    of that family: standard pivoted/rank-revealing Cholesky, which at each
    step seeds the next column from whichever of the L already-known rows
    currently has the largest remaining diagonal residual, instead of
    blindly taking them in time order.

    This is NOT the same as the earlier-tried cholesky_causal_deferred()
    (which was found NOT to help): that function deferred a column's
    establishment to a FUTURE, not-yet-seen row, which forces every row
    before the eventual pivot to assert an unverifiable "exactly zero"
    coupling in that direction. Pivoting WITHIN the fully-known L x L block
    never asserts anything about unseen data -- it only reorders among rows
    already fully in hand, so there is no analogous "confidently wrong
    zero" failure mode.

    Returns (V, pivot_rows): V as in cholesky_causal(); pivot_rows[k] is the
    actual time-row index (1-based, matching the rest of this module's
    convention) chosen as the pivot for column k.
    """
    nT = lam.shape[0] - 1
    V = np.zeros((nT + 1, L), dtype=complex)

    def lam_at(n, j):
        if j <= n:
            return lam[n, j]
        return np.conj(lam[j, n])

    # Rows 1..L (row 0 is degenerate/skipped, matching cholesky_causal()).
    rows = list(range(1, L + 1))  # actual time-row indices available to seed from
    M = np.array([[lam_at(rows[i], rows[j]) for j in range(L)] for i in range(L)],
                 dtype=complex)

    Vblock = np.zeros((L, L), dtype=complex)  # Vblock[original block row i, column k]
    diag_resid = np.array([M[i, i].real for i in range(L)])
    remaining = list(range(L))
    pivot_rows = []

    for k in range(L):
        best = max(remaining, key=lambda i: diag_resid[i])
        remaining.remove(best)
        pivot_val = max(diag_resid[best], 0.0)
        Vblock[best, k] = np.sqrt(pivot_val)
        if verbose:
            print(f"  column {k}: pivot = row n={rows[best]} (t={rows[best]}), "
                  f"residual={diag_resid[best]:.3e}")
        if abs(Vblock[best, k]) >= 1e-14:
            for i in remaining:
                val = M[i, best]
                for kk in range(k):
                    val -= Vblock[i, kk] * np.conj(Vblock[best, kk])
                Vblock[i, k] = val / np.conj(Vblock[best, k])
            for i in remaining:
                diag_resid[i] -= abs(Vblock[i, k]) ** 2
        pivot_rows.append(rows[best])

    for i in range(L):
        V[rows[i], :] = Vblock[i, :]

    for n in range(L + 1, nT + 1):
        a = np.array([lam_at(n, m) for m in range(1, n)])
        Q = V[1:n, :]
        d = lam_at(n, n).real
        V[n, :] = _solve_optimal_update(Q, a, d, L)

    return V, pivot_rows


def cholesky_causal_deferred(lam, L, rtol=1e-6, verbose=False):
    """
    Causally rank-revealing variant of cholesky_causal(): identical algebra,
    but does NOT force column p to be seeded at row p+1 unconditionally.
    Instead it tracks how many columns k (<=L) have been genuinely
    established so far, and at each row n (while k<L) computes the EXACT
    forward-substituted residual left over after explaining this row's
    diagonal with the k already-established columns. A new column is only
    opened (pivot_row.append(n); k += 1) if that residual exceeds
    `rtol * (largest diagonal magnitude seen so far)`; otherwise the row is
    left with zeros in the not-yet-established columns and processing moves
    to row n+1, deferring establishment of column k to whenever a row
    actually offers enough new signal. Once k reaches L, subsequent rows
    fall back to the ordinary Eq. 63 optimal-update solve exactly as in
    cholesky_causal(), using ALL prior rows (their already-established
    entries, retroactively-zero not-yet-established ones included) as the
    reference set -- causality is preserved throughout: only already-seen
    rows are ever used to decide anything.

    Motivation (found 2026-07-09, see project_gbek_solveOptimalUpdateJoint_bug
    and the follow-on rank-scan/seeding-fragility investigation): the plain
    causal seeding recursion computes each successive column's residual as
    a small difference of O(1)-scale terms (observed: 3.6e-3, 1.4e-5,
    6.5e-6, 1.4e-7, 2.1e-8 for L=5's first 5 rows on the true-atomic-limit
    target) -- i.e. the LOCAL rank of the seeding window is genuinely much
    smaller than L, so the later "columns" are seeded almost entirely from
    numerical noise (confirmed via an independent 50-digit mpmath replay of
    the identical subtraction chain: it reproduces the double-precision
    residuals almost exactly, meaning the arithmetic itself is not the
    problem -- the true mathematical residual, given the available
    precision of the input Lambda, really is that small and that
    unreliable). Baking a noise-dominated direction in as a permanent basis
    column then measurably hurts short-time accuracy (see err^step results)
    and, in the worst case (L=5, a residual of 2e-8), can flip sign between
    two otherwise-equivalent implementations and collapse an entire
    column to exactly zero for the rest of the run. Deferring establishment
    until real signal appears avoids both failure modes without needing
    higher-precision arithmetic (verified not to be the bottleneck here).

    rtol: relative threshold (fraction of the largest diagonal magnitude
      seen so far) below which a would-be new column's residual is treated
      as unresolvable noise rather than real structure. Default 1e-6 sits
      comfortably between this target's 3rd (6.5e-6) and 4th (1.4e-7)
      column residuals -- i.e. it accepts the first 3 genuine directions
      and defers the rest. Tune per-target if the natural residual gap
      isn't this clean.
    verbose: print each column-establishment decision (row, residual,
      accepted/deferred) for diagnostic use.

    Returns (V, pivot_row): V as in cholesky_causal(); pivot_row is the
    list of row indices at which each established column was actually
    seeded (len(pivot_row) <= L -- if less, the target's genuine local rank
    never reached L anywhere in [0, t_max], which is a valid finding, not
    an error).
    """
    nT = lam.shape[0] - 1
    V = np.zeros((nT + 1, L), dtype=complex)

    def lam_at(n, j):
        if j <= n:
            return lam[n, j]
        return np.conj(lam[j, n])

    pivot_row = []
    k = 0
    max_diag_seen = 0.0

    for n in range(1, nT + 1):
        d_nn = lam_at(n, n).real
        max_diag_seen = max(max_diag_seen, abs(d_nn))

        if k < L:
            for p in range(k):
                val = lam_at(n, pivot_row[p])
                for q in range(p):
                    val -= V[n, q] * np.conj(V[pivot_row[p], q])
                piv = V[pivot_row[p], p]
                V[n, p] = 0 if abs(piv) < 1e-14 else val / np.conj(piv)

            r = d_nn - sum(abs(V[n, p]) ** 2 for p in range(k))
            threshold = rtol * max(max_diag_seen, 1e-300)
            if r > threshold:
                V[n, k] = np.sqrt(max(r, 0.0))
                pivot_row.append(n)
                if verbose:
                    print(f"  establish column {k} at row n={n}: residual={r:.3e} "
                          f"(threshold={threshold:.3e})")
                k += 1
            elif verbose:
                print(f"  defer column {k} at row n={n}: residual={r:.3e} "
                      f"<= threshold={threshold:.3e}")
        else:
            a = np.array([lam_at(n, m) for m in range(1, n)])
            Q = V[1:n, :]
            d = d_nn
            V[n, :] = _solve_optimal_update(Q, a, d, L)

    if verbose and k < L:
        print(f"  WARNING: only {k}/{L} columns ever established -- target's "
              f"genuine local rank never reached L in [0, t_max={nT}*dt]")

    return V, pivot_row


def cholesky_causal_linear_only(lam, L):
    """
    Growing reference window (bug 1 fixed), but still missing the diagonal
    constraint (bug 2 NOT fixed): solves only the linear part of Eq. 63,
    Q q ~ a, ignoring the q^dagger q ~ d term entirely. This was
    cholesky_causal()'s behavior between the two fixes -- kept here purely
    for regression/contrast (see check_cholesky_step.py and the "second,
    smaller bug" discussion in README.md). Systematically undershoots the
    true diagonal by ~5-15% on realistic full-rank targets, growing worse
    with n, but nowhere near bug 1's near-total collapse.
    """
    nT = lam.shape[0] - 1
    V = np.zeros((nT + 1, L), dtype=complex)

    def lam_at(n, j):
        if j <= n:
            return lam[n, j]
        return np.conj(lam[j, n])

    for n in range(1, nT + 1):
        if n <= L:
            col = n - 1
            for p in range(col):
                val = lam_at(n, p + 1)
                for k in range(p):
                    val -= V[n, k] * np.conj(V[p + 1, k])
                pivot = V[p + 1, p]
                V[n, p] = 0 if abs(pivot) < 1e-14 else val / np.conj(pivot)
            d = lam_at(n, n).real
            for k in range(col):
                d -= abs(V[n, k]) ** 2
            V[n, col] = np.sqrt(max(d, 0.0))
        else:
            a = np.array([lam_at(n, k) for k in range(1, n)])
            Q = V[1:n, :]
            G = Q.conj().T @ Q
            b = Q.conj().T @ a
            V[n, :] = np.linalg.lstsq(G, b, rcond=None)[0]
    return V


def cholesky_causal_buggy_fixed_window(lam, L):
    """
    The original bug (bug 1), preserved here for regression/contrast purposes
    only -- do not use this for anything other than demonstrating why the
    growing-reference-window fix in cholesky_causal() matters.

    Hardcodes its least-squares reference set to the first L seeding rows
    forever, instead of growing it to all s = n-1 previously-determined
    rows. This was NeqBathDecomposition.h::choleskyOptimalUpdate's actual
    behavior before any of this module's fixes.
    """
    nT = lam.shape[0] - 1
    V = np.zeros((nT + 1, L), dtype=complex)

    def lam_at(n, j):
        if j <= n:
            return lam[n, j]
        return np.conj(lam[j, n])

    for n in range(1, nT + 1):
        if n <= L:
            col = n - 1
            for p in range(col):
                val = lam_at(n, p + 1)
                for k in range(p):
                    val -= V[n, k] * np.conj(V[p + 1, k])
                pivot = V[p + 1, p]
                V[n, p] = 0 if abs(pivot) < 1e-14 else val / np.conj(pivot)
            d = lam_at(n, n).real
            for k in range(col):
                d -= abs(V[n, k]) ** 2
            V[n, col] = np.sqrt(max(d, 0.0))
        else:
            a = np.array([lam_at(n, k + 1) for k in range(L)])
            G = np.zeros((L, L), dtype=complex)
            for p in range(L):
                for pp in range(L):
                    G[p, pp] = sum(V[k, p] * np.conj(V[k, pp]) for k in range(1, L + 1))
            b = np.zeros(L, dtype=complex)
            for p in range(L):
                b[p] = sum(V[k, p] * a[k - 1] for k in range(1, L + 1))
            q = np.linalg.lstsq(G, b, rcond=None)[0]
            V[n, :] = q
    return V


def reconstruct(V):
    return V @ V.conj().T


if __name__ == "__main__":
    # Mirror the C++ Catch2 tests: synthetic rank-r PSD target matrices,
    # check reconstruction error decreases with L and L=r reconstructs exactly.
    def make_psd_target(nT, rank):
        N = nT + 1
        v = np.zeros((rank, N))
        for r in range(rank):
            decay = 0.3 + 0.2 * r
            freq = 0.5 * (r + 1)
            for n in range(1, N):
                v[r, n] = np.exp(-decay * n) * np.sin(freq * n)
        mat = np.zeros((N, N))
        for r in range(rank):
            mat += np.outer(v[r], v[r])
        return mat

    nT = 10
    target2 = make_psd_target(nT, 2)

    V1 = cholesky_causal(target2.astype(complex), 1)
    err1 = np.max(np.abs(reconstruct(V1).real - target2))
    V2 = cholesky_causal(target2.astype(complex), 2)
    err2 = np.max(np.abs(reconstruct(V2).real - target2))
    print("L=1 error on rank-2 target:", err1)
    print("L=2 error on rank-2 target:", err2)
    assert err2 < 1e-8, "L=2 should reconstruct a rank-2 matrix exactly"
    assert err2 < err1 and err1 > 1e-4, "L=2 must be strictly better than L=1"
    assert np.max(np.abs(V1[0, :])) < 1e-15, "row 0 must be exactly zero"
    print("OK: Cholesky decomposition matches expected C++ behavior")

    # Demonstrate both bugs on a genuinely full-rank target: an exact-rank
    # target cannot distinguish any of these algorithms (residual is zero
    # either way once L == rank), which is exactly why both bugs went
    # undetected by the pre-existing tests.
    print()
    dt, decay = 0.2, 0.15
    N = 31
    full_rank_target = np.zeros((N, N))
    for n in range(N):
        for j in range(N):
            full_rank_target[n, j] = 0.5 * np.exp(-decay * abs(n - j) * dt)
    full_rank_target[0, :] = 0.0
    full_rank_target[:, 0] = 0.0

    V_fixed = cholesky_causal(full_rank_target.astype(complex), 3)
    V_linear = cholesky_causal_linear_only(full_rank_target.astype(complex), 3)
    V_buggy = cholesky_causal_buggy_fixed_window(full_rank_target.astype(complex), 3)
    diag_fixed = reconstruct(V_fixed).real.diagonal()
    diag_linear = reconstruct(V_linear).real.diagonal()
    diag_buggy = reconstruct(V_buggy).real.diagonal()
    print(f"{'n':>3} {'t':>5} {'target':>8} {'fixed':>8} {'linear-only':>12} {'buggy':>8}")
    for n in [1, 5, 10, 15, 20, 25, 30]:
        print(f"{n:3d} {n*dt:5.2f} {full_rank_target[n,n]:8.5f} "
              f"{diag_fixed[n]:8.5f} {diag_linear[n]:12.5f} {diag_buggy[n]:8.5f}")
