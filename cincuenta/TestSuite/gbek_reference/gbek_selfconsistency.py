"""
Full GBEK atomic-limit self-consistency loop (Gramsch, Balzer, Eckstein,
Kollar, PRB 88, 235106 (2013), Sec. VI), producing the "exact" reference
-i*Lambda^+_<(t,t') = -i*Delta^+_<(t,t') used for Fig. 3's top-left panel.

Physics:
  - Atomic-limit start: Delta^- = 0 identically (no first bath), so
    Delta = Delta^+ entirely.
  - Hopping (bandwidth) quench: v(t) cosine ramp 0 -> 1 over [0, tq], then
    constant.  U fixed throughout.
  - Bethe-lattice self-consistency: Delta(t,t') = hop(t) G(t,t') hop(t'),
    hop(t) = t_star_f * v(t).
  - Second bath: L pairs of bath sites (one initially doubly-occupied, one
    initially empty), coupled via a rank-L causal Cholesky decomposition of
    -i*Delta^<(t,t') -- see gbek_cholesky.py, implemented directly from GBEK
    Eq. 56-63 (NeqBathDecomposition.h::choleskyOptimalUpdate was fixed to
    match it after this independent implementation exposed a bug there;
    see that file's module docstring and the README's "The real bug").
  - Gsigma(t,t') = 1/2 (Galpha_0sigma + Gbeta_0sigma).  By the up<->down /
    alpha<->beta symmetry of this setup, Gbeta_0sigma = Galpha_0,-sigma, so
    only the alpha (impurity seeded spin-up) trajectory needs to be computed;
    averaging its up- and down-spin correlators gives the same result as
    averaging the alpha and beta seeds.
"""
import numpy as np

from provenance import write_provenance
import scipy.sparse as sp

from gbek_ed import build_basis, c_operator, number_operator, double_occupation_operator
from gbek_dynamics import hop_operator, propagate, compute_g_lesser, compute_diag_observable
from gbek_cholesky import cholesky_causal, eigenvector_decompose


def v_ramp(t, tq):
    if tq <= 0:
        return 1.0 if t > 0 else 0.0
    if t <= 0:
        return 0.0
    if t >= tq:
        return 1.0
    w0 = np.pi / tq
    return 0.5 * (1 - np.cos(w0 * t))


def build_templates(nsites, basis, index, U, L):
    """
    U_term: onsite interaction (particle-hole symmetric, mu_imp = -U/2 built in).
    occ_templates[p] / empty_templates[p]: RAW (not pre-symmetrized) sparse
      matrix c0^dag c_site (summed over spin) for the occupied/empty partner
      of bath pair p, so that
        H(t) = U_term
             + sum_p [V[n,p]      * occ_templates[p]   + conj(V[n,p])      * occ_templates[p]^dag
                      + conj(V[n,p]) * empty_templates[p] + V[n,p]          * empty_templates[p]^dag]
      The occupied bath site couples via -i*Delta^<_+ = sum_p V_p(t)V_p(t')^*
      (Eq. 52a); the empty site couples via i*Delta^>_+, and Delta^<_+ =
      (Delta^>_+)^* (Eq. 49/particle-hole relation) means the empty site's
      hopping is conj(V[n,p]), NOT the same V[n,p] as the occupied site.
      Using the same V for both (as an earlier version of this script did)
      is invisible for real V but breaks particle-hole symmetry -- and hence
      the exactly-0.5 diagonal -- as soon as V goes complex.
    Site layout: site 0 = impurity; for pair p (0-indexed), occ site = 1+2p,
    empty site = 2+2p.
    """
    n0up_diag = number_operator(basis, 0, 0).diagonal()
    n0dn_diag = number_operator(basis, 0, 1).diagonal()
    U_term = sp.diags(U * (n0up_diag - 0.5) * (n0dn_diag - 0.5)).tocsr().astype(complex)

    occ_templates, empty_templates = [], []
    for p in range(L):
        occ_site = 1 + 2 * p
        empty_site = 2 + 2 * p
        occ_term = sp.csr_matrix((len(basis), len(basis)), dtype=complex)
        empty_term = sp.csr_matrix((len(basis), len(basis)), dtype=complex)
        for spin in (0, 1):
            occ_term = occ_term + hop_operator(basis, index, 0, occ_site, spin)
            empty_term = empty_term + hop_operator(basis, index, 0, empty_site, spin)
        occ_templates.append(occ_term)
        empty_templates.append(empty_term)
    return U_term, occ_templates, empty_templates


def make_h_seq(U_term, occ_templates, empty_templates, V):
    """V: (N, L) complex hopping values for steps n=0..N-1 (time t_n, held for dt)."""
    N, L = V.shape
    H_seq = []
    for n in range(N):
        H = U_term.copy()
        for p in range(L):
            if V[n, p] != 0:
                v, vc = V[n, p], np.conj(V[n, p])
                H = H + v * occ_templates[p] + vc * occ_templates[p].conj().T
                H = H + vc * empty_templates[p] + v * empty_templates[p].conj().T
        H_seq.append(H)
    return H_seq


def initial_state(nsites, L, basis, index):
    up_word = 1  # impurity spin-up occupied
    dn_word = 0
    for p in range(L):
        occ_site = 1 + 2 * p
        up_word |= (1 << occ_site)
        dn_word |= (1 << occ_site)
        # empty_site stays 0 in both words
    psi0 = np.zeros(len(basis), dtype=complex)
    psi0[index[(up_word, dn_word)]] = 1.0
    return psi0


def run_self_consistency(L, N, dt, U, tstar_f, tq, n_iterations, tol=1e-6,
                          mode="cholesky", verbose=True, record_history=False):
    """
    mode: "cholesky" (default, causal rank-L decomposition of Lambda every
      iteration) or "eigenvector" (non-causal top-L spectral decomposition
      of the FULL Lambda every iteration) -- see gbek_cholesky.py's
      cholesky_causal()/eigenvector_decompose() for the algorithms. Used to
      reproduce GBEK Fig. 7/8's top-vs-bottom-panel comparison.
    record_history: if True, additionally computes and records the
      double occupation d_it(t) = <n0up(t) n0dn(t)> after EVERY iteration
      (not just the converged one) in the returned docc_history list, needed
      for Fig. 7's per-iteration curves. d(t) is computed directly from the
      alpha trajectory phi_states with no beta-sector averaging: an
      empirical check (check_alpha_beta_docc_symmetry.py) confirms
      d_alpha(t) == d_beta(t) to numerical precision, consistent with this
      module's existing G_avg alpha/beta-symmetry argument (docstring above).
    """
    nsites = 1 + 2 * L
    ts = np.arange(N + 1) * dt
    v = np.array([v_ramp(t, tq) for t in ts])
    hop_t = tstar_f * v

    nup_ext, ndown_ext = 1 + L, L
    basis_N, index_N = build_basis(nsites, nup_ext, ndown_ext)
    basis_up, index_up = build_basis(nsites, nup_ext - 1, ndown_ext)
    basis_dn, index_dn = build_basis(nsites, nup_ext, ndown_ext - 1)
    if verbose:
        print(f"L={L} nsites={nsites} dims: N={len(basis_N)} "
              f"up-1={len(basis_up)} dn-1={len(basis_dn)}")

    Uterm_N, occ_N, empty_N = build_templates(nsites, basis_N, index_N, U, L)
    Uterm_up, occ_up, empty_up = build_templates(nsites, basis_up, index_up, U, L)
    Uterm_dn, occ_dn, empty_dn = build_templates(nsites, basis_dn, index_dn, U, L)

    c0up = c_operator(nsites, basis_N, index_N, basis_up, index_up, site=0, spin=0)
    c0dn = c_operator(nsites, basis_N, index_N, basis_dn, index_dn, site=0, spin=1)

    psi0 = initial_state(nsites, L, basis_N, index_N)
    docc_op = double_occupation_operator(basis_N, site=0) if record_history else None

    Lambda = np.zeros((N + 1, N + 1), dtype=complex)  # -i*Delta^<(t_n,t_j), j<=n
    history = []
    docc_history = [] if record_history else None
    for it in range(n_iterations):
        if mode == "cholesky":
            V = cholesky_causal(Lambda, L)  # (N+1, L); V[0,:] = 0 by construction
        elif mode == "eigenvector":
            V = eigenvector_decompose(Lambda, L)
        else:
            raise ValueError(f"unknown mode {mode!r}")
        H_seq_N = make_h_seq(Uterm_N, occ_N, empty_N, V[:N, :])
        H_seq_up = make_h_seq(Uterm_up, occ_up, empty_up, V[:N, :])
        H_seq_dn = make_h_seq(Uterm_dn, occ_dn, empty_dn, V[:N, :])

        phi_states = propagate(psi0, H_seq_N, dt)
        G_up = compute_g_lesser(phi_states, c0up, H_seq_up, dt)
        G_dn = compute_g_lesser(phi_states, c0dn, H_seq_dn, dt)
        G_avg = 0.5 * (G_up + G_dn)

        if record_history:
            docc_history.append(compute_diag_observable(phi_states, docc_op))

        Lambda_new = np.zeros_like(Lambda)
        for n in range(N + 1):
            for j in range(n + 1):
                Lambda_new[n, j] = hop_t[n] * hop_t[j] * (-1j) * G_avg[n, j]

        diff = np.max(np.abs(Lambda_new - Lambda))
        history.append(diff)
        if verbose:
            diag05 = Lambda_new[N, N].real
            print(f"  iter {it}: max|dLambda|={diff:.3e}  "
                  f"Lambda(tmax,tmax).real={diag05:.5f} (expect ~{0.5*v[N]**2:.5f})")
        Lambda = Lambda_new
        if diff < tol:
            break

    return Lambda, V, hop_t, ts, history, docc_history


def dump_lesser(path, Lambda, dt):
    """
    Write -i*Delta^<(t,t') in the same "t t' Re Im" format used by
    KadanoffBaym::dump (see cincuenta/TestSuite/compare_neq_delta_lesser.py),
    so the reference curve can be plotted with the existing tooling.

    Lambda = -i*Delta^< is HERMITIAN (Lambda(t,t')^* = Lambda(t',t)), not
    anti-Hermitian: the "-i" flips the parity of the underlying anti-Hermitian
    relation Delta^<(t,t')^* = -Delta^<(t',t) (the same relation cincuenta's
    C++ correctly uses for gimp.lesser() itself, i.e. for Delta^</G^< without
    the "-i" factor). Using "-conj" here instead of "+conj" produced a sign
    flip exactly at t=t' in the reconstructed upper triangle.
    """
    N = Lambda.shape[0] - 1
    with open(path, "w") as fh:
        for n in range(N + 1):
            for j in range(N + 1):
                val = Lambda[n, j] if j <= n else np.conj(Lambda[j, n])
                fh.write(f"{n*dt:.6f} {j*dt:.6f} {val.real:.8e} {val.imag:.8e}\n")


def _probe():
    # Compute-discipline probe: small N, few iterations, to check convergence
    # rate and basic sanity checks before scaling up.  Try a couple of L
    # values to see whether the diagonal deviation from 0.5 and the PSD
    # violation shrink with increasing rank (truncation effect) or persist
    # (implementation bug).
    dt, U, tstar_f, tq = 0.04, 2.0, 1.0, 0.25
    N = 15

    for L in (1, 2, 3):
        print(f"\n=== L={L} ===")
        Lambda, V, hop_t, ts, history, _ = run_self_consistency(
            L, N, dt, U, tstar_f, tq, n_iterations=10, tol=1e-8, verbose=False)
        print("Convergence history:", [f"{h:.2e}" for h in history])

        diag = np.array([Lambda[n, n].real for n in range(N + 1)])
        expect_diag = 0.5 * (tstar_f * np.array([v_ramp(t, tq) for t in ts])) ** 2
        print("Diagonal -i*Lambda(t,t) vs 0.5*hop(t)^2:")
        for n in [0, 6, 7, 10, N]:
            print(f"  t={ts[n]:.3f}  Lambda_diag={diag[n]:.5f}  expect={expect_diag[n]:.5f}")

        full = np.zeros_like(Lambda)
        for n in range(N + 1):
            for j in range(N + 1):
                full[n, j] = Lambda[n, j] if j <= n else np.conj(Lambda[j, n])
        eigs = np.linalg.eigvalsh(full)
        print("min eigenvalue of -i*Lambda (should be >= 0):", eigs.min())


def main():
    import argparse
    import time

    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--probe", action="store_true",
                    help="run the small L=1,2,3 / N=15 sanity probe instead of a full grid")
    p.add_argument("--L", type=int, default=3, help="second-bath rank (Lbath = 2L)")
    p.add_argument("--N", type=int, default=100, help="number of real-time steps")
    p.add_argument("--dt", type=float, default=0.04)
    p.add_argument("--U", type=float, default=2.0)
    p.add_argument("--tstar-f", type=float, default=1.0, help="final Bethe-lattice hopping")
    p.add_argument("--tq", type=float, default=0.25, help="cosine ramp duration")
    p.add_argument("--iterations", type=int, default=30,
                    help="max DMFT iterations; the L=3,N=100 paper grid needs ~15 to reach 1e-8")
    p.add_argument("--tol", type=float, default=1e-8)
    p.add_argument("--out", default="gbek-atomic-limit-exact-lesser",
                    help="output file, same 't t\\' Re Im' format as compare_neq_delta_lesser.py")
    args = p.parse_args()

    if args.probe:
        _probe()
        return

    t0 = time.time()
    Lambda, V, hop_t, ts, history, _ = run_self_consistency(
        args.L, args.N, args.dt, args.U, args.tstar_f, args.tq,
        n_iterations=args.iterations, tol=args.tol, verbose=True)
    print(f"\nWall time: {time.time()-t0:.1f} s")
    print("Convergence history:", [f"{h:.2e}" for h in history])

    dump_lesser(args.out, Lambda, args.dt)
    print(f"Wrote {args.out}")
    prov = write_provenance(
        args.out,
        notes=f"L={args.L} N={args.N} dt={args.dt} U={args.U} tq={args.tq} "
              f"iterations={args.iterations} tol={args.tol} "
              f"final_diff={history[-1]:.3e}")
    print(f"Wrote {prov}")


if __name__ == "__main__":
    main()
