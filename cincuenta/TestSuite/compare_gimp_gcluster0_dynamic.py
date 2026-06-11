#!/usr/bin/env python3
"""
Compare G_imp(iw_n) from TargetingDynamic continued-fraction data with the
exact G_cluster0(iw_n) for a non-interacting AIM (U=0).

Usage:
    python3 compare_gimp_gcluster0_dynamic.py <hd5_type0> <hd5_type1> \\
        --nbath <N> --v <V> --eps_bath <e1 e2 ...> [--beta <beta>] \\
        --tol <tol>

Reads CF data from /Def/FinalPsi/CF{Eigs,Intensities,Energy,Isign,Weight}.
Applies physical weight alpha0=(1-n_imp), alpha1=n_imp from non-interacting SP.
Checks max |G_imp(iw_n) - G_cluster0(iw_n)| / |G_cluster0(iw_n)| <= tol.

Exit code: 0 if within tolerance, 1 otherwise.
"""
import argparse
import sys
import h5py
import numpy as np


def sp_hamiltonian(eps_imp, eps_bath, V_bath):
    n = len(eps_bath)
    H = np.zeros((n + 1, n + 1))
    H[0, 0] = eps_imp
    for i, (V, e) in enumerate(zip(V_bath, eps_bath), 1):
        H[i, i] = e
        H[0, i] = H[i, 0] = V
    return H


def noninteracting_nimp(eps_imp, eps_bath, V_bath, n_occ):
    """Impurity occupancy from non-interacting SP spectrum (per spin)."""
    H = sp_hamiltonian(eps_imp, eps_bath, V_bath)
    ev, evec = np.linalg.eigh(H)
    return float(sum(evec[0, l] ** 2 for l in range(n_occ)))


def gcluster0_matsubara(wn, eps_imp, eps_bath, V_bath):
    """Exact local Green's function at site 0 for complex z = i*wn."""
    H = sp_hamiltonian(eps_imp, eps_bath, V_bath)
    ev, evec = np.linalg.eigh(H)
    z = complex(0, wn)
    return sum(evec[0, l] ** 2 / (z - ev[l]) for l in range(len(ev)))


def eval_cf(z, eigs, intensities, Eg, weight, isign):
    """Evaluate continued-fraction G(z) = weight * sum_l I_l / (z - isign*(Eg-eig_l))."""
    return weight * sum(
        i / (z - isign * (Eg - e)) for e, i in zip(eigs, intensities)
    )


def read_cf(hd5_path):
    with h5py.File(hd5_path, "r") as f:
        grp = f["Def"]["FinalPsi"]
        eigs = grp["CFEigs"][:]
        ints = grp["CFIntensities"][:]
        Eg = float(grp["CFEnergy"][0])
        weight = float(grp["CFWeight"][0])
        isign = int(grp["CFIsign"][0])
    return eigs, ints, Eg, weight, isign


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("hd5_type0", help="HDF5 output from DynamicDmrgType=0 run")
    p.add_argument("hd5_type1", help="HDF5 output from DynamicDmrgType=1 run")
    p.add_argument("--nbath", type=int, required=True, help="Number of bath sites")
    p.add_argument("--v", type=float, default=0.5, help="Hybridization V (uniform)")
    p.add_argument(
        "--eps_bath",
        type=float,
        nargs="+",
        required=True,
        help="Bath on-site energies",
    )
    p.add_argument("--eps_imp", type=float, default=0.0, help="Impurity on-site energy")
    p.add_argument(
        "--nelec_up",
        type=int,
        required=True,
        help="Number of up-spin electrons in the GS (TargetElectronsUp)",
    )
    p.add_argument("--beta", type=float, default=20.0, help="Inverse temperature")
    p.add_argument("--nmat", type=int, default=20, help="Number of Matsubara frequencies")
    p.add_argument("--tol", type=float, default=1e-6, help="Max relative error tolerance")
    args = p.parse_args()

    eps_bath = args.eps_bath
    V_bath = [args.v] * args.nbath

    n_imp_up = noninteracting_nimp(args.eps_imp, eps_bath, V_bath, args.nelec_up)
    alpha0 = 1.0 - n_imp_up   # weight for Type 0 (addition, c†)
    alpha1 = n_imp_up          # weight for Type 1 (removal, c)

    eigs0, ints0, Eg0, w0, isign0 = read_cf(args.hd5_type0)
    eigs1, ints1, Eg1, w1, isign1 = read_cf(args.hd5_type1)

    print(f"n_imp(up)={n_imp_up:.6f}, alpha0={alpha0:.6f}, alpha1={alpha1:.6f}")
    print(f"nPoles(type0)={len(eigs0)}, nPoles(type1)={len(eigs1)}")
    print(f"sumI(type0)={sum(ints0):.8f}, sumI(type1)={sum(ints1):.8f}")

    errors = []
    print(f"\n{'n':>4}  {'wn':>8}  {'Gimp.imag':>14}  {'Gc.imag':>14}  {'rel_err':>10}")
    for n in range(args.nmat):
        wn = (2 * n + 1) * np.pi / args.beta
        z = complex(0, wn)
        g0 = eval_cf(z, eigs0, ints0, Eg0, w0, isign0)
        g1 = eval_cf(z, eigs1, ints1, Eg1, w1, isign1)
        gimp = alpha0 * g0 + alpha1 * g1
        gc = gcluster0_matsubara(wn, args.eps_imp, eps_bath, V_bath)
        rel_err = abs(gimp - gc) / abs(gc)
        errors.append(rel_err)
        print(f"{n:4d}  {wn:8.4f}  {gimp.imag:14.8f}  {gc.imag:14.8f}  {rel_err:10.3e}")

    max_err = max(errors)
    print(f"\nmax rel_err = {max_err:.3e}  (tol={args.tol:.0e})")

    if max_err > args.tol:
        print(f"FAIL: max relative error {max_err:.3e} exceeds tolerance {args.tol:.0e}")
        sys.exit(1)
    else:
        print(f"PASS: max relative error {max_err:.3e} <= tolerance {args.tol:.0e}")
        sys.exit(0)


if __name__ == "__main__":
    main()
