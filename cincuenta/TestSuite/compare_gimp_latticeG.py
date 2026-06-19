#!/usr/bin/env python3
"""Compare G_imp(iw_n) against G_cluster_0(iw_n) for all Matsubara frequencies.

Physical interpretation
-----------------------
For U=0 (no electron-electron interactions) the self-energy Sigma=0 exactly.
The DMFT self-consistency condition then requires G_imp = G_cluster_0 (the local
lattice Green's function) up to bath-discretization error.  Failure of this
condition indicates a bug in the impurity solver or the DMFT loop.

File format
-----------
Both files are written by cincuenta's DmftSolver and have the same layout:
    omega  Re(G)  Im(G)
one line per Matsubara frequency (negative frequencies listed first, then positive).

Usage
-----
    python3 compare_gimp_latticeG.py gimp_exactdiag.txt latticeG_exactdiag.txt \\
        [--tolerance TOL]

Exit code 0 on success, 1 if max deviation exceeds --tolerance.
"""

import argparse
import math
import sys


def read_matsubara_file(path):
    """Return list of (omega, re, im) tuples from a 'omega Re Im' file."""
    data = []
    with open(path) as fh:
        for line in fh:
            cols = line.split()
            if len(cols) < 3:
                continue
            try:
                data.append((float(cols[0]), float(cols[1]), float(cols[2])))
            except ValueError:
                continue
    if not data:
        print(f"ERROR: no data read from {path}", file=sys.stderr)
        sys.exit(1)
    return data


def main():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("gimp",     metavar="GIMP_FILE",
                   help="gimp_<solver>.txt written by DmftSolver after Matsubara solve")
    p.add_argument("latticeG", metavar="LATTICEG_FILE",
                   help="latticeG_<solver>.txt written alongside gimp (G_cluster_0)")
    p.add_argument("--tolerance", default=None, type=float, metavar="TOL",
                   help="Fail if max|G_imp - G_cluster_0| exceeds TOL (default: print only)")
    args = p.parse_args()

    gimp_data = read_matsubara_file(args.gimp)
    latt_data = read_matsubara_file(args.latticeG)

    if len(gimp_data) != len(latt_data):
        print(
            f"ERROR: length mismatch: gimp has {len(gimp_data)} lines, "
            f"latticeG has {len(latt_data)} lines",
            file=sys.stderr,
        )
        sys.exit(1)

    errors = []
    for (og, rg, ig), (ol, rl, il) in zip(gimp_data, latt_data):
        if abs(og - ol) > 1e-10:
            print(
                f"WARNING: omega mismatch at gimp={og:.6g} vs latticeG={ol:.6g}",
                file=sys.stderr,
            )
        err = math.sqrt((rg - rl) ** 2 + (ig - il) ** 2)
        errors.append(err)

    max_err  = max(errors)
    mean_err = sum(errors) / len(errors)
    rms_err  = math.sqrt(sum(e ** 2 for e in errors) / len(errors))

    # Only print positive-frequency summary (last half of the list)
    n_pos = len(gimp_data) // 2
    pos_gimp = gimp_data[n_pos:]
    pos_latt = latt_data[n_pos:]
    pos_errs = errors[n_pos:]

    header = (
        "\nU=0 check: G_imp(iw_n) vs G_cluster_0(iw_n) — should agree up to bath fitting\n"
        f"  Files: {args.gimp}  vs  {args.latticeG}\n"
    )
    print(header)
    print(f"{'omega':>12}  {'Re Gimp':>12}  {'Im Gimp':>12}  "
          f"{'Re Glatt':>12}  {'Im Glatt':>12}  {'|err|':>10}")
    print("-" * 76)
    for (og, rg, ig), (_, rl, il), e in zip(pos_gimp, pos_latt, pos_errs):
        print(f"{og:>12.6f}  {rg:>12.6f}  {ig:>12.6f}  "
              f"{rl:>12.6f}  {il:>12.6f}  {e:>10.6f}")

    print(f"\nAll {len(errors)} Matsubara frequencies:")
    print(f"  Max  |G_imp - G_cluster_0| = {max_err:.6f}")
    print(f"  Mean |G_imp - G_cluster_0| = {mean_err:.6f}")
    print(f"  RMS  |G_imp - G_cluster_0| = {rms_err:.6f}")

    if args.tolerance is not None:
        if max_err > args.tolerance:
            print(
                f"\nFAIL: max |G_imp - G_cluster_0| = {max_err:.6f} "
                f"exceeds tolerance {args.tolerance}",
                file=sys.stderr,
            )
            sys.exit(1)
        print(f"\nPASS: max deviation {max_err:.6f} <= tolerance {args.tolerance}")


if __name__ == "__main__":
    main()
