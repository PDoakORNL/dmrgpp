#!/usr/bin/env python3
"""Compare G_imp(t,t') against the Weiss field G_0(t,t') for a no-quench (U_i=U_f=0) run.

Physical interpretation
-----------------------
When U_i = U_f = 0 (no interaction quench) the self-energy Sigma_imp = 0 at all times.
The Dyson equation G_imp^{-1} = G_0^{-1} - Sigma_imp then gives G_imp = G_0 exactly.
Failure of this condition indicates a bug in the neq impurity solver or the Volterra
integro-differential solver for the Weiss field G_0.

Files read
----------
  green-retarded      : G_imp^R(t,t')  written by NeqDmftSolver::dumpGreenFunctions
  weiss-green-retarded: G_0^R(t,t')    same call, from NeqLatticeGf::g0()
  green-lesser        : G_imp^<(t,t')
  weiss-green-lesser  : G_0^<(t,t')

Format: "t  t'  Re(G)  Im(G)"  one row per (t,t') pair.

Usage
-----
    python3 compare_gimp_weiss.py \\
        --retarded-gimp  green-retarded \\
        --retarded-weiss weiss-green-retarded \\
        --lesser-gimp    green-lesser \\
        --lesser-weiss   weiss-green-lesser \\
        [--tolerance TOL]

Exit code 0 on success, 1 if max deviation exceeds --tolerance.
"""

import argparse
import math
import sys


def read_kb_file(path):
    """Return dict {(t, t'): (re, im)} from a KadanoffBaym dump file."""
    data = {}
    with open(path) as fh:
        for line in fh:
            cols = line.split()
            if len(cols) < 4:
                continue
            try:
                t  = round(float(cols[0]), 8)
                tp = round(float(cols[1]), 8)
                re = float(cols[2])
                im = float(cols[3])
            except ValueError:
                continue
            data[(t, tp)] = (re, im)
    if not data:
        print(f"ERROR: no data read from {path}", file=sys.stderr)
        sys.exit(1)
    return data


def compare_component(label, gimp_file, weiss_file, tolerance):
    """Read both files, compute max deviation, print table. Returns (max_err, passed)."""
    gimp  = read_kb_file(gimp_file)
    weiss = read_kb_file(weiss_file)

    common = sorted(set(gimp) & set(weiss))
    if not common:
        print(f"ERROR: no common (t,t') points between {gimp_file} and {weiss_file}",
              file=sys.stderr)
        sys.exit(1)

    errors = []
    print(f"\n{label}: G_imp vs G_0 (Weiss field) — should agree exactly at U=0")
    print(f"  gimp : {gimp_file}")
    print(f"  weiss: {weiss_file}")
    print(f"\n{'t':>10}  {'t_prime':>10}  {'Re Gimp':>10}  {'Im Gimp':>10}  "
          f"{'Re G0':>10}  {'Im G0':>10}  {'|err|':>8}")
    print("-" * 74)

    for key in common:
        t, tp = key
        rg, ig = gimp[key]
        rw, iw = weiss[key]
        err = math.sqrt((rg - rw) ** 2 + (ig - iw) ** 2)
        errors.append(err)
        print(f"{t:>10.4f}  {tp:>10.4f}  {rg:>10.5f}  {ig:>10.5f}  "
              f"{rw:>10.5f}  {iw:>10.5f}  {err:>8.5f}")

    max_err  = max(errors)
    mean_err = sum(errors) / len(errors)
    rms_err  = math.sqrt(sum(e ** 2 for e in errors) / len(errors))
    print(f"\n  Points compared : {len(errors)}")
    print(f"  Max  |G_imp - G_0| = {max_err:.6f}")
    print(f"  Mean |G_imp - G_0| = {mean_err:.6f}")
    print(f"  RMS  |G_imp - G_0| = {rms_err:.6f}")

    passed = (tolerance is None) or (max_err <= tolerance)
    if not passed:
        print(f"\n  FAIL: max deviation {max_err:.6f} exceeds tolerance {tolerance}",
              file=sys.stderr)
    else:
        if tolerance is not None:
            print(f"\n  PASS: max deviation {max_err:.6f} <= tolerance {tolerance}")

    return max_err, passed


def main():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--retarded-gimp",  default="green-retarded",
                   metavar="FILE", help="G_imp^R file (green-retarded)")
    p.add_argument("--retarded-weiss", default="weiss-green-retarded",
                   metavar="FILE", help="G_0^R file  (weiss-green-retarded)")
    p.add_argument("--lesser-gimp",    default="green-lesser",
                   metavar="FILE", help="G_imp^< file (green-lesser)")
    p.add_argument("--lesser-weiss",   default="weiss-green-lesser",
                   metavar="FILE", help="G_0^< file  (weiss-green-lesser)")
    p.add_argument("--tolerance", default=None, type=float, metavar="TOL",
                   help="Fail if max|G_imp - G_0| exceeds TOL for either component")
    args = p.parse_args()

    print("U=0 no-quench check: G_imp must equal Weiss field G_0 when Sigma_imp=0\n")

    _, ok_r = compare_component("G^R (retarded)",
                                args.retarded_gimp, args.retarded_weiss, args.tolerance)
    _, ok_l = compare_component("G^< (lesser)",
                                args.lesser_gimp, args.lesser_weiss, args.tolerance)

    if not (ok_r and ok_l):
        sys.exit(1)


if __name__ == "__main__":
    main()
