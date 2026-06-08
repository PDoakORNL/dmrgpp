#!/usr/bin/env python3
"""Compare G_imp(iw_n) against G_cluster0(iw_n) for a U=0 DMFT run.

Both files have the format written by DmftSolver::logDebug / writeGcluster0ForDebug:
    omega  Re(G)  Im(G)
one line per Matsubara frequency.

At U=0 the self-energy is identically zero, so the impurity solver must reproduce
the non-interacting cluster propagator exactly (to solver precision).  Any
discrepancy larger than the tolerance indicates a bug in how bath parameters are
passed to or used by the impurity solver -- notably the particleholesymmetric
bath-parameter layout introduced in the Anderson_symPH fix (PR #195).

Usage:
    python3 compare_gimp_gcluster0.py gimp_exactdiag.txt gcluster0_exactdiag.txt \\
        [--tolerance 1e-4]

Exit 0 if max|G_imp - G_cluster0| <= tolerance, else 1.
"""

import argparse
import sys


def read_freq_file(path):
    """Return dict {omega: (Re, Im)} from a 'omega Re Im' file."""
    data = {}
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split()
            if len(parts) < 3:
                continue
            omega, re, im = float(parts[0]), float(parts[1]), float(parts[2])
            data[omega] = (re, im)
    return data


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('gimp',      help='gimp_<solver>.txt')
    parser.add_argument('gcluster0', help='gcluster0_<solver>.txt')
    parser.add_argument('--tolerance', type=float, default=1e-4,
                        help='max allowed |G_imp - G_cluster0| (default 1e-4)')
    args = parser.parse_args()

    gimp      = read_freq_file(args.gimp)
    gcluster0 = read_freq_file(args.gcluster0)

    common = sorted(set(gimp) & set(gcluster0))
    if not common:
        print(f"ERROR: no common frequencies between {args.gimp} and {args.gcluster0}",
              file=sys.stderr)
        sys.exit(1)

    max_diff = 0.0
    worst_omega = None
    for omega in common:
        re_diff = abs(gimp[omega][0] - gcluster0[omega][0])
        im_diff = abs(gimp[omega][1] - gcluster0[omega][1])
        diff = max(re_diff, im_diff)
        if diff > max_diff:
            max_diff = diff
            worst_omega = omega

    print(f"max|G_imp - G_cluster0| = {max_diff:.3e}  (tolerance {args.tolerance:.3e})"
          f"  worst omega = {worst_omega}")

    if max_diff > args.tolerance:
        print(f"FAIL: difference {max_diff:.3e} exceeds tolerance {args.tolerance:.3e}",
              file=sys.stderr)
        sys.exit(1)

    print("PASS")


if __name__ == '__main__':
    main()
