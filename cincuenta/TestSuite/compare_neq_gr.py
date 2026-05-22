#!/usr/bin/env python3
"""Compare G^>(t,0) and G^R(t,0) between a tDMRG run and an exact-diag run.

Usage:
    python3 compare_neq_gr.py \\
        --tdmrg-retarded  green-retarded-tdmrg \\
        --ed-retarded     green-retarded-ed \\
        --ed-lesser       green-lesser-ed

The tDMRG stores G^R(t,0) ≈ G^>(t,0) (particle sector only, no hole
correction).  The exact-diag stores the full G^R(t,0) = G^>-G^<.  This
script compares both the particle sector (G^>) and the full retarded
function, and prints per-time errors plus summary statistics.

File format (KadanoffBaym::dump output):
    t  t'  Re(G)  Im(G)
one row per (t,t') pair, all pairs written.
"""

import argparse
import math
import sys


def read_t0_slice(path):
    """Return dict {t: (re, im)} for all rows where t'=0."""
    data = {}
    with open(path) as f:
        for line in f:
            cols = line.split()
            if len(cols) < 4:
                continue
            try:
                t, tp = float(cols[0]), float(cols[1])
                re, im = float(cols[2]), float(cols[3])
            except ValueError:
                continue
            if abs(tp) < 1e-9:
                data[round(t, 8)] = (re, im)
    return data


def add(a, b):
    return (a[0] + b[0], a[1] + b[1])


def absdiff(a, b):
    return (abs(a[0] - b[0]), abs(a[1] - b[1]))


def print_table(header, rows, col_headers):
    print(header)
    widths = [max(len(h), 6) for h in col_headers]
    fmt_h = "  ".join(f"{{:>{w}}}" for w in widths)
    fmt_r = "  ".join(f"{{:>{w}.5f}}" for w in widths[1:])
    print(fmt_h.format(*col_headers))
    print("-" * (sum(widths) + 2 * (len(widths) - 1)))
    for row in rows:
        t_str = f"{row[0]:>5.1f}"
        vals = "  ".join(f"{v:>{widths[i+1]}.5f}" for i, v in enumerate(row[1:]))
        print(f"{t_str}  {vals}")


def summarise(label, errors_re, errors_im):
    print(f"\n{label}")
    print(f"  Max  |dRe|={max(errors_re):.5f}  |dIm|={max(errors_im):.5f}")
    print(f"  Mean |dRe|={sum(errors_re)/len(errors_re):.5f}"
          f"  |dIm|={sum(errors_im)/len(errors_im):.5f}")
    rms_re = math.sqrt(sum(e**2 for e in errors_re) / len(errors_re))
    rms_im = math.sqrt(sum(e**2 for e in errors_im) / len(errors_im))
    print(f"  RMS  |dRe|={rms_re:.5f}  |dIm|={rms_im:.5f}")


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--tdmrg-retarded", default="green-retarded-tdmrg",
                   metavar="FILE", help="tDMRG green-retarded file (stores G^>)")
    p.add_argument("--ed-retarded",    default="green-retarded-ed",
                   metavar="FILE", help="Exact-diag green-retarded file (stores G^R)")
    p.add_argument("--ed-lesser",      default="green-lesser-ed",
                   metavar="FILE", help="Exact-diag green-lesser file (stores G^<)")
    args = p.parse_args()

    td_ggt  = read_t0_slice(args.tdmrg_retarded)   # tDMRG: G^>(t,0)
    ed_gr   = read_t0_slice(args.ed_retarded)       # ED:    G^R(t,0)
    ed_gl   = read_t0_slice(args.ed_lesser)         # ED:    G^<(t,0)

    # Reconstruct ED G^>(t,0) = G^R(t,0) + G^<(t,0)
    ed_ggt = {t: add(ed_gr[t], ed_gl[t]) for t in ed_gr if t in ed_gl}

    times = sorted(set(td_ggt) & set(ed_gr) & set(ed_ggt))
    if not times:
        print("ERROR: no common time points found -- check file paths", file=sys.stderr)
        sys.exit(1)

    # ---- G^>(t,0) comparison ------------------------------------------------
    cols = ["t", "Re tDMRG", "Im tDMRG", "Re ED G^>", "Im ED G^>", "|dRe|", "|dIm|"]
    rows_ggt, dre_ggt, dim_ggt = [], [], []
    for t in times:
        r_td, i_td = td_ggt[t]
        r_ed, i_ed = ed_ggt[t]
        dr, di = abs(r_td - r_ed), abs(i_td - i_ed)
        dre_ggt.append(dr); dim_ggt.append(di)
        rows_ggt.append((t, r_td, i_td, r_ed, i_ed, dr, di))

    print_table("\nG^>(t,0): tDMRG vs exact-diag (particle sector, single-spin)",
                rows_ggt, cols)
    summarise("G^>(t,0) errors", dre_ggt, dim_ggt)

    # ---- G^R(t,0) comparison ------------------------------------------------
    # tDMRG approximates G^R ≈ G^> (no hole-sector correction yet).
    cols2 = ["t", "Re tDMRG", "Im tDMRG", "Re ED G^R", "Im ED G^R", "|dRe|", "|dIm|"]
    rows_gr, dre_gr, dim_gr = [], [], []
    for t in times:
        r_td, i_td = td_ggt[t]
        r_ed, i_ed = ed_gr[t]
        dr, di = abs(r_td - r_ed), abs(i_td - i_ed)
        dre_gr.append(dr); dim_gr.append(di)
        rows_gr.append((t, r_td, i_td, r_ed, i_ed, dr, di))

    print_table("\nG^R(t,0): tDMRG (=G^>, no G^< correction) vs exact-diag",
                rows_gr, cols2)
    summarise("G^R(t,0) errors (includes missing G^< off-diagonal)", dre_gr, dim_gr)

    # ---- G^< missing term ---------------------------------------------------
    # Show G^<(t,0) from ED to quantify the hole-sector correction needed.
    cols3 = ["t", "Re ED G^<", "Im ED G^<", "|G^<|"]
    rows_gl = []
    for t in times:
        if t not in ed_gl:
            continue
        r, im = ed_gl[t]
        rows_gl.append((t, r, im, math.hypot(r, im)))

    print_table("\nG^<(t,0) from exact-diag (hole-sector correction missing from tDMRG)",
                rows_gl, cols3)


if __name__ == "__main__":
    main()
