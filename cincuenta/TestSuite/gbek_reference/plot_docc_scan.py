#!/usr/bin/env python3
"""
Plot GBEK PRB 88, 235106 (2013) Fig. 7 or Fig. 8 (double occupation vs.
decomposition method), from the .npz files produced by run_fig7_scan.py /
run_fig8_scan.py, and print quantitative pass/fail diagnostics against the
paper's actual qualitative claims (not just "a plot got produced").

Fig. 7 (--figure 7): one line per DMFT iteration (1..7), top panel =
eigenvector decomposition, bottom panel = Cholesky decomposition, d(t) vs t,
for each L in --L (comma-separated, default "2"; L=4/Lbath=8 data is the
2026-07-10 follow-up investigating the Fig. 8 combo-A magnitude gap, see
investigate_L4_tmax4.py and run_fig7_scan.py --L 4). The paper's claim:
Cholesky is CAUSAL (an iteration only changes the curve beyond a growing
time t_n*, converging by t_max), eigenvector is NOT (changes occur at all t
on every iteration).

Fig. 8 (--figure 8): converged d(t) for the paper's 4 (L_bath, t_max)
combinations (A/B/C/D), same two-panel layout, PLUS combo E (L=4/Lbath=8,
t_max=4 -- not one of the paper's original 4, an investigative extra added
2026-07-10 to check whether combo A's large late-time growth is plain
rank-insufficiency; included automatically if its .npz files are present).
The paper's claim: eigenvector can converge to a visibly WRONG answer when
L_bath is undersized for t_max (combo A vs D); Cholesky agrees closely
instead, with the agreement window growing with L_bath.

Usage:
    uv run --with numpy --with matplotlib python3 plot_docc_scan.py --figure 7
    uv run --with numpy --with matplotlib python3 plot_docc_scan.py --figure 7 --L 2,4
    uv run --with numpy --with matplotlib python3 plot_docc_scan.py --figure 8
"""
import argparse
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from provenance import write_provenance

MODES = ("eigenvector", "cholesky")  # top panel, bottom panel, in this order
COMBOS = ["A", "B", "C", "D", "E"]  # E is an investigative extra, included if present


def causal_break_time(d_n, d_np1, ts, tol):
    """Largest t such that |d_n(t)-d_{n+1}(t)| < tol for ALL t' <= t (a
    "breaks and stays broken" scan, not a pointwise check at one t -- see
    plan Step 7). Returns (t_star, persistent) where persistent is False if
    the disagreement doesn't actually stay broken afterward (an anomaly
    worth reporting, not a coincidental one-off dip)."""
    diff = np.abs(d_n - d_np1)
    bad = diff >= tol
    if not bad.any():
        return ts[-1], True
    first_bad = int(np.argmax(bad))
    t_star = ts[first_bad - 1] if first_bad > 0 else ts[0] - (ts[1] - ts[0])
    # "stays broken": what fraction of points from first_bad onward are still bad?
    persistent_frac = bad[first_bad:].mean()
    return t_star, persistent_frac > 0.8


def plot_figure7(args):
    L_list = [int(x) for x in args.L.split(",")]
    data = {(L, mode): np.load(f"fig7_docc_L{L}_{mode}.npz")
            for L in L_list for mode in MODES}

    fig, axes = plt.subplots(len(L_list), 2, figsize=(11, 4 * len(L_list)),
                              sharex=True, squeeze=False)
    extra_files = ["plot_docc_scan.py"]

    for row, L in enumerate(L_list):
        n_iter = data[(L, "cholesky")]["docc_history"].shape[0]
        # Match the paper's own Fig. 7 convention: earlier (not-yet-converged)
        # iterations are dashed, in a color gradient; the final, most-converged
        # iteration is a solid red line that stands out from the rest.
        colors = plt.cm.viridis(np.linspace(0.1, 0.9, n_iter - 1))
        for col, mode in enumerate(MODES):
            ax = axes[row][col]
            ts = data[(L, mode)]["ts"]
            hist = data[(L, mode)]["docc_history"]
            for it in range(n_iter):
                is_final = it == n_iter - 1
                color = "red" if is_final else colors[it]
                linestyle = "-" if is_final else "--"
                lw = 2.2 if is_final else 1.3
                ax.plot(ts, hist[it], color=color, linestyle=linestyle, lw=lw,
                        label=f"iter {it+1}")
            ax.set_ylabel("d(t)")
            ax.set_title(f"{mode}, L={L} (Lbath={2*L})")
            ax.grid(alpha=0.2)
            if args.ylim is not None:
                ax.set_ylim(*args.ylim)
            if row == 0 and col == 0:
                ax.legend(fontsize=7, ncol=4, loc="upper left")
            extra_files.append(f"fig7_docc_L{L}_{mode}.npz")
    for col in range(2):
        axes[-1][col].set_xlabel("t")
    fig.suptitle(args.title)
    fig.tight_layout()
    fig.savefig(args.out, dpi=150)
    print(f"Wrote {args.out}")

    prov = write_provenance(args.out, extra_files=extra_files,
                             notes=f"figure=7 L_list={L_list} tol_docc={args.tol_docc}")
    print(f"Wrote {prov}")

    print("\n--- Fig. 7 physics criteria ---")
    for L in L_list:
        for mode in MODES:
            ts = data[(L, mode)]["ts"]
            hist = data[(L, mode)]["docc_history"]
            n_iter = hist.shape[0]
            scale = np.max(np.abs(hist))
            tol = args.tol_docc if args.tol_docc is not None else 0.02 * scale
            print(f"\nL={L} mode={mode}  (dynamic range max|d|={scale:.4f}, tol_docc={tol:.4f})")
            t_stars = []
            for it in range(n_iter - 1):
                t_star, persistent = causal_break_time(hist[it], hist[it + 1], ts, tol)
                t_stars.append(t_star)
                flag = "" if persistent else "  <-- NOT persistent (breaks then re-agrees; anomaly)"
                print(f"  iter {it+1}->{it+2}: t* = {t_star:.3f}{flag}")
            t_stars = np.array(t_stars)
            non_decreasing = np.all(np.diff(t_stars) >= -1e-9)
            print(f"  t* sequence non-decreasing across iterations: {non_decreasing}")
            print(f"  t*(final pair) = {t_stars[-1]:.3f} vs t_max = {ts[-1]:.3f} "
                  f"({'near t_max -- CAUSAL, converging as expected' if t_stars[-1] > 0.8*ts[-1] else 'well short of t_max'})")

    print("\nExpected per the paper: Cholesky's t* sequence is non-decreasing and "
          "reaches near t_max by iteration 7; eigenvector's does NOT (changes persist "
          "at small t even on later iterations). A result where BOTH modes show the "
          "same causal structure, or NEITHER does, contradicts the paper's central "
          "claim and should be investigated, not smoothed over.")


def plot_figure8(args):
    combos = []
    for combo in COMBOS:
        if all(Path(f"fig8_docc_{combo}_{mode}.npz").exists() for mode in MODES):
            combos.append(combo)
        else:
            print(f"NOTE: skipping combo {combo} -- .npz file(s) not found")

    data = {}
    for combo in combos:
        for mode in MODES:
            data[(combo, mode)] = np.load(f"fig8_docc_{combo}_{mode}.npz")

    fig, axes = plt.subplots(2, 1, figsize=(7, 8), sharex=True)

    # Color by L_bath (so same-L_bath combos at different t_max share a
    # color and become directly comparable), linestyle by t_max (so the
    # causal-property question -- does the same L_bath's t_max=4 curve lie
    # exactly on top of its t_max=2 curve over their shared range -- is
    # actually visible instead of one line silently occluding the other).
    L_values = sorted({int(data[(c, mode)]["L"]) for c in combos})
    L_colors = dict(zip(L_values, plt.cm.viridis(np.linspace(0.15, 0.85, len(L_values)))))
    t_max_styles = {2.0: "-", 4.0: "--"}

    for ax, mode in zip(axes, MODES):
        for combo in combos:
            d = data[(combo, mode)]
            L, t_max = int(d["L"]), float(d["t_max"])
            ax.plot(d["ts"], d["d_t"], color=L_colors[L],
                    linestyle=t_max_styles.get(t_max, ":"), lw=2.0,
                    label=f"{combo}: L_bath={2*L}, t_max={t_max:g}")
        ax.set_ylabel("d(t)")
        ax.set_title(mode)
        ax.grid(alpha=0.2)
        if args.ylim is not None:
            ax.set_ylim(*args.ylim)
    axes[0].legend(fontsize=8, loc="upper left")
    axes[-1].set_xlabel("t")
    fig.suptitle(args.title)
    fig.tight_layout()
    fig.savefig(args.out, dpi=150)
    print(f"Wrote {args.out}")

    extra_files = ["plot_docc_scan.py"] + [f"fig8_docc_{c}_{m}.npz" for c in combos for m in MODES]
    prov = write_provenance(args.out, extra_files=extra_files, notes=f"figure=8 combos={combos}")
    print(f"Wrote {prov}")

    print("\n--- Fig. 8 physics criteria ---")
    print("Note: the paper's claim is about EARLY-time agreement across L_bath, not "
          "whole-range agreement -- different L_bath are expected to diverge from each "
          "other at LATE t (that's the whole point: smaller L_bath's converged answer "
          "is only correct up to its own t*). So the metric below is an agreement "
          "HORIZON (largest t up to which two combos agree, breaks-and-stays-broken), "
          "not a max-over-the-whole-range discrepancy, which would conflate expected "
          "late-time spread with genuine early-time disagreement. Each pair is compared "
          "only over their OWN overlapping time range (the shorter of the two combos' "
          "own t_max), not always combo D's t_max<=2 range.")

    def interp_to(ts_target, ts_src, d_src):
        return np.interp(ts_target, ts_src, d_src)

    def agreement_horizon(ts, d_x, d_y, tol):
        """Same 'breaks and stays broken' idea as Fig. 7's causal_break_time,
        applied here across combos (different L_bath) instead of across
        iterations."""
        return causal_break_time(d_x, d_y, ts, tol)

    all_pairs = [("A", "D"), ("B", "C"), ("B", "D"), ("C", "D"), ("A", "E"), ("D", "E")]
    pairs = [(x, y) for x, y in all_pairs if x in combos and y in combos]

    for mode in MODES:
        scale = max(np.ptp(data[(c, mode)]["d_t"]) for c in combos)
        tol = args.tol_docc if args.tol_docc is not None else 0.02 * scale
        print(f"\nmode={mode}  (dynamic range ~{scale:.4f}, tol={tol:.4f})")
        for x, y in pairs:
            ts_x, ts_y = data[(x, mode)]["ts"], data[(y, mode)]["ts"]
            ts_common = ts_x if ts_x[-1] <= ts_y[-1] else ts_y  # shorter of the two combos' own range
            d_x = interp_to(ts_common, ts_x, data[(x, mode)]["d_t"])
            d_y = interp_to(ts_common, ts_y, data[(y, mode)]["d_t"])
            t_star, persistent = agreement_horizon(ts_common, d_x, d_y, tol)
            flag = "" if persistent else "  <-- NOT persistent (re-agrees later; anomaly)"
            print(f"  {x} vs {y} (t<={ts_common[-1]:g}): agreement horizon t* = {t_star:.3f}{flag}")

    print("\nExpected per the paper: Cholesky's agreement horizons should be "
          "substantially LARGER than eigenvector's for the same combo pairs (Cholesky "
          "stays correct/self-consistent across L_bath choices up to a growing t*; "
          "eigenvector can diverge across L_bath even at short times because its "
          "non-causal decomposition mixes in information from beyond t_max). If "
          "eigenvector's horizons are NOT smaller than Cholesky's, that contradicts "
          "the paper's claim and should be investigated, not smoothed over.")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--figure", choices=["7", "8"], required=True)
    ap.add_argument("--out", default=None)
    ap.add_argument("--title", default=None)
    ap.add_argument("--L", default="2",
                    help="Fig. 7 only: comma-separated L values to plot, one row per L "
                         "(default \"2\", the paper's own combo; pass \"2,4\" to also "
                         "include the L=4/Lbath=8 follow-up investigation)")
    ap.add_argument("--tol-docc", type=float, default=None,
                    help="absolute tolerance for the causal-break/agreement-horizon scan; "
                         "default is 2%% of that mode's own |d(t)| dynamic range")
    ap.add_argument("--ylim", default="paper",
                    help="'paper' (default) uses the actual paper axis range for direct "
                         "visual comparison (Fig. 7: 0-0.15, Fig. 8: 0-0.1, matching "
                         "page_12.png); 'auto' lets matplotlib autoscale; or pass "
                         "'ymin,ymax' explicitly")
    args = ap.parse_args()

    paper_ylim = {"7": (0.0, 0.15), "8": (0.0, 0.1)}[args.figure]
    if args.ylim == "paper":
        args.ylim = paper_ylim
    elif args.ylim == "auto":
        args.ylim = None
    else:
        args.ylim = tuple(float(x) for x in args.ylim.split(","))

    if args.figure == "7":
        args.out = args.out or "fig7_docc.png"
        args.title = args.title or "cf. GBEK Fig. 7: double occupation per DMFT iteration"
        plot_figure7(args)
    else:
        args.out = args.out or "fig8_docc.png"
        args.title = args.title or "cf. GBEK Fig. 8: converged double occupation vs. L_bath/t_max"
        plot_figure8(args)


if __name__ == "__main__":
    main()
