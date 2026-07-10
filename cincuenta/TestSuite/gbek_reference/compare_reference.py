#!/usr/bin/env python3
"""
Compare the exact atomic-limit reference -i*Delta^+_<(t,t') (produced by
gbek_selfconsistency.py) against cincuenta's own rank-L Cholesky
approximation (dumped by ImpuritySolverNeqGBEK::dumpPlusBath, e.g.
"gebk-fig3-L3-plus-bath-lesser").

Usage:
    python3 compare_reference.py gbek-atomic-limit-exact-lesser \\
        gebk-fig3-L3-plus-bath-lesser --tmax 4.0
"""
import argparse
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from provenance import write_provenance


def read_lesser_file(path):
    rows = []
    with open(path) as fh:
        for line in fh:
            cols = line.split()
            if len(cols) < 4:
                continue
            rows.append([float(c) for c in cols[:4]])
    rows = np.array(rows)
    ts = np.unique(rows[:, 0])
    N = len(ts)
    re = np.zeros((N, N))
    im = np.zeros((N, N))
    idx = {t: i for i, t in enumerate(ts)}
    for t, tp, r, i in rows:
        n, j = idx[t], idx[tp]
        re[n, j] = r
        im[n, j] = i
    return ts, re, im


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("reference", help="exact reference file (gbek_selfconsistency.py output)")
    ap.add_argument("approx", help="cincuenta's rank-L Cholesky output (dumpPlusBath)")
    ap.add_argument("--tmax", type=float, default=None)
    ap.add_argument("--rank", type=int, default=None,
                    help="Cholesky rank L used for --approx, shown in its plot title")
    ap.add_argument("--out", default="gbek_reference_comparison.png")
    args = ap.parse_args()

    ts_ref, re_ref, im_ref = read_lesser_file(args.reference)
    ts_app, re_app, im_app = read_lesser_file(args.approx)

    fig, axes = plt.subplots(2, 3, figsize=(15, 8))

    vmax = max(np.abs(re_ref).max(), np.abs(re_app).max())
    im0 = axes[0, 0].imshow(re_ref, origin="lower", extent=[ts_ref[0], ts_ref[-1]] * 2,
                             vmin=-vmax, vmax=vmax, cmap="RdBu_r")
    axes[0, 0].set_title("Exact reference: Re[-i Lambda^<_+(t,t')]")
    plt.colorbar(im0, ax=axes[0, 0])

    im1 = axes[0, 1].imshow(re_app, origin="lower", extent=[ts_app[0], ts_app[-1]] * 2,
                             vmin=-vmax, vmax=vmax, cmap="RdBu_r")
    rank_label = f"rank-{args.rank}" if args.rank is not None else "rank-L"
    axes[0, 1].set_title(f"cincuenta {rank_label} Cholesky approx")
    plt.colorbar(im1, ax=axes[0, 1])

    # Interpolate approx onto the reference grid for a residual map, if grids differ
    if not np.allclose(ts_ref, ts_app):
        from scipy.interpolate import RegularGridInterpolator
        interp = RegularGridInterpolator((ts_app, ts_app), re_app,
                                          bounds_error=False, fill_value=np.nan)
        grid_n, grid_j = np.meshgrid(ts_ref, ts_ref, indexing="ij")
        re_app_interp = interp((grid_n, grid_j))
    else:
        re_app_interp = re_app

    resid = re_ref - re_app_interp
    # Symmetric about 0 (so white == no error) and on the SAME scale as the
    # exact/approx plots to its left, so the residual's size relative to the
    # actual values is directly visible instead of auto-scaling to whatever
    # (typically much smaller) range the residual itself spans.
    im2 = axes[0, 2].imshow(resid, origin="lower", extent=[ts_ref[0], ts_ref[-1]] * 2,
                             vmin=-vmax, vmax=vmax, cmap="RdBu_r")
    axes[0, 2].set_title("Residual (exact - approx)")
    plt.colorbar(im2, ax=axes[0, 2])

    # Diagonal comparison
    diag_ref = np.diag(re_ref)
    diag_app = np.diag(re_app)
    axes[1, 0].plot(ts_ref, diag_ref, label="exact reference", lw=2)
    axes[1, 0].plot(ts_app, diag_app, "--", label="cincuenta approx", lw=2)
    axes[1, 0].axhline(0.5, color="gray", ls=":", lw=1, label="exact p-h value (0.5)")
    axes[1, 0].set_xlabel("t")
    axes[1, 0].set_ylabel("Re[-i Lambda(t,t)]")
    axes[1, 0].set_title("Diagonal")
    axes[1, 0].legend(fontsize=8)

    # A couple of off-diagonal (t, t'=fixed) slices
    for tp_fixed, ax_idx in [(0.4, 1), (1.5, 2)]:
        ax = axes[1, ax_idx]
        j_ref = np.argmin(np.abs(ts_ref - tp_fixed))
        j_app = np.argmin(np.abs(ts_app - tp_fixed))
        ax.plot(ts_ref, re_ref[:, j_ref], label="exact reference", lw=2)
        ax.plot(ts_app, re_app[:, j_app], "--", label="cincuenta approx", lw=2)
        ax.set_xlabel("t")
        ax.set_title(f"Slice at t'={tp_fixed}")
        ax.legend(fontsize=8)

    for ax in axes.flat:
        if args.tmax:
            ax.set_xlim(0, args.tmax)
            if hasattr(ax, "set_ylim") and ax in axes[0, :]:
                ax.set_ylim(0, args.tmax)

    fig.tight_layout()
    fig.savefig(args.out, dpi=150)
    print(f"Wrote {args.out}")
    prov = write_provenance(args.out, extra_files=[args.reference, args.approx],
                             notes=f"reference={args.reference} approx={args.approx}")
    print(f"Wrote {prov}")

    print("\nSummary:")
    print(f"  max|diag_ref - 0.5| (post t=0.4): "
          f"{np.max(np.abs(diag_ref[ts_ref > 0.4] - 0.5)):.4f}")
    print(f"  max|diag_app - 0.5| (post t=0.4): "
          f"{np.max(np.abs(diag_app[ts_app > 0.4] - 0.5)):.4f}")
    print(f"  max|residual|: {np.nanmax(np.abs(resid)):.4f}")


if __name__ == "__main__":
    main()
