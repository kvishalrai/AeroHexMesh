#!/usr/bin/env python3
"""
One-off assembly script: stitches the 4 raw digitized-point files for the
higher-fidelity OAT15 refinement (2 segments x top/bottom each) into a
single, ordered, full-precision airfoil .dat ready for Construct2D (and
for xcut_pipeline's step1_build_curve.py).

Usage (defaults match this repo's own layout):
    python3 assemble_oat15_refined.py \\
        --source-dir ../sample_airfoils/oat15_refined_source \\
        --out ../sample_airfoils/oat15_refined.dat
"""
import argparse
import os

import numpy as np

from xcut_common import write_airfoil_dat


def load(fn):
    with open(fn) as f:
        lines = [l for l in f if l.strip() and not l.strip().startswith('#')]
    return np.array([[float(v) for v in l.split()] for l in lines])


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--source-dir', required=True,
                     help='directory containing the 4 puntos_xy_segmento*.dat files')
    ap.add_argument('--out', required=True, help='output airfoil .dat path')
    args = ap.parse_args()

    s1a = load(os.path.join(args.source_dir, 'puntos_xy_segmento1_abajo.dat'))   # bottom, x<=0.07
    s1r = load(os.path.join(args.source_dir, 'puntos_xy_segmento1_arriba.dat'))  # top,    x<=0.07
    s2a = load(os.path.join(args.source_dir, 'puntos_xy_segmento2_abajo.dat'))   # bottom, x>=0.07
    s2r = load(os.path.join(args.source_dir, 'puntos_xy_segmento2_arriba.dat'))  # top,    x>=0.07

    assert np.array_equal(s1a[0], s1r[0]), "seg1 bottom/top don't share the same LE point"

    # bottom path: TE -> LE  (seg2_abajo is LE->TE, seg1_abajo is LE->x=0.07; reverse both)
    bottom = np.vstack([s2a[::-1], s1a[::-1]])
    # top path: LE -> TE, excluding the LE point already supplied by `bottom`'s tail
    top = np.vstack([s1r[1:], s2r])

    curve = np.vstack([bottom, top])
    print(f"bottom path: {len(bottom)} pts (TE -> LE), top path: {len(top)} pts (LE -> TE, "
          f"LE point deduplicated)")
    print(f"assembled curve: {len(curve)} pts total")
    print(f"start (bottomTE): {curve[0]}   end (topTE): {curve[-1]}")
    print(f"LE point (junction): {bottom[-1]}")

    # monotonicity / sanity checks
    d_bottom_x = np.diff(bottom[:, 0])
    d_top_x = np.diff(top[:, 0])
    print(f"bottom path x monotonic decreasing: {np.all(d_bottom_x <= 0)} "
          f"(min step {d_bottom_x.max():.2e}, i.e. smallest |decrease|)")
    print(f"top path x monotonic increasing: {np.all(d_top_x >= 0)}")

    # arc-length point-spacing analysis along the WHOLE assembled curve
    seg = np.hypot(np.diff(curve[:, 0]), np.diff(curve[:, 1]))
    s = np.concatenate([[0], np.cumsum(seg)])
    print()
    print(f"total arc length: {s[-1]:.6f}")
    print(f"point spacing: min={seg.min():.6e}  max={seg.max():.6e}  "
          f"mean={seg.mean():.6e}  max/min ratio={seg.max()/seg.min():.2f}")

    # where is spacing finest/coarsest, and how does it vary segment to segment
    regions = [
        ("bottom seg2 (TE..x=0.07)", 0, len(s2a) - 1),
        ("bottom seg1 (x=0.07..LE)", len(s2a) - 1, len(bottom) - 1),
        ("top seg1 (LE..x=0.07)", len(bottom) - 1, len(bottom) - 1 + len(top[:len(s1r)-1])),
        ("top seg2 (x=0.07..TE)", len(bottom) - 1 + len(s1r) - 1, len(curve) - 1),
    ]
    for name, i0, i1 in regions:
        sub = seg[i0:i1]
        if len(sub):
            print(f"  {name:28s} n={len(sub):4d}  spacing min={sub.min():.3e} "
                  f"max={sub.max():.3e} mean={sub.mean():.3e}")

    # spacing right at the two segment-boundary "joins" (already measured above but
    # repeat explicitly here since these are the only points not from a single
    # uniform digitization)
    i_join_bottom = len(s2a) - 1
    i_join_top = len(bottom) - 1 + len(s1r) - 1
    print()
    print(f"spacing AT the seg1/seg2 boundary joins (bottom): {seg[i_join_bottom]:.6e} "
          f"(neighbors: {seg[i_join_bottom-1]:.3e}, {seg[i_join_bottom+1]:.3e})")
    print(f"spacing AT the seg1/seg2 boundary joins (top):    {seg[i_join_top]:.6e} "
          f"(neighbors: {seg[i_join_top-1]:.3e}, {seg[i_join_top+1]:.3e})")

    # fold/self-intersection sanity: since this is a single open curve (not a grid),
    # just check no zero-length or backtracking segments beyond the monotonicity
    # check already done, and that consecutive segment length doesn't jump by an
    # extreme factor (sign of a duplicated or dropped point)
    ratios = seg[1:] / seg[:-1]
    bad = np.where((ratios > 5) | (ratios < 0.2))[0]
    print()
    print(f"segments with >5x jump in spacing vs neighbor: {len(bad)}"
          + (f" at curve indices {bad.tolist()}" if len(bad) else ""))

    write_airfoil_dat(args.out, "OAT15_refined", curve[:, 0], curve[:, 1])
    print()
    print(f"wrote {args.out}")


if __name__ == '__main__':
    main()
