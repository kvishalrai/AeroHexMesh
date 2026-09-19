#!/usr/bin/env python3
"""
Step 4 of the blunt-TE XCUT meshing pipeline: build the thin
wake-gap-filler "box" that closes the blunt TE face, matching the
resolution of the surrounding wake-arm rows.

Geometry (validated interactively before being scripted): the trimmed
mesh's own wall row (j=1) for i=1..i_bot_te (bottom wake arm, x: Lfar
-> 1) and i=i_top_te..imax (top wake arm, x: 1 -> Lfar) are re-used
AS-IS as the box's two long edges -- they are already fold-free,
already share nodes with the main mesh by construction, and (being
prescribed straight y=const lines) form a simple, non-twisting thin
strip rather than the main mesh's diverging radial O/C-grid columns.

The number of points spanning the gap (across the TE face) is chosen
so the spacing is UNIFORM and matches the mesh's own first radial
(wall-normal) grid spacing near the TE, per the user's requirement --
NOT a fixed count.

Usage:
    python3 step4_build_te_box.py \\
        --trimmed-p3d oat15_xcut_trimmed.p3d \\
        --meta oat15_xcut_trimmed.meta.json \\
        --out oat15_te_box.p3d
"""
import argparse
import os

import numpy as np

from xcut_common import (read_p3d, write_p3d, read_meta, write_meta,
                          compute_quality_stats, write_p3d_function)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--trimmed-p3d', required=True)
    ap.add_argument('--meta', required=True, help='step3 .meta.json sidecar')
    ap.add_argument('--out', required=True)
    args = ap.parse_args()

    meta = read_meta(args.meta)
    i_bot_te, i_top_te = meta['i_bot_te'], meta['i_top_te']  # 1-based

    (x, y), = read_p3d(args.trimmed_p3d)
    imax, jmax = x.shape
    print(f"trimmed mesh: {imax} x {jmax}, bottomTE i={i_bot_te}, topTE i={i_top_te}")

    # 0-based slices: bottom wake arm i=1..i_bot_te (1-based) -> [0:i_bot_te]
    xb = x[0:i_bot_te, 0]
    yb = y[0:i_bot_te, 0]
    # top wake arm i=i_top_te..imax (1-based) -> [i_top_te-1:imax]
    xt = x[i_top_te - 1:imax, 0]
    yt = y[i_top_te - 1:imax, 0]

    xb_r = xb[::-1]
    yb_r = yb[::-1]
    if len(xb_r) != len(xt):
        raise SystemExit(f"bottom wake arm ({len(xb_r)} pts) and top wake arm "
                          f"({len(xt)} pts) have different point counts -- "
                          f"expected them to match (both built with the same "
                          f"--nwake in step 1)")
    dx_mismatch = np.max(np.abs(xb_r - xt))
    if dx_mismatch > 1e-6:
        print(f"WARNING: bottom/top wake x-distributions differ by up to "
              f"{dx_mismatch:.3g} -- box will still be built (bilinear in y "
              f"at each i) but edges won't be perfectly straight")

    N_I = len(xt)
    x_col = 0.5 * (xb_r + xt)  # shared x at each i (== either side if symmetric)

    # gap half-width at each i (varies very slightly along the strip)
    half_gap = 0.5 * (yt - yb_r)
    print(f"gap half-width range along strip: {half_gap.min():.6g} .. {half_gap.max():.6g}")

    # first radial (wall-normal) spacing of the surrounding mesh, measured
    # at the two TE-corner columns (0-based index i_bot_te-1 / i_top_te-1)
    dy_bot = y[i_bot_te - 1, 1] - y[i_bot_te - 1, 0]
    dy_top = y[i_top_te - 1, 1] - y[i_top_te - 1, 0]
    d_wall = 0.5 * (abs(dy_bot) + abs(dy_top))
    print(f"first radial (wall-normal) spacing near TE: bottom={dy_bot:.6g} "
          f"top={dy_top:.6g} -> using d_wall={d_wall:.6g}")

    gap_height = yt[0] - yb_r[0]  # full gap at the TE face itself (x=1)
    n_j = max(2, int(round(gap_height / d_wall)) + 1)
    print(f"TE-face gap height={gap_height:.6g} -> {n_j} points across the gap "
          f"(uniform spacing {gap_height / (n_j - 1):.6g}, target was {d_wall:.6g})")

    x_box = np.tile(x_col[:, None], (1, n_j))
    y_box = np.zeros((N_I, n_j))
    for i in range(N_I):
        y_box[i, :] = np.linspace(yb_r[i], yt[i], n_j)

    # fold check
    bad = 0
    for i in range(N_I - 1):
        for j in range(n_j - 1):
            v1x, v1y = x_box[i + 1, j] - x_box[i, j], y_box[i + 1, j] - y_box[i, j]
            v2x, v2y = x_box[i, j + 1] - x_box[i, j], y_box[i, j + 1] - y_box[i, j]
            if v1x * v2y - v1y * v2x <= 0:
                bad += 1
    total = (N_I - 1) * (n_j - 1)
    print(f"fold check: {bad} / {total} folded or zero-area cells")
    if bad:
        raise SystemExit("box has folded cells -- refusing to write output; "
                          "inspect the wake-arm geometry before proceeding")

    write_p3d(args.out, x_box, y_box)
    meta_out = os.path.splitext(args.out)[0] + '.meta.json'
    write_meta(meta_out, dict(meta, stage='step4',
                               trimmed_p3d=args.trimmed_p3d,
                               box_imax=N_I, box_jmax=n_j,
                               box_path=args.out))
    print(f"wrote {args.out} ({N_I} x {n_j}, 0 folded cells)")
    print(f"wrote {meta_out}")

    # Quality stats: the box never goes through Construct2D (it's built
    # directly here), so it never gets a *_stats.p3d from Construct2D the
    # way the main mesh does. Compute the same skew-angle / xi-growth /
    # eta-growth metrics ourselves (Python port of Construct2D's own
    # compute_quality_stats, verified against real Construct2D output --
    # see xcut_common.compute_quality_stats) so the box's quality is
    # visible too, in the same format as any other *_stats.p3d.
    skewang, growthz, growthn = compute_quality_stats(x_box, y_box)
    stats_out = os.path.splitext(args.out)[0] + '_stats.p3d'
    write_p3d_function(stats_out,
                        ['#Grid quality information',
                         '#skew angle, xi-growth, eta-growth'],
                        N_I, n_j, 1, [skewang, growthz, growthn])
    print(f"wrote {stats_out} (max skew={skewang.max():.4f} deg, "
          f"max |xi-growth|={np.abs(growthz).max():.4f}, "
          f"max |eta-growth|={np.abs(growthn).max():.4f})")


if __name__ == '__main__':
    main()
