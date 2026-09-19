#!/usr/bin/env python3
"""
Step 5 of the blunt-TE XCUT meshing pipeline: merge the trimmed main
mesh (step 3) and the TE-face gap-filler box (step 4) into a final
2-block p3d + nmf.

Connectivity (validated interactively before being scripted):
  - Block 1 (main) i=1 and i=imax planes: FARFIELD (wake tips, x=Lfar)
  - Block 1 j=1 plane, split three ways:
      i=1..i_bot_te            <-> Block 2 j=1 row (ONE_TO_ONE, reversed
                                    direction -- main's bottom wake row
                                    runs Lfar->1, box's own i runs 1->Lfar)
      i=i_bot_te..i_top_te      VISCOUS (the real airfoil surface)
      i=i_top_te..imax         <-> Block 2 j=jmax row (ONE_TO_ONE, same
                                    direction -- both run 1->Lfar)
  - Block 1 j=jmax plane: FARFIELD (outer boundary)
  - Block 2 i=1 plane (the TE face itself): VISCOUS
  - Block 2 i=imax plane (far end of the gap strip, x=Lfar): FARFIELD,
    consistent with Block 1's own wake-tip treatment at the same x

Usage:
    python3 step5_merge_blocks.py \\
        --trimmed-p3d oat15_xcut_trimmed.p3d \\
        --box-p3d oat15_te_box.p3d \\
        --meta oat15_te_box.meta.json \\
        --out-prefix oat15_final
"""
import argparse
import os

import numpy as np

from xcut_common import (read_p3d, write_p3d_multiblock, read_meta,
                          read_p3d_function, write_p3d_function_multiblock)


NMF_HEADER = """# ==================== Neutral Map File (xcut_pipeline step 5, final 2-block) ====================
# ==================== ============================================================ ====================
# Block#   IDIM   JDIM   KDIM
# -----------------------------------------------------------------------------------
       2

       1 {imax1:6d} {jmax1:6d}      1
       2 {imax2:6d} {jmax2:6d}      1

# ===================================================================================
# Type         B1  F1     S1   E1     S2   E2    B2  F2     S1   E1     S2   E2  Swap
#                                                              Compute forces (walls)
# -----------------------------------------------------------------------------------
"""


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--trimmed-p3d', required=True)
    ap.add_argument('--box-p3d', required=True)
    ap.add_argument('--meta', required=True, help='step4 .meta.json sidecar')
    ap.add_argument('--out-prefix', required=True)
    args = ap.parse_args()

    meta = read_meta(args.meta)
    i_bot_te, i_top_te = meta['i_bot_te'], meta['i_top_te']

    (x1, y1), = read_p3d(args.trimmed_p3d)
    (x2, y2), = read_p3d(args.box_p3d)
    imax1, jmax1 = x1.shape
    imax2, jmax2 = x2.shape
    print(f"block 1 (main): {imax1} x {jmax1}")
    print(f"block 2 (box):  {imax2} x {jmax2}")

    # coordinate-verify both ONE_TO_ONE connections before writing anything
    bot1 = np.column_stack([x1[0:i_bot_te, 0], y1[0:i_bot_te, 0]])
    bot2 = np.column_stack([x2[:, 0], y2[:, 0]])[::-1]
    err_bot = np.max(np.abs(bot1 - bot2))

    top1 = np.column_stack([x1[i_top_te - 1:imax1, 0], y1[i_top_te - 1:imax1, 0]])
    top2 = np.column_stack([x2[:, jmax2 - 1], y2[:, jmax2 - 1]])
    err_top = np.max(np.abs(top1 - top2))

    print(f"ONE_TO_ONE coordinate check: bottom max mismatch={err_bot:.3g}, "
          f"top max mismatch={err_top:.3g}")
    tol = 1e-6
    if err_bot > tol or err_top > tol:
        raise SystemExit("ONE_TO_ONE connection coordinates do not match "
                          "within tolerance -- refusing to write a mesh with "
                          "an inconsistent interface; check --meta / box inputs")

    p3d_out = args.out_prefix + '.p3d'
    nmf_out = args.out_prefix + '.nmf'
    write_p3d_multiblock(p3d_out, [(x1, y1), (x2, y2)])

    with open(nmf_out, 'w') as f:
        f.write(NMF_HEADER.format(imax1=imax1, jmax1=jmax1, imax2=imax2, jmax2=jmax2))
        f.write(f"FARFIELD        1   1      1 {jmax1:4d}      1    1\n")
        f.write(f"FARFIELD        1   2      1 {jmax1:4d}      1    1\n")
        f.write(f"ONE_TO_ONE      1   3      1 {i_bot_te:4d}      1    1"
                f"     2   3  {imax2:4d}    1      1    1 FALSE\n")
        f.write(f"VISCOUS         1   3  {i_bot_te:4d} {i_top_te:4d}      1    1"
                f"                                    TRUE\n")
        f.write(f"ONE_TO_ONE      1   3  {i_top_te:4d} {imax1:4d}      1    1"
                f"     2   4      1 {imax2:4d}      1    1 FALSE\n")
        f.write(f"FARFIELD        1   4      1 {imax1:4d}      1    1\n")
        f.write(f"VISCOUS         2   1      1 {jmax2:4d}      1    1"
                f"                                    TRUE\n")
        f.write(f"FARFIELD        2   2      1 {jmax2:4d}      1    1\n")

    print(f"wrote {p3d_out}")
    print(f"wrote {nmf_out}")

    # Opportunistically merge quality stats too, if a *_stats.p3d sits
    # next to BOTH inputs (block 1's from Construct2D itself via step 3's
    # --stats-p3d passthrough, block 2's from step 4's own quality-stats
    # computation) -- into one 2-block stats file alongside the 2-block
    # mesh, so both blocks' quality is visible together (e.g. in
    # postpycess.py). Silently skipped if either is missing -- this was
    # always optional and older runs of this pipeline won't have block
    # 2's stats file.
    stats1_path = os.path.splitext(args.trimmed_p3d)[0] + '_stats.p3d'
    stats2_path = os.path.splitext(args.box_p3d)[0] + '_stats.p3d'
    have1, have2 = os.path.exists(stats1_path), os.path.exists(stats2_path)
    if have1 and have2:
        comments, si1, sj1, sk1, arrays1 = read_p3d_function(stats1_path)
        _, si2, sj2, sk2, arrays2 = read_p3d_function(stats2_path)
        if (si1, sj1) != (imax1, jmax1):
            print(f"NOTE: {stats1_path} dims ({si1}x{sj1}) don't match "
                  f"block 1 ({imax1}x{jmax1}) -- skipping merged quality stats")
        elif (si2, sj2) != (imax2, jmax2):
            print(f"NOTE: {stats2_path} dims ({si2}x{sj2}) don't match "
                  f"block 2 ({imax2}x{jmax2}) -- skipping merged quality stats")
        else:
            stats_out = args.out_prefix + '_stats.p3d'
            write_p3d_function_multiblock(
                stats_out, comments,
                [(imax1, jmax1, sk1), (imax2, jmax2, sk2)],
                [arrays1, arrays2])
            print(f"wrote {stats_out} (merged quality stats for both blocks, "
                  f"from {stats1_path} + {stats2_path})")
    else:
        missing = [p for p, have in [(stats1_path, have1), (stats2_path, have2)] if not have]
        print(f"no merged quality stats written (missing: {', '.join(missing)})")


if __name__ == '__main__':
    main()
