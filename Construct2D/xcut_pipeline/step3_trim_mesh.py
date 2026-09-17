#!/usr/bin/env python3
"""
Step 3 of the blunt-TE XCUT meshing pipeline: trim the degenerate
wake-tip-closure elements off the raw CGRD+NWKE=0+BUFF mesh (the ones
adjacent to the artificial (Lfar, 0.0) closure point from step 1), and
regenerate a matching single-block p3d + nmf.

Usage:
    python3 step3_trim_mesh.py \\
        --p3d oat15_xcut.p3d --nmf oat15_xcut.nmf \\
        --meta oat15_xcut.meta.json \\
        --nelm-clean 1 \\
        --out-prefix oat15_xcut_trimmed

`--nelm-clean N` removes the N elements adjacent to i=1 AND the N
elements adjacent to i=imax (i.e. 2N columns total). N=1 matches what
was needed for the OAT15/Lfar=20/Nwake=50 case validated so far; a
different geometry or closure-point placement may need more -- inspect
the raw mesh near the wake tips first if unsure.

Writes <out-prefix>.p3d, <out-prefix>.nmf, and <out-prefix>.meta.json
(carrying forward the updated bottomTE/topTE column indices that step
4/5 need -- never re-detected from geometry).
"""
import argparse
import os

from xcut_common import (
    read_p3d, write_p3d, read_meta, write_meta,
    read_p3d_function, write_p3d_function,
)


NMF_HEADER = """# ==================== Neutral Map File (xcut_pipeline step 3) ====================
# ==================== ============================================= ====================
# Block#   IDIM   JDIM   KDIM
# -----------------------------------------------------------------------------------
       1

       1 {imax:6d} {jmax:6d}      1

# ===================================================================================
# Type         B1  F1     S1   E1     S2   E2    B2  F2     S1   E1     S2   E2  Swap
#                                                              Compute forces (walls)
# -----------------------------------------------------------------------------------
"""


def write_single_block_nmf(fn, imax, jmax, i_bot_te, i_top_te):
    with open(fn, 'w') as f:
        f.write(NMF_HEADER.format(imax=imax, jmax=jmax))
        f.write(f"FARFIELD        1   1      1  {jmax:3d}      1    1\n")
        f.write(f"FARFIELD        1   2      1  {jmax:3d}      1    1\n")
        f.write(f"VISCOUS         1   3   {i_bot_te:4d} {i_top_te:4d}      1    1"
                f"                                    TRUE\n")
        f.write(f"FARFIELD        1   4      1  {imax:3d}      1    1\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--p3d', required=True, help='raw p3d from Construct2D (CGRD+NWKE=0+BUFF)')
    ap.add_argument('--nmf', help='raw nmf from Construct2D (read for reference/sanity only)')
    ap.add_argument('--stats-p3d', help='raw *_stats.p3d mesh-quality file from Construct2D '
                                         '(optional -- trimmed the same way as --p3d if given, '
                                         'so quality can still be inspected after trimming)')
    ap.add_argument('--meta', required=True, help='step1 .meta.json sidecar')
    ap.add_argument('--nelm-clean', type=int, required=True,
                     help='columns to remove from EACH end (i=1..N and imax-N+1..imax)')
    ap.add_argument('--out-prefix', required=True)
    args = ap.parse_args()

    meta = read_meta(args.meta)
    (x, y), = read_p3d(args.p3d)
    imax, jmax = x.shape
    print(f"raw mesh: {imax} x {jmax}")

    if imax != meta['n_total']:
        print(f"WARNING: raw p3d imax ({imax}) != step1 curve point count "
              f"({meta['n_total']}) -- check you're passing the matching "
              f"meta.json / p3d pair")

    n = args.nelm_clean
    if n < 0 or 2 * n >= imax:
        raise SystemExit(f"--nelm-clean {n} is out of range for imax={imax}")

    x_new = x[n:imax - n, :]
    y_new = y[n:imax - n, :]
    imax_new = x_new.shape[0]
    print(f"trimmed {n} column(s) off each end -> {imax_new} x {jmax}")

    i_bot_te_new = meta['i_bot_te_raw'] - n
    i_top_te_new = meta['i_top_te_raw'] - n
    if not (1 <= i_bot_te_new < i_top_te_new <= imax_new):
        raise SystemExit(
            f"trimmed bottomTE/topTE indices ({i_bot_te_new}, {i_top_te_new}) "
            f"fall outside the trimmed mesh (1..{imax_new}) -- nelm-clean too large?")
    print(f"new bottomTE index={i_bot_te_new}, topTE index={i_top_te_new} (1-based)")

    p3d_out = args.out_prefix + '.p3d'
    nmf_out = args.out_prefix + '.nmf'
    meta_out = args.out_prefix + '.meta.json'

    write_p3d(p3d_out, x_new, y_new)
    write_single_block_nmf(nmf_out, imax_new, jmax, i_bot_te_new, i_top_te_new)

    meta_extra = dict(raw_p3d=args.p3d, nelm_clean=n,
                       imax=imax_new, jmax=jmax,
                       i_bot_te=i_bot_te_new, i_top_te=i_top_te_new,
                       p3d_path=p3d_out)

    stats_out = None
    if args.stats_p3d:
        comments, s_imax, s_jmax, s_kmax, arrays = read_p3d_function(args.stats_p3d)
        if (s_imax, s_jmax) != (imax, jmax):
            raise SystemExit(f"--stats-p3d dims ({s_imax}x{s_jmax}) don't match "
                              f"--p3d dims ({imax}x{jmax}) -- mismatched pair?")
        arrays_new = [a[n:imax - n, :] for a in arrays]
        stats_out = args.out_prefix + '_stats.p3d'
        write_p3d_function(stats_out, comments, imax_new, jmax, s_kmax, arrays_new)
        print(f"trimmed stats file ({len(arrays)} field(s): "
              f"{comments[-1].lstrip('#') if len(comments) > 1 else ''})")
        meta_extra['stats_p3d_path'] = stats_out

    write_meta(meta_out, dict(meta, stage='step3', **meta_extra))
    print(f"wrote {p3d_out}")
    print(f"wrote {nmf_out}")
    if stats_out:
        print(f"wrote {stats_out}")
    print(f"wrote {meta_out}")


if __name__ == '__main__':
    main()
