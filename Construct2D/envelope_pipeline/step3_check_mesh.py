#!/usr/bin/env python3
"""
Step 3 of the 30P-30N envelope pipeline: check the Construct2D mesh from
step 2 against the quality gates in PLAN.md (C8). Nothing is trimmed: the
downstream end is open and there is no upstream cap, so there are no
degenerate columns to remove.

Reported per curve region (the i-ranges from the step-1 meta), so a problem
can be traced to a bridge, a surface or a line:

  - folded / zero-area cells          (gate: none)
  - skew: max deviation from 90 deg   (gate: min cell angle >= --min-angle)
  - xi and eta growth: neighbour size ratio (gate: <= --max-growth)

The quality fields are computed with xcut_common.compute_quality_stats (a
port of Construct2D's own) and written to <out-prefix>_stats.p3d at full
precision, plus a <out-prefix>.meta.json carrying the earlier meta forward.

Usage:
    python3 step3_check_mesh.py --meta 30p30n_envelope_v1_c2d.meta.json \\
        --out-prefix 30p30n_envelope_v1_check
"""
import argparse
import numpy as np

from env_common import (read_meta, write_meta, read_p3d, refuse_existing,
                        compute_quality_stats)
from step2_run_construct2d import cell_areas


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--meta', required=True, help='step-2 .meta.json')
    ap.add_argument('--out-prefix', required=True)
    ap.add_argument('--min-angle', type=float, default=20.0)
    ap.add_argument('--max-growth', type=float, default=1.2)
    ap.add_argument('--jcheck', type=int, default=None,
                    help='only gate the first N j-levels (default: all)')
    args = ap.parse_args()

    out_stats = args.out_prefix + '_stats.p3d'
    out_meta = args.out_prefix + '.meta.json'
    refuse_existing(out_stats, out_meta)

    meta = read_meta(args.meta)
    (x, y), = read_p3d(meta['p3d'])
    imax, jmax = x.shape
    jg = jmax if args.jcheck is None else min(args.jcheck, jmax)

    skew, gz, gn = compute_quality_stats(x, y)
    A = cell_areas(x, y)
    folded = (A * np.sign(np.median(A))) <= 0

    # growth as a ratio (compute_quality_stats gives (L1-L0)/min(L0,L1))
    rz, rn = 1.0 + np.abs(gz), 1.0 + np.abs(gn)
    minang = 90.0 - skew

    regions = {k[2:]: v for k, v in meta.items()
               if k.startswith('i_') and not k.startswith('i_if_')
               and isinstance(v, list)}
    print(f"Mesh {imax} x {jmax}; gating j = 1..{jg}")
    print(f"{'region':20s} {'i-range':>11s} {'folded':>7s} {'min angle':>10s} "
          f"{'max xi-gr':>10s} {'max eta-gr':>10s}")
    rows, failures = {}, []
    for name, (a, b) in sorted(regions.items(), key=lambda kv: kv[1][0]):
        sl = slice(a - 1, b)
        r = dict(folded=int(folded[a - 1:min(b, imax - 1), :jg - 1].sum()),
                 min_angle=float(minang[sl, :jg].min()),
                 max_xi_growth=float(rz[sl, :jg].max()),
                 max_eta_growth=float(rn[sl, :jg].max()))
        rows[name] = r
        flag = []
        if r['folded']:
            flag.append('folded')
        if r['min_angle'] < args.min_angle:
            flag.append('angle')
        if max(r['max_xi_growth'], r['max_eta_growth']) > args.max_growth:
            flag.append('growth')
        if flag:
            failures.append(f"{name}: {', '.join(flag)}")
        print(f"{name:20s} {a:5d}-{b:<5d} {r['folded']:7d} {r['min_angle']:10.2f} "
              f"{r['max_xi_growth']:10.3f} {r['max_eta_growth']:10.3f}"
              + ("   <-- " + ", ".join(flag) if flag else ""))

    worst = np.unravel_index(np.argmin(minang[:, :jg]), minang[:, :jg].shape)
    print(f"\nWorst cell angle {minang[worst]:.2f} deg at i={worst[0] + 1}, j={worst[1] + 1}"
          f" ({x[worst]:.4f}, {y[worst]:.4f})")

    with open(out_stats, 'w') as f:
        # same header as Construct2D's own *_stats.p3d
        f.write("#Grid quality information\n#skew angle, xi-growth, eta-growth\n")
        f.write(f"{imax} {jmax} 1 3\n")
        for arr in (skew, gz, gn):
            for v in arr.flatten(order='F'):
                f.write(f"{v:25.16E}\n")

    meta.update({'stage': 'step3', 'step2_meta': args.meta,
                 'gates': {'min_angle': args.min_angle, 'max_growth': args.max_growth,
                           'j_levels': jg},
                 'region_quality': rows, 'quality_failures': failures,
                 'check_stats_p3d': out_stats})
    write_meta(out_meta, meta)
    print(f"\nWrote {out_stats}\nWrote {out_meta}")
    if failures:
        print("\nGATES FAILED:")
        for f_ in failures:
            print("  - " + f_)
        raise SystemExit(1)
    print("\nAll gates passed.")


if __name__ == '__main__':
    main()
