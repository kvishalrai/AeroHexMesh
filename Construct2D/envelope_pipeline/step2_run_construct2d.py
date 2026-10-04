#!/usr/bin/env python3
"""
Step 2 of the 30P-30N envelope pipeline: run Construct2D on the step-1
envelope curve, non-interactively and reproducibly.

Needs the XCUT build of Construct2D (this repo's src/, PLAN.md step 0):
NWKE=0 keeps the curve's coordinates as given, accepts the open downstream
end, and treats i=1/imax as the two cut columns.

Every setting goes through a grid_options.in namelist (Construct2D reads it
at startup), and the only keystrokes piped in are y (accept CGRD for an open curve) / GRID /
BUFF / QUIT. The run
happens in a fresh temporary folder, so nothing next to the curve can be
picked up or overwritten by accident. Outputs are then copied to
<out-prefix>.p3d / .nmf / _stats.p3d plus a .meta.json carrying the step-1
meta forward with the exact settings used.

Checks after the run: the mesh's j=1 row is bit-identical to the input
curve, and no cell has a non-positive area.

Usage:
    python3 step2_run_construct2d.py --curve 30p30n_envelope_v1.dat \\
        --construct2d ../construct2d --out-prefix 30p30n_envelope_v1_c2d
"""
import argparse
import os
import shutil
import subprocess
import tempfile
import numpy as np

from env_common import (read_curve_dat, read_meta, write_meta, refuse_existing,
                        read_p3d)


def cell_areas(x, y):
    """Signed area of every quad cell (i, j) -> (i+1, j+1)."""
    x00, y00 = x[:-1, :-1], y[:-1, :-1]
    x10, y10 = x[1:, :-1], y[1:, :-1]
    x11, y11 = x[1:, 1:], y[1:, 1:]
    x01, y01 = x[:-1, 1:], y[:-1, 1:]
    return 0.5 * ((x11 - x00) * (y01 - y10) - (y11 - y00) * (x01 - x10))


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--curve', required=True, help='step-1 curve .dat')
    ap.add_argument('--construct2d', required=True,
                    help='path to a construct2d binary built from this repo (XCUT)')
    ap.add_argument('--out-prefix', required=True, help='output path prefix (new files)')
    ap.add_argument('--jmax', type=int, default=100)
    ap.add_argument('--radi', type=float, default=15.0,
                    help='HYPR marches about 1.2 x RADI from the curve')
    ap.add_argument('--ypls', type=float, default=1.0)
    ap.add_argument('--recd', type=float, default=1.7e6,
                    help='Reynolds number for the first-cell height (30P-30N '
                         'BANC case: 1.71e6 on the stowed chord)')
    ap.add_argument('--cfrc', type=float, default=0.5,
                    help='chord fraction used in the y+ wall-distance estimate')
    ap.add_argument('--alfa', type=float, default=1.0)
    ap.add_argument('--epsi', type=float, default=15.0)
    ap.add_argument('--epse', type=float, default=0.0)
    ap.add_argument('--funi', type=float, default=0.2)
    ap.add_argument('--asmt', type=int, default=20)
    ap.add_argument('--xdamp', type=float, default=0.0,
                    help='XCUT pull of the cut columns toward the TE x-location')
    args = ap.parse_args()

    out = {k: args.out_prefix + ext for k, ext in
           (('p3d', '.p3d'), ('nmf', '.nmf'), ('stats', '_stats.p3d'),
            ('meta', '.meta.json'), ('log', '.log'))}
    refuse_existing(*out.values())

    curve = read_curve_dat(args.curve)
    meta_in = os.path.splitext(args.curve)[0] + '.meta.json'
    meta = read_meta(meta_in) if os.path.exists(meta_in) else {}
    exe = os.path.abspath(args.construct2d)
    name = 'env'

    nml = f"""&SOPT
  nwke = 0
  radi = {args.radi!r}
  xdamp = {args.xdamp!r}
/
&VOPT
  name = '{name}'
  jmax = {args.jmax}
  slvr = 'HYPR'
  topo = 'CGRD'
  ypls = {args.ypls!r}
  recd = {args.recd!r}
  cfrc = {args.cfrc!r}
  alfa = {args.alfa!r}
  epsi = {args.epsi!r}
  epse = {args.epse!r}
  funi = {args.funi!r}
  asmt = {args.asmt}
/
&OOPT
  gdim = 2
/
"""
    # 'y' answers Construct2D's "C-grid for a blunt TE?" question, which it
    # asks because the envelope's downstream end is open.
    keys = 'y\nGRID\nBUFF\nQUIT\n'

    with tempfile.TemporaryDirectory() as tmp:
        shutil.copy(args.curve, os.path.join(tmp, name + '.dat'))
        with open(os.path.join(tmp, 'grid_options.in'), 'w') as f:
            f.write(nml)
        res = subprocess.run([exe, name + '.dat'], cwd=tmp, input=keys,
                             capture_output=True, text=True)
        with open(out['log'], 'w') as f:
            f.write(res.stdout + res.stderr)
        produced = {k: os.path.join(tmp, name + ext) for k, ext in
                    (('p3d', '.p3d'), ('nmf', '.nmf'), ('stats', '_stats.p3d'))}
        missing = [p for p in produced.values() if not os.path.exists(p)]
        if missing:
            raise SystemExit(f"Construct2D did not produce {missing}; see {out['log']}")
        for k, p in produced.items():
            shutil.copy(p, out[k])

    (x, y), = read_p3d(out['p3d'])
    imax, jmax = x.shape
    print(f"Mesh: imax={imax} jmax={jmax} (curve has {len(curve)} points)")

    failures = []
    if imax != len(curve):
        failures.append(f"imax {imax} != curve points {len(curve)}")
    elif not (np.array_equal(x[:, 0], curve[:, 0]) and np.array_equal(y[:, 0], curve[:, 1])):
        err = max(np.abs(x[:, 0] - curve[:, 0]).max(), np.abs(y[:, 0] - curve[:, 1]).max())
        failures.append(f"wall row is not bit-identical to the curve (max diff {err:.3e})")
    else:
        print("Wall row: bit-identical to the input curve")

    A = cell_areas(x, y)
    # The curve runs clockwise around the body, so positive orientation is
    # whichever sign the bulk of the cells have; any cell of the other sign
    # (or zero) is folded.
    sgn = np.sign(np.median(A))
    bad = np.argwhere(A * sgn <= 0)
    print(f"Cells: {A.size}, folded or zero-area: {len(bad)}")
    if len(bad):
        jfirst = np.bincount(bad[:, 1]).nonzero()[0][:5]
        print(f"  first affected j-levels: {jfirst.tolist()}; "
              f"i range {bad[:, 0].min()}..{bad[:, 0].max()}")

    meta.update({
        'stage': 'step2',
        'step1_meta': meta_in if os.path.exists(meta_in) else None,
        'construct2d': exe,
        'grid_options_in': nml,
        'keystrokes': keys,
        'imax': int(imax), 'jmax': int(jmax),
        'n_folded_cells': int(len(bad)),
        'first_cell_height': float(np.hypot(x[:, 1] - x[:, 0], y[:, 1] - y[:, 0]).min()),
        'p3d': out['p3d'], 'nmf': out['nmf'], 'stats_p3d': out['stats'],
        'failures': failures,
    })
    write_meta(out['meta'], meta)
    for k in ('p3d', 'nmf', 'stats', 'meta', 'log'):
        print(f"Wrote {out[k]}")
    if failures:
        print("\nFAILED:")
        for f_ in failures:
            print("  - " + f_)
        raise SystemExit(1)


if __name__ == '__main__':
    main()
