#!/usr/bin/env python3
"""
Step 1 variant of the blunt-TE XCUT pipeline with a GRADED wake: the
wake-point growth ratio changes with downstream distance, instead of one
ratio fitted to a fixed --nwake.

Identical to step1_build_curve.py (same airfoil extraction, same wake
shape, same curve ordering, same meta keys, so steps 3-5 work unchanged)
except for how points are placed along the two wake lines:

    TE  <= x < --x1 : each step is (1 + p1/100) times the previous one
    --x1 <= x < --x2 : growth p2 %
    --x2 <= x <= Lfar: growth p3 %

The first step is the airfoil's own TE spacing d0 (as in step1). The
region a step belongs to is decided by the x-coordinate of the point the
step starts from. In the last region the ratio is re-solved slightly
(never above p3) so that the final point lands exactly on Lfar.

The number of wake points is therefore an OUTPUT, not an input; it is
printed and recorded in the meta file as `nwake` (per side, TE corner
included), like step1's.

Example (x = 1 %, y = 4 %, z = 20 %):
    python3 xcut_pipeline/step1_build_curve_graded_wake.py \\
        --airfoil sample_airfoils/oat15_refined.dat \\
        --lfar 40.0 --wake-type spline \\
        --growth 1 4 20 \\
        --out sample_airfoils/oat15_slc_xcut_spline_graded_1_4_20.dat
"""
import argparse
import os
import numpy as np

from xcut_common import (read_airfoil_dat, write_airfoil_dat, write_meta,
                          geometric_growth_rate, _hermite_points)
from step1_build_curve import extract_airfoil_surface


def graded_stations(x_of_s, L, d0, ratios, bounds):
    """Arc-length stations 0 = s_0 < ... < s_n = L along a wake line.

    x_of_s: callable, x-coordinate of the line at arc length s.
    ratios: (r1, r2, r3); bounds: (x1, x2). Step k+1 = step k * r, where r
    is chosen from the x of the step's start point.
    """
    r1, r2, r3 = ratios
    x1, x2 = bounds
    s, step = [0.0], d0
    # regions 1 and 2: march with the exact ratios
    while True:
        x = x_of_s(s[-1])
        if x >= x2:
            break
        if s[-1] + step >= L:
            raise SystemExit("wake ends before the last growth region starts -- "
                             "increase --lfar or lower --x2")
        s.append(s[-1] + step)
        step *= r1 if x < x1 else r2
    # region 3: fixed first step `step`, n intervals, last point exactly at L
    L3 = L - s[-1]
    if step >= L3:
        s.append(L)
        return np.array(s), r3
    n = int(np.ceil(np.log(1 + L3 * (r3 - 1) / step) / np.log(r3)))  # intervals at ratio r3
    n = max(n, 1)
    r3_fit = geometric_growth_rate(L3, step, n + 1) if n > 1 else 1.0
    for _ in range(n):
        s.append(s[-1] + step)
        step *= r3_fit
    s = np.array(s)
    s[-1] = L
    return s, r3_fit


def print_curve_stats(curve, n_af, nwake, is_open, x1, x2):
    """Print statistics of the assembled 1D curve (point spacing and growth),
    for the wake, the airfoil surface, their junctions and the whole curve.
    Growth = ratio of consecutive segment lengths, shown in % (|ratio-1|)."""
    def segs(P):
        return np.hypot(*np.diff(P, axis=0).T)

    def growth(ds):
        g = ds[1:] / ds[:-1]
        return 100 * (np.maximum(g, 1 / g) - 1)    # symmetric: shrink counts too

    off = 0 if is_open else 1
    i_bte, i_tte = off + nwake - 1, off + nwake - 1 + n_af - 1   # 0-based TE corners
    af = curve[i_bte:i_tte + 1]
    top = curve[i_tte:i_tte + nwake]                  # TE corner -> far tip
    bot = curve[off:i_bte + 1][::-1]                  # TE corner -> far tip
    dsa, dst, dsb = segs(af), segs(top), segs(bot)
    w = 72
    print("=" * w)
    print("1D curve statistics")
    print("=" * w)
    print(f"total points: {len(curve)}   (airfoil {n_af}, wake {nwake}/side incl. TE corner"
          f"{'' if is_open else ', closure point at both ends'})")
    print(f"\nairfoil surface: {n_af} pts, {n_af-1} segments")
    print(f"  spacing   min {dsa.min():.4e}  max {dsa.max():.4e}  mean {dsa.mean():.4e}")
    ga = growth(dsa)
    print(f"  growth    max {ga.max():.2f} %  (at airfoil point {i_bte + 2 + int(np.argmax(ga))}, "
          f"1-based curve index)")
    print(f"  TE spacing: bottom {dsa[0]:.4e}, top {dsa[-1]:.4e}")
    print(f"\nwake (per side; top shown, bottom differs by < {np.abs(dst-dsb).max():.1e}):")
    s = np.concatenate([[0], np.cumsum(dst)])
    print(f"  {nwake} pts, {nwake-1} segments, arc length {s[-1]:.4f}, "
          f"x from {top[0,0]:.4f} to {top[-1,0]:.4f}")
    print(f"  spacing   first {dst[0]:.4e}  last {dst[-1]:.4e}  min {dst.min():.4e}  max {dst.max():.4e}")
    gw = growth(dst)
    print(f"  growth    max {gw.max():.2f} %  over the whole wake")
    print(f"  {'region':<16}{'segments':>9}{'first ds':>12}{'last ds':>12}{'growth min':>12}{'growth max':>12}")
    xs = top[:-1, 0]
    for name, lo, hi in [(f'TE <= x < {x1:g}', -np.inf, x1), (f'{x1:g} <= x < {x2:g}', x1, x2),
                         (f'x >= {x2:g}', x2, np.inf)]:
        k = np.where((xs >= lo) & (xs < hi))[0]
        if len(k) == 0:
            print(f"  {name:<16}{0:>9}")
            continue
        g = 100 * (dst[k[1:]] / dst[k[:-1]] - 1) if len(k) > 1 else np.array([np.nan])
        print(f"  {name:<16}{len(k):>9}{dst[k[0]]:>12.4e}{dst[k[-1]]:>12.4e}"
              f"{np.nanmin(g):>11.2f}%{np.nanmax(g):>11.2f}%")
    print("\nTE junctions (first wake segment / last airfoil segment):")
    print(f"  bottom {dsb[0] / dsa[0]:.4f}   top {dst[0] / dsa[-1]:.4f}")
    if not is_open:
        dc = np.hypot(*(curve[1] - curve[0]))
        print(f"\nclosure: far tips are {dc:.4e} from the (Lfar, 0) closure point "
              f"(degenerate cells there are removed by step 3)")
        body = curve[1:-1]
    else:
        body = curve
    db = segs(body)
    gb = growth(db)
    print(f"\nwhole curve{'' if is_open else ' (excluding the closure point)'}: "
          f"spacing min {db.min():.4e} max {db.max():.4e}, max growth {gb.max():.2f} %")
    print("=" * w)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--airfoil', required=True)
    ap.add_argument('--lfar', type=float, required=True, help='absolute x of the wake far end')
    ap.add_argument('--growth', type=float, nargs=3, required=True, metavar=('P1', 'P2', 'P3'),
                    help='growth in percent for the three regions, e.g. 1 4 20')
    ap.add_argument('--x1', type=float, default=2.0, help='end of region 1, absolute x (default 2)')
    ap.add_argument('--x2', type=float, default=3.0, help='end of region 2, absolute x (default 3)')
    ap.add_argument('--wake-type', choices=['straight', 'spline'], default='spline')
    ap.add_argument('--pull-length', type=float, default=1.0,
                    help='(spline only) TE-bisector tangent reach, chords (default 1.0)')
    ap.add_argument('--open', action='store_true',
                    help='leave out the (Lfar, 0) closure point (needs a Construct2D build with the open-cut fix)')
    ap.add_argument('--out', required=True)
    args = ap.parse_args()

    ratios = tuple(1.0 + p / 100.0 for p in args.growth)
    if min(ratios) < 1.0:
        raise SystemExit("growth percentages must be >= 0")

    header, x_raw, y_raw = read_airfoil_dat(args.airfoil)
    xc, yc = extract_airfoil_surface(x_raw, y_raw)
    n_af = len(xc)
    y_mid = 0.5 * (yc[0] + yc[-1])
    h = 0.5 * (yc[-1] - yc[0])
    d0 = 0.5 * (np.hypot(xc[1] - xc[0], yc[1] - yc[0]) + np.hypot(xc[-1] - xc[-2], yc[-1] - yc[-2]))
    x_te = xc[0]
    if not (x_te < args.x1 < args.x2 < args.lfar):
        raise SystemExit(f"need TE x ({x_te}) < x1 < x2 < lfar")
    print(f"airfoil: {n_af} pts, TE x={x_te}, half-gap h={h:.7f}, d0={d0:.6g}")

    # wake centerline, sampled densely and parametrized by arc length
    if args.wake_type == 'straight':
        line = np.column_stack([np.linspace(x_te, args.lfar, 20001), np.full(20001, y_mid)])
        y_far = y_mid
    else:
        t_bot = np.array([xc[0] - xc[1], yc[0] - yc[1]])
        t_top = np.array([xc[-1] - xc[-2], yc[-1] - yc[-2]])
        b = t_bot / np.linalg.norm(t_bot) + t_top / np.linalg.norm(t_top)
        b /= np.linalg.norm(b)
        if b[0] < 0:
            raise SystemExit("TE bisector points upstream -- check curve orientation")
        y_far = 0.0
        line = _hermite_points(np.linspace(0, 1, 20001), np.array([x_te, y_mid]),
                               np.array([args.lfar, y_far]), b * args.pull_length,
                               np.array([1.0, 0.0]) * args.pull_length)
    arc = np.concatenate([[0.0], np.cumsum(np.hypot(*np.diff(line, axis=0).T))])
    L = arc[-1]
    x_of_s = lambda s: np.interp(s, arc, line[:, 0])

    s, r3_fit = graded_stations(x_of_s, L, d0, ratios, (args.x1, args.x2))
    cl = np.column_stack([np.interp(s, arc, line[:, 0]), np.interp(s, arc, line[:, 1])])
    cl[0] = (x_te, y_mid)
    cl[-1] = (args.lfar, y_far)
    bot = np.column_stack([cl[:, 0], cl[:, 1] - h])
    top = np.column_stack([cl[:, 0], cl[:, 1] + h])
    bot[0], top[0] = (xc[0], yc[0]), (xc[-1], yc[-1])   # exact, unmodified TE corners
    nwake = len(cl)

    print(f"wake: {nwake} pts/side; last-region growth re-solved to {100*(r3_fit-1):.3f} % "
          f"(requested {args.growth[2]} %) so the wake ends exactly at x={args.lfar}")

    curve = [bot[::-1][:-1], np.column_stack([xc, yc]), top[1:]]
    if not args.open:
        close = np.array([[args.lfar, 0.0]])
        curve = [close] + curve + [close]
    curve = np.vstack(curve)
    seg = np.hypot(*np.diff(curve, axis=0).T)
    if seg.min() <= 1e-12:
        raise SystemExit("zero-length segment in assembled curve")
    shift = 1 if args.open else 0
    i_bot_te, i_top_te = nwake + 1 - shift, nwake + n_af - shift

    print_curve_stats(curve, n_af, nwake, args.open, args.x1, args.x2)
    write_airfoil_dat(args.out, header, curve[:, 0], curve[:, 1])
    write_meta(os.path.splitext(args.out)[0] + '.meta.json', dict(
        stage='step1', airfoil_source=args.airfoil, lfar=args.lfar, nwake=nwake,
        wake_type=args.wake_type,
        pull_length=(args.pull_length if args.wake_type == 'spline' else None),
        wake_grading=dict(growth_percent=list(args.growth), x1=args.x1, x2=args.x2,
                          last_region_growth_percent_fitted=100 * (r3_fit - 1)),
        n_af=n_af, h=h, d0=d0, n_total=len(curve), open=args.open,
        i_bot_te_raw=i_bot_te, i_top_te_raw=i_top_te, curve_path=args.out))
    print(f"assembled curve: {len(curve)} pts (bottomTE i={i_bot_te}, topTE i={i_top_te}, 1-based)")
    print(f"wrote {args.out}")


if __name__ == '__main__':
    main()
