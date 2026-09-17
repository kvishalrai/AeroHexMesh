#!/usr/bin/env python3
"""
Step 1 of the blunt-TE XCUT meshing pipeline: build the "whole curve"
.dat file that Construct2D loads for a CGRD + NWKE=0 + BUFF run.

Given the user's raw airfoil geometry, this:
  - extracts the true airfoil surface (bottomTE corner -> ... -> topTE
    corner), auto-stripping any pre-existing TE-closure artifact the
    source file may carry (see xcut_common.strip pipeline for details);
  - symmetrizes the TE half-gap (h) and the local TE point spacing (d0)
    between top and bottom, so the two new wake lines are mirror images
    of each other;
  - appends straight, constant-y wake lines (geometric point growth,
    first step = the airfoil's own local TE spacing) from x=1 out to
    x=Lfar on both sides;
  - closes the whole curve with a single shared point at (Lfar, 0.0)
    as both the first and last point of the array, which is what makes
    Construct2D see a *closed* curve (tegap=False) and unlocks BUFF.

Usage:
    python3 step1_build_curve.py \\
        --airfoil sample_airfoils/oat15.dat \\
        --lfar 20.0 --nwake 50 \\
        --out sample_airfoils/oat15_xcut.dat

`--lfar` is an ABSOLUTE x-coordinate (chord normalized to 1.0, matching
this project's established convention -- e.g. --lfar 20.0 means the
wake extends to 20 chords downstream of the LE, i.e. 19 chords beyond
the TE). If you actually meant "length beyond the TE", pass 1.0+that
length instead.

`--nwake` is the number of points PER SIDE (bottom and top each get
this many, including the TE corner point itself as the first one).

Writes the .dat curve plus a JSON metadata sidecar (<out>.meta.json)
that later pipeline steps consume -- never guesses indices from
geometry again once this is recorded.
"""
import argparse
import os
import numpy as np

from xcut_common import (
    read_airfoil_dat, write_airfoil_dat, strip_te_closure_artifact,
    build_wake_line, build_wake_centerline_spline, write_meta,
)


def extract_airfoil_surface(x, y, tol=1e-6):
    """
    Strip any TE-closure artifact (x[0]==x[-1]) and any partial
    TE-face resolution points, returning the curve trimmed down to
    exactly [bottomTE corner .. topTE corner], plus the two corner
    y-values. bottomTE/topTE corner = the extremal-|y| point among
    the run of points sitting at x == max(x) at each end of the
    array.
    """
    x, y, stripped = strip_te_closure_artifact(x, y)
    if stripped:
        print("  stripped a shared (x[0]==x[-1]) TE-closure artifact point")

    x_te = x.max()

    # leading run: points at x==x_te starting from index 0
    i = 0
    while i < len(x) and abs(x[i] - x_te) <= tol:
        i += 1
    lead_run = slice(0, i)
    if i == 0:
        raise ValueError("airfoil curve does not start at the TE (x=max(x))")
    i_bot = lead_run.start + int(np.argmax(np.abs(y[lead_run])))

    # trailing run: points at x==x_te ending at the last index
    j = len(x)
    k = j
    while k > 0 and abs(x[k - 1] - x_te) <= tol:
        k -= 1
    trail_run = slice(k, j)
    if k == j:
        raise ValueError("airfoil curve does not end at the TE (x=max(x))")
    i_top = trail_run.start + int(np.argmax(np.abs(y[trail_run])))

    n_stripped_lead = i_bot - lead_run.start
    n_stripped_trail = trail_run.stop - 1 - i_top
    if n_stripped_lead or n_stripped_trail:
        print(f"  discarded {n_stripped_lead} leading / {n_stripped_trail} "
              f"trailing TE-face point(s) inside of the true corners "
              f"(these are superseded by the new wake-line + box treatment)")

    xc = x[i_bot:i_top + 1]
    yc = y[i_bot:i_top + 1]
    return xc, yc


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--airfoil', required=True, help='input airfoil .dat (sample_airfoils/...)')
    ap.add_argument('--lfar', type=float, required=True, help='absolute x-coordinate of wake far end')
    ap.add_argument('--nwake', type=int, required=True, help='wake points per side (TE corner inclusive)')
    ap.add_argument('--wake-type', choices=['straight', 'spline'], default='straight',
                     help="'straight' (default): constant-y lines, matches every run "
                          "validated so far. 'spline': a single smooth Hermite centerline "
                          "from the TE-bisector direction to horizontal far downstream, "
                          "offset by +-h for the two sides -- exactly parallel and "
                          "equidistant by construction, since both sides share the same "
                          "underlying curve shifted by a constant.")
    ap.add_argument('--pull-length', type=float, default=1.0,
                     help="(--wake-type spline only) how far the TE-bisector tangent's "
                          "influence reaches, as an ABSOLUTE length in chords (NOT a "
                          "fraction of the wake length) -- ~1.0 gives a gentle sag with "
                          "no overshoot; increase for a longer, more gradual turn, "
                          "decrease for a tighter one hugging the TE angle")
    ap.add_argument('--out', required=True, help='output .dat path')
    args = ap.parse_args()

    header, x_raw, y_raw = read_airfoil_dat(args.airfoil)
    print(f"read {len(x_raw)} points from {args.airfoil}")

    xc, yc = extract_airfoil_surface(x_raw, y_raw)
    n_af = len(xc)
    print(f"airfoil surface (corner-to-corner): {n_af} points")

    # HARD REQUIREMENT: the airfoil coordinates (including the TE corner
    # points themselves) are never modified -- xc/yc are used exactly as
    # extracted from here on, with no overwriting. Equidistant/parallel
    # wake lines are instead achieved by centering the construction on the
    # TRUE corner midpoint (y_mid, generally not exactly 0) rather than by
    # nudging the corners onto a forced +-h -- see below.
    y_bot_raw, y_top_raw = yc[0], yc[-1]
    y_mid = 0.5 * (y_bot_raw + y_top_raw)
    h = 0.5 * (y_top_raw - y_bot_raw)
    print(f"TE corners (unmodified): bottom y={y_bot_raw:.7f}  top y={y_top_raw:.7f}"
          f"  -> half-gap h={h:.7f} (about centerline y_mid={y_mid:.7f})")

    d0_bot = np.hypot(xc[1] - xc[0], yc[1] - yc[0])
    d0_top = np.hypot(xc[-1] - xc[-2], yc[-1] - yc[-2])
    d0 = 0.5 * (d0_bot + d0_top)
    print(f"local TE point spacing: bottom={d0_bot:.6g}  top={d0_top:.6g}"
          f"  -> shared d0={d0:.6g}")

    if args.lfar <= xc[0]:
        raise SystemExit(f"--lfar ({args.lfar}) must be > TE x-coordinate ({xc[0]})")

    if abs(xc[0] - xc[-1]) > 1e-9:
        print(f"NOTE: bottomTE x ({xc[0]:.9f}) and topTE x ({xc[-1]:.9f}) differ -- "
              f"wake lines attach to each corner's own x exactly, unchanged")

    if args.wake_type == 'straight':
        # attach directly to each corner's own unmodified (x, y) -- gap is
        # y_top_raw - y_bot_raw everywhere by construction (both constant-y),
        # exactly equidistant with zero coordinate changes needed
        bot_wake, r_bot, L = build_wake_line(xc[0], y_bot_raw, d0, args.nwake, args.lfar)
        top_wake, r_top, _ = build_wake_line(xc[-1], y_top_raw, d0, args.nwake, args.lfar)
        print(f"wake line (straight): {args.nwake} pts/side, "
              f"growth ratio={r_bot:.5f}, length={L:.4f}")
    else:
        # TE-bisector direction: average of the two surface tangents as they
        # exit the TE corners (pointing downstream), from the ACTUAL
        # (possibly asymmetric) local geometry -- using one shared bisector
        # for both sides is what guarantees the two offset wake lines stay
        # exactly parallel/equidistant, rather than fitting two splines
        # independently and hoping they stay close.
        t_bot = np.array([xc[0] - xc[1], yc[0] - yc[1]])
        t_top = np.array([xc[-1] - xc[-2], yc[-1] - yc[-2]])
        bisector = t_bot / np.linalg.norm(t_bot) + t_top / np.linalg.norm(t_top)
        if bisector[0] < 0:
            raise SystemExit("computed TE-bisector direction points upstream -- "
                              "check the airfoil curve orientation")
        angle_deg = np.degrees(np.arctan2(bisector[1], bisector[0]))
        print(f"TE-bisector direction: {angle_deg:.3f} deg from +x")

        # centerline starts at the TRUE corner midpoint y_mid (not forced to
        # 0) and ends at y=0 far downstream; offsetting by +-h from this
        # centerline reproduces y_bot_raw/y_top_raw EXACTLY at the TE end
        # (since y_mid-h == y_bot_raw and y_mid+h == y_top_raw algebraically)
        # -- no coordinate nudging needed, and the two sides stay exactly
        # parallel/equidistant since they share one offset curve.
        centerline, r_bot, L = build_wake_centerline_spline(
            xc[0], y_mid, bisector, args.lfar, 0.0, (1.0, 0.0),
            d0, args.nwake, pull_length=args.pull_length)
        bot_wake = np.column_stack([centerline[:, 0], centerline[:, 1] - h])
        top_wake = np.column_stack([centerline[:, 0], centerline[:, 1] + h])
        assert abs(bot_wake[0, 1] - y_bot_raw) < 1e-12
        assert abs(top_wake[0, 1] - y_top_raw) < 1e-12
        print(f"wake line (spline): {args.nwake} pts/side, pull_length={args.pull_length}, "
              f"growth ratio={r_bot:.5f}, arc length={L:.4f}, "
              f"peak centerline offset={np.max(np.abs(centerline[:,1])):.5f}")
        gap = top_wake[:, 1] - bot_wake[:, 1]
        print(f"gap check: min={gap.min():.8f} max={gap.max():.8f} "
              f"(constant by construction, should equal 2h={2*h:.8f})")
        print(f"TE-attachment check: bottom wake starts at y={bot_wake[0,1]:.10f} "
              f"(corner y={y_bot_raw:.10f}), top wake starts at y={top_wake[0,1]:.10f} "
              f"(corner y={y_top_raw:.10f}) -- exact match required and asserted above")

    assert np.allclose(bot_wake[:, 0], top_wake[:, 0])

    close_pt = np.array([[args.lfar, 0.0]])

    # bottom wake ordered far->near (matches this project's established
    # curve orientation); EXCLUDE its near-TE point (shared with xc[0])
    bot_far_to_near = bot_wake[::-1][:-1]
    # top wake ordered near->far; EXCLUDE its near-TE point (shared with xc[-1])
    top_near_to_far = top_wake[1:]

    curve = np.vstack([close_pt, bot_far_to_near,
                        np.column_stack([xc, yc]),
                        top_near_to_far, close_pt])

    n_total = len(curve)
    i_bot_te = args.nwake + 1          # 1-based index of bottomTE in the raw curve
    i_top_te = args.nwake + n_af       # 1-based index of topTE in the raw curve
    print(f"assembled curve: {n_total} points "
          f"(bottomTE at i={i_bot_te}, topTE at i={i_top_te}, 1-based)")

    write_airfoil_dat(args.out, header, curve[:, 0], curve[:, 1])
    meta_path = os.path.splitext(args.out)[0] + '.meta.json'
    write_meta(meta_path, dict(
        stage='step1',
        airfoil_source=args.airfoil,
        lfar=args.lfar,
        nwake=args.nwake,
        wake_type=args.wake_type,
        pull_length=(args.pull_length if args.wake_type == 'spline' else None),
        n_af=n_af,
        h=h,
        d0=d0,
        n_total=n_total,
        i_bot_te_raw=i_bot_te,
        i_top_te_raw=i_top_te,
        curve_path=args.out,
    ))
    print(f"wrote {args.out}")
    print(f"wrote {meta_path}")
    print()
    print("Next: load this curve in Construct2D and run CGRD + NWKE=0 + BUFF "
          "to generate the raw p3d/nmf (step 2, done by you interactively).")


if __name__ == '__main__':
    main()
