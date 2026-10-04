#!/usr/bin/env python3
"""
Step 1 of the 30P-30N envelope pipeline: build the single open curve that
Construct2D (XCUT build: HYPR, CGRD, NWKE=0, BUFF) extrudes, plus a
.meta.json describing every piece of it. See PLAN.md next to this file.

Curve order (i increasing, the XCUT convention):

    lower downstream line (far end -> flap lower TE corner, corner excluded)
    flap lower surface      F[0 .. flap_lower_landing]
    lower flap bridge       (interior points only)
    main lower surface      M[57 .. main_lower_landing]
    lower slat bridge       (interior points only)
    slat bottom + outer     S[75 .. 199], S[0]  (hook tip -> LE -> trailing tip)
    upper slat bridge       (interior points only)
    main upper surface      M[main_upper_landing .. 219], M[0]
    upper flap bridge       (interior points only)
    flap upper surface      F[flap_upper_landing .. 241]
    upper downstream line   (flap upper TE corner excluded -> far end)

Every raw surface point on the curve is written bit-identical to the source
file. Surfaces are refined by insertion only (PLAN.md rule C2): points are
added on per-element splines (split at sharp corners) so the spacing starts
at the wall spacing y0 at every junction with a bridge or downstream line
and grows at no more than --wall-growth. Bridges and downstream lines start
from the same spacing, so the slot blocks resolve their walls and match the
Construct2D mesh across every interface (PLAN.md C3/C4). y0 is computed
exactly as Construct2D does from --ypls/--recd/--cfrc, which step 2 reuses. The two downstream ends stay open (no closure point). By default
(--far-end free) the downstream lines leave along the flap TE bisector and
turn onto a straight far section at --far-angle; this needs the open-cut
fix in this repo's Construct2D (hyperbolic_surface_grid.f90). With
--far-end mirror they instead turn back onto y=0 and end as exact mirror
images about y=0, which older XCUT builds require.

Usage:
    python3 step1_build_envelope.py --out ../sample_airfoils/30p30n_envelope_v1.dat
"""
import argparse
import os
import numpy as np

from env_common import (ElementSpline, raw_walk, refine_walk, wall_distance,
                        read_raw_element, write_curve_dat, write_meta,
                        refuse_existing, unit, heading, turn_deg, seg_lengths,
                        vinokur, choose_bridge_count, sample_by_arclength,
                        hermite_dense, arc_length, geometric_growth_rate,
                        polyline_crossings, point_to_polyline_distance,
                        max_neighbour_ratio)

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, '..', 'sample_airfoils', '30p30n_source')

# Fixed features of the raw files (0-based, after dropping the duplicate
# closing point of slat and main) -- PLAN.md "Geometry facts".
SLAT_TIP, SLAT_HOOK, SLAT_LE = 0, 75, 143
MAIN_TE, MAIN_COVE_TOP, MAIN_COVE_BOT, MAIN_LE = 0, 40, 57, 122
FLAP_LOWER_TE, FLAP_LE = 0, 97


def landing_for_angle(P, surf, lo, hi, theta):
    """Index k in [lo, hi] where a straight bridge P -> surf[k] meets the
    surface (central-difference tangent, increasing index) at the signed
    angle closest to theta (deg, ccw > 0, i.e. concave on the envelope).
    theta = 0 is the tangent landing. Returns (k, actual angle)."""
    best, best_k, best_a = None, None, None
    for k in range(lo, hi + 1):
        t = surf[(k + 1) % len(surf)] - surf[k - 1]
        a = turn_deg(surf[k] - P, t)
        if best is None or abs(a - theta) < best:
            best, best_k, best_a = abs(a - theta), k, a
    return best_k, best_a


def rotate_deg(v, deg):
    c, s_ = np.cos(np.radians(deg)), np.sin(np.radians(deg))
    return np.array([c * v[0] - s_ * v[1], s_ * v[0] + c * v[1]])


def straight_bridge(P0, P1, ds0, ds1, n):
    s = vinokur(np.hypot(*(P1 - P0)), n, ds0, ds1)
    t = s / s[-1]
    return P0[None, :] + t[:, None] * (P1 - P0)[None, :]


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--out', required=True, help='output curve .dat (new file)')
    ap.add_argument('--main-lower-landing', type=int, default=103,
                    help='main index where the lower slat bridge lands (default 103)')
    ap.add_argument('--flap-lower-landing', type=int, default=66,
                    help='flap index where the lower flap bridge lands (default 66)')
    ap.add_argument('--main-upper-landing', type=int, default=None,
                    help='main index for the upper slat bridge (default: found '
                         'from the landing angle)')
    ap.add_argument('--flap-upper-landing', type=int, default=None,
                    help='flap index for the upper flap bridge (default: found '
                         'from the landing angle)')
    ap.add_argument('--landing-angle', type=float, default=22.0,
                    help='angle (deg) at which every bridge meets the surface it '
                         'lands on. It is both the concave kink in the envelope '
                         'and the slot block\'s corner angle there. 0 = tangent '
                         '(smooth envelope, zero-angle block corner). Default 22: '
                         'block corners clear the 20 deg gate, and the '
                         'Construct2D mesh is unaffected up to about 23 deg '
                         '(PLAN.md, "Landing points")')
    for b in ('slat-upper', 'flap-upper', 'slat-lower', 'flap-lower'):
        ap.add_argument(f'--theta-{b}', type=float, default=None,
                        help=f'landing angle for the {b} bridge (default: --landing-angle)')
    ap.add_argument('--lower-arrival', choices=['chord', 'tangent'], default='chord',
                    help="lower bridges' end direction: along the bridge chord "
                         "(small convex kink at the hook tip / cove corner, no "
                         "S-turn) or tangent to the next surface (default chord)")
    ap.add_argument('--pull', type=float, default=1.0,
                    help='Hermite tangent length of the lower bridges, as a '
                         'fraction of the bridge chord (default 1.0)')
    ap.add_argument('--ypls', type=float, default=1.0,
                    help='y+ of the first cell (Construct2D YPLS; default 1)')
    ap.add_argument('--recd', type=float, default=1.7e6,
                    help='Reynolds number for y0 (Construct2D RECD; default 1.7e6, '
                         'the 30P-30N BANC case)')
    ap.add_argument('--cfrc', type=float, default=0.5,
                    help='reference length for y0 (Construct2D CFRC; default 0.5)')
    ap.add_argument('--wall-growth', type=float, default=1.15,
                    help='max spacing growth along surfaces away from a junction '
                         '(default 1.15)')
    ap.add_argument('--refine', type=int, default=1,
                    help='split every raw interval into at least this many parts '
                         '(default 1: raw density away from junctions)')
    ap.add_argument('--no-refine', action='store_true',
                    help='raw surface points only (the old coarse envelope)')
    ap.add_argument('--rmax', type=float, default=1.2,
                    help='max neighbour spacing ratio on bridges (default 1.2)')
    ap.add_argument('--dsmax', type=float, default=0.005,
                    help='max point spacing on bridges (default 0.005)')
    ap.add_argument('--match', choices=['slat', 'both', 'none'], default='slat',
                    help='which slots get equal point counts on their two bridges '
                         '(needed when one H-block spans the slot). Default slat: '
                         'the main-flap bridges differ ~10x in length, so that '
                         'slot will be split into sub-blocks (PLAN.md) and its '
                         'counts are left independent until that layout is fixed')
    ap.add_argument('--xfar', type=float, default=20.0,
                    help='x where the downstream lines end (default 20)')
    ap.add_argument('--far-end', choices=['free', 'mirror'], default='free',
                    help="free: straight far section at --far-angle, lines are "
                         "exact translates (needs this repo's open-cut "
                         "Construct2D fix). mirror: far ends at (xfar, -/+h), "
                         "for XCUT builds without the fix (default free)")
    ap.add_argument('--far-angle', type=float, default=0.0,
                    help='free mode: direction of the straight far section, '
                         'deg (default 0)')
    ap.add_argument('--turn-length', type=float, default=None,
                    help='length over which the downstream centerline turns from '
                         'the flap TE bisector to the far direction (default '
                         '1.0 free, 2.0 mirror)')
    ap.add_argument('--line-growth', type=float, default=1.15,
                    help='max geometric growth on the downstream lines (default 1.15)')
    ap.add_argument('--max-concave', type=float, default=30.0,
                    help='refuse if any single-point concave turn exceeds this (deg)')
    ap.add_argument('--clearance', type=float, default=0.2,
                    help='warn if the middle half of a bridge passes closer to a '
                         'block wall than this fraction of the slot gap (default '
                         '0.2). A warning, not a refusal: a tangent landing runs '
                         'close to its wall by design, and how close is part of '
                         'the landing trade-off (PLAN.md)')
    args = ap.parse_args()

    out = args.out
    meta_fn = os.path.splitext(out)[0] + '.meta.json'
    refuse_existing(out, meta_fn)

    S = read_raw_element(os.path.join(SRC, '30P-30N-Slat-raw.dat'), closed=True)
    M = read_raw_element(os.path.join(SRC, '30P-30N-Main-raw.dat'), closed=True)
    F = read_raw_element(os.path.join(SRC, '30P-30N-Flap-raw.dat'), closed=False)
    nS, nM, nF = len(S), len(M), len(F)
    FLAP_UPPER_TE = nF - 1

    # ---- landing points --------------------------------------------------
    ml_land = args.main_lower_landing
    fl_land = args.flap_lower_landing
    def th(name):
        v = getattr(args, 'theta_' + name.replace('-', '_'))
        return args.landing_angle if v is None else v
    theta = {b: th(b) for b in ('slat-upper', 'flap-upper', 'slat-lower', 'flap-lower')}

    # Upper landings: search the surface from just aft of the LE back to the
    # tangent point (a landing aft of it would cut into the element).
    mu_tan, _ = landing_for_angle(S[SLAT_TIP], M, MAIN_LE + 1, nM - 1, 0.0)
    fu_tan, _ = landing_for_angle(M[MAIN_TE], F, FLAP_LE + 1, nF - 2, 0.0)
    if args.main_upper_landing is not None:
        mu_land = args.main_upper_landing
    else:
        mu_land, _ = landing_for_angle(S[SLAT_TIP], M, MAIN_LE + 1, mu_tan,
                                       theta['slat-upper'])
    if args.flap_upper_landing is not None:
        fu_land = args.flap_upper_landing
    else:
        fu_land, _ = landing_for_angle(M[MAIN_TE], F, FLAP_LE + 1, fu_tan,
                                       theta['flap-upper'])
    print(f"Upper slat landing: main {mu_land} (tangent point {mu_tan}); "
          f"upper flap landing: flap {fu_land} (tangent point {fu_tan})")
    if not (MAIN_COVE_BOT < ml_land < MAIN_LE):
        raise SystemExit("--main-lower-landing must lie on the main lower surface")
    if not (FLAP_LOWER_TE < fl_land < FLAP_LE):
        raise SystemExit("--flap-lower-landing must lie on the flap lower surface")

    def sp(a, i, j):
        return float(np.hypot(*(a[j] - a[i])))

    # ---- refined surfaces ---------------------------------------------------
    y0 = wall_distance(args.ypls, args.recd, args.cfrc)
    esS = ElementSpline(S, True, [SLAT_TIP, SLAT_HOOK])
    esM = ElementSpline(M, True, [MAIN_TE, MAIN_COVE_TOP, MAIN_COVE_BOT])
    esF = ElementSpline(F, False, [])

    def surface(es, idx):
        if args.no_refine:
            return es.pts[idx], list(range(len(idx)))
        return refine_walk(es, idx, y0, args.wall_growth, args.refine)

    surf = {
        'flap_lower': ('flap', list(range(0, fl_land + 1)), esF),
        'main_lower': ('main', raw_walk(nM, MAIN_COVE_BOT, ml_land, True), esM),
        'slat_bottom_outer': ('slat', raw_walk(nS, SLAT_HOOK, SLAT_TIP, True), esS),
        'main_upper': ('main', raw_walk(nM, mu_land, MAIN_TE, True), esM),
        'flap_upper': ('flap', list(range(fu_land, nF)), esF),
    }
    spts, spos = {}, {}
    for name, (elem, idx, es) in surf.items():
        spts[name], spos[name] = surface(es, idx)
    print(f"Wall spacing y0 = {y0:.3e} (y+ {args.ypls}, Re {args.recd:g}, "
          f"cfrc {args.cfrc}); surfaces " + ("raw only" if args.no_refine else
          f"refined: " + ", ".join(f"{k} {len(surf[k][1])}->{len(v)}"
                                  for k, v in spts.items())))

    def first_ds(name):
        return float(np.hypot(*(spts[name][1] - spts[name][0])))

    def last_ds(name):
        return float(np.hypot(*(spts[name][-1] - spts[name][-2])))

    # ---- bridges -----------------------------------------------------------
    # Upper bridges: straight, landing at the chosen angle.
    up_slat = dict(P0=S[SLAT_TIP], P1=M[mu_land],
                   ds0=last_ds('slat_bottom_outer'), ds1=first_ds('main_upper'))
    up_flap = dict(P0=M[MAIN_TE], P1=F[fu_land],
                   ds0=last_ds('main_upper'), ds1=first_ds('flap_upper'))

    # Lower bridges: Hermite, leaving the surface they start on at the chosen
    # angle (rotated ccw from the surface tangent, i.e. away from the nose
    # the slot block wraps around).
    def lower_dense(P0, t0, P1, t_next):
        chord = P1 - P0
        t1 = chord if args.lower_arrival == 'chord' else t_next
        L = np.hypot(*chord)
        return hermite_dense(P0, P1, t0, t1, args.pull * L, args.pull * L)

    lo_slat_dense = lower_dense(M[ml_land],
                                rotate_deg(M[ml_land] - M[ml_land - 1], theta['slat-lower']),
                                S[SLAT_HOOK], S[SLAT_HOOK + 1] - S[SLAT_HOOK])
    lo_flap_dense = lower_dense(F[fl_land],
                                rotate_deg(F[fl_land] - F[fl_land - 1], theta['flap-lower']),
                                M[MAIN_COVE_BOT], M[MAIN_COVE_BOT + 1] - M[MAIN_COVE_BOT])
    lo_slat = dict(L=arc_length(lo_slat_dense), ds0=last_ds('main_lower'),
                   ds1=first_ds('slat_bottom_outer'))
    lo_flap = dict(L=arc_length(lo_flap_dense), ds0=last_ds('flap_lower'),
                   ds1=first_ds('main_lower'))
    for b in (up_slat, up_flap):
        b['L'] = float(np.hypot(*(b['P1'] - b['P0'])))

    for b in (up_slat, up_flap, lo_slat, lo_flap):
        b['n'] = choose_bridge_count(b['L'], b['ds0'], b['ds1'], args.rmax, args.dsmax)
    pairs = {'slat': [(up_slat, lo_slat)], 'none': [],
             'both': [(up_slat, lo_slat), (up_flap, lo_flap)]}[args.match]
    if pairs:
        # Raising one bridge's count to match its partner changes its
        # distribution, so search upward until both pass the limits at the
        # same count.
        for a, b in pairs:
            n = max(a['n'], b['n'])
            while True:
                na = choose_bridge_count(a['L'], a['ds0'], a['ds1'], args.rmax,
                                         args.dsmax, nmin=n)
                nb = choose_bridge_count(b['L'], b['ds0'], b['ds1'], args.rmax,
                                         args.dsmax, nmin=n)
                if na == nb == n:
                    break
                n = max(na, nb)
            a['n'] = b['n'] = n

    up_slat_pts = straight_bridge(up_slat['P0'], up_slat['P1'],
                                  up_slat['ds0'], up_slat['ds1'], up_slat['n'])
    up_flap_pts = straight_bridge(up_flap['P0'], up_flap['P1'],
                                  up_flap['ds0'], up_flap['ds1'], up_flap['n'])
    lo_slat_pts = sample_by_arclength(lo_slat_dense,
                                      vinokur(lo_slat['L'], lo_slat['n'],
                                              lo_slat['ds0'], lo_slat['ds1']))
    lo_flap_pts = sample_by_arclength(lo_flap_dense,
                                      vinokur(lo_flap['L'], lo_flap['n'],
                                              lo_flap['ds0'], lo_flap['ds1']))

    # ---- downstream lines --------------------------------------------------
    # One centerline from the flap TE midpoint along the TE bisector, a
    # Hermite turn onto the far direction, then straight to x = xfar. Both
    # lines share the centerline's point distribution exactly.
    #
    # free:   lines = centerline -/+ g/2 (exact translates; g = TE gap
    #         vector). Needs the open-cut fix in this repo's Construct2D,
    #         which marches i=1 and i=imax independently.
    # mirror: the far direction is +x along y=0, and the offset rotates
    #         smoothly from g/2 to (0, |g|/2), so the far ends are exact
    #         mirror images (xfar, -/+h). Without the open-cut fix the C-grid
    #         marcher sets column imax to the mirror of column 1 about y=0 at
    #         every level, so only this shape meshes cleanly there.
    g = F[FLAP_UPPER_TE] - F[FLAP_LOWER_TE]
    h = 0.5 * float(np.hypot(*g))
    Pm = 0.5 * (F[FLAP_LOWER_TE] + F[FLAP_UPPER_TE])
    t0 = unit(unit(F[FLAP_LOWER_TE] - F[1]) + unit(F[FLAP_UPPER_TE] - F[FLAP_UPPER_TE - 1]))
    mirror = args.far_end == 'mirror'
    Lturn = args.turn_length if args.turn_length is not None else (2.0 if mirror else 1.0)
    if mirror:
        efar = np.array([1.0, 0.0])
        P1 = np.array([Pm[0] + Lturn, 0.0])
    else:
        efar = np.array([np.cos(np.radians(args.far_angle)),
                         np.sin(np.radians(args.far_angle))])
        P1 = Pm + 0.5 * Lturn * (t0 + efar)
    if args.xfar <= P1[0] or efar[0] <= 0:
        raise SystemExit("--xfar must lie downstream of the turn section")
    Lt = float(np.hypot(*(P1 - Pm)))
    turn = hermite_dense(Pm, P1, t0, efar, Lt, Lt)
    Lstraight = (args.xfar - P1[0]) / efar[0]
    straight = P1[None, :] + np.linspace(0.0, Lstraight, 2000)[:, None] * efar[None, :]
    center_dense = np.vstack([turn, straight[1:]])
    Lturn_arc = arc_length(turn)
    Lc = arc_length(center_dense)
    d0 = min(first_ds('flap_lower'), last_ds('flap_upper'))
    nline = 2
    while geometric_growth_rate(Lc, d0, nline + 1) > args.line_growth:
        nline += 1
    rline = geometric_growth_rate(Lc, d0, nline + 1)
    s = np.concatenate([[0.0], np.cumsum(d0 * rline ** np.arange(nline))])
    s[-1] = Lc
    center = sample_by_arclength(center_dense, s)
    if mirror:
        in_straight = s >= Lturn_arc
        center[in_straight, 1] = 0.0          # exactly on y=0 (interpolation noise)
        center[-1] = [args.xfar, 0.0]
        w = np.clip(s / Lturn_arc, 0.0, 1.0)
        w = w * w * (3.0 - 2.0 * w)            # smoothstep
        odir = ((1.0 - w)[:, None] * unit(g)[None, :]
                + w[:, None] * np.array([0.0, 1.0])[None, :])
        odir /= np.hypot(*odir.T)[:, None]
        off = h * odir
        off[in_straight] = [0.0, h]
    else:
        off = np.tile(0.5 * g, (len(center), 1))
    lower_line = center - off   # starts (to round-off) at the flap lower TE corner
    upper_line = center + off   # starts at the flap upper TE corner
    if mirror:
        lower_line[-1] = [args.xfar, -h]
        upper_line[-1] = [args.xfar, h]
        if lower_line[-1][0] != upper_line[-1][0] or lower_line[-1][1] != -upper_line[-1][1]:
            raise SystemExit("internal error: far ends are not exact mirror images")
    dc = np.diff(center, axis=0)
    gap_n = np.abs(dc[:, 0] * (upper_line - lower_line)[1:, 1]
                   - dc[:, 1] * (upper_line - lower_line)[1:, 0]) / seg_lengths(center)
    print(f"Downstream lines ({args.far_end}): {nline} intervals each, growth "
          f"{rline:.4f}, first step {d0:.3e}; centerline y range "
          f"{center[:, 1].min():.4f}..{center[:, 1].max():.4f}; normal gap "
          f"{gap_n.min():.5f}..{gap_n.max():.5f}; far ends "
          f"({lower_line[-1][0]:.4f}, {lower_line[-1][1]:.4f}) / "
          f"({upper_line[-1][0]:.4f}, {upper_line[-1][1]:.4f})")

    # ---- assemble ---------------------------------------------------------
    pieces = []   # (name, points, source element, source index list or None)

    def add(name, pts, elem=None, idx=None):
        pieces.append((name, np.asarray(pts, float), elem, idx))

    add('lower_line', lower_line[1:][::-1])
    def add_surf(name):
        add(name, spts[name], surf[name][0], surf[name][1])

    add_surf('flap_lower')
    add('bridge_flap_lower', lo_flap_pts[1:-1])
    add_surf('main_lower')
    add('bridge_slat_lower', lo_slat_pts[1:-1])
    add_surf('slat_bottom_outer')
    add('bridge_slat_upper', up_slat_pts[1:-1])
    add_surf('main_upper')
    add('bridge_flap_upper', up_flap_pts[1:-1])
    add_surf('flap_upper')
    add('upper_line', upper_line[1:])

    curve = np.vstack([p[1] for p in pieces])
    ranges, start = {}, 1
    for name, pts, _, _ in pieces:
        ranges[name] = [start, start + len(pts) - 1]
        start += len(pts)
    ntot = len(curve)

    # Bit-identical check of every raw surface point.
    for name, pts, elem, idx in pieces:
        if elem is None:
            continue
        src = {'slat': S, 'main': M, 'flap': F}[elem]
        if not np.array_equal(pts[spos[name]], src[idx]):
            raise SystemExit(f"internal error: {name} is not bit-identical to source")

    # Bridge index ranges including their two surface endpoints (these are
    # the interfaces shared with the slot blocks).
    def with_ends(name):
        a, b = ranges[name]
        return [a - 1, b + 1]

    # ---- checks -------------------------------------------------------------
    failures, warnings = [], []

    hits = polyline_crossings(curve)
    if hits:
        failures.append(f"curve crosses itself at segment pairs {hits[:5]}")

    walls = {
        'slat_cove': S[SLAT_TIP:SLAT_HOOK + 1],
        'main_nose': M[ml_land:mu_land + 1],
        'main_cove': M[MAIN_TE:MAIN_COVE_BOT + 1],
        'flap_nose': F[fl_land:fu_land + 1],
    }
    for wname, w in walls.items():
        h = polyline_crossings(curve, w)
        if h:
            failures.append(f"curve cuts through block wall {wname} ({len(h)} crossings)")

    slot_gap = {
        'slat': float(np.min(np.hypot(*(M - S[SLAT_TIP]).T))),
        'flap': float(np.min(np.hypot(*(F - M[MAIN_TE]).T))),
    }
    clear = {}
    for bname, pts, slot, wnames in (
            ('bridge_slat_upper', up_slat_pts, 'slat', ('slat_cove', 'main_nose')),
            ('bridge_slat_lower', lo_slat_pts, 'slat', ('slat_cove', 'main_nose')),
            ('bridge_flap_upper', up_flap_pts, 'flap', ('main_cove', 'flap_nose')),
            ('bridge_flap_lower', lo_flap_pts, 'flap', ('main_cove', 'flap_nose'))):
        arc = np.concatenate([[0.0], np.cumsum(seg_lengths(pts))])
        mid = pts[(arc > 0.25 * arc[-1]) & (arc < 0.75 * arc[-1])]
        d = min(point_to_polyline_distance(mid, walls[w]).min() for w in wnames)
        clear[bname] = float(d)
        if d < args.clearance * slot_gap[slot]:
            warnings.append(f"{bname} passes {d:.4f} from a block wall "
                            f"(< {args.clearance} x slot gap {slot_gap[slot]:.4f})")

    # Turn angle at every interior point. The curve runs clockwise around the
    # body (body on the right), so a ccw (positive) turn is concave.
    tang = np.diff(curve, axis=0)
    turns = np.array([turn_deg(tang[k - 1], tang[k]) for k in range(1, len(tang))])
    worst = np.argsort(-turns)[:3]
    if turns.max() > args.max_concave:
        failures.append(f"concave turn {turns.max():.1f} deg at curve point "
                        f"{int(np.argmax(turns)) + 2}")

    ds = seg_lengths(curve)
    ratios = np.maximum(ds[1:] / ds[:-1], ds[:-1] / ds[1:])

    def turn_at(i1):          # 1-based curve index -> turn there
        return float(turns[i1 - 2])

    junctions = {
        'flap_lower_landing': ranges['flap_lower'][1],
        'cove_bottom_corner': ranges['main_lower'][0],
        'main_lower_landing': ranges['main_lower'][1],
        'hook_tip': ranges['slat_bottom_outer'][0],
        'slat_tip': ranges['slat_bottom_outer'][1],
        'main_upper_landing': ranges['main_upper'][0],
        'main_te': ranges['main_upper'][1],
        'flap_upper_landing': ranges['flap_upper'][0],
    }

    # ---- report -------------------------------------------------------------
    print(f"\nCurve: {ntot} points")
    for name in ranges:
        print(f"  {name:20s} {ranges[name][0]:5d} .. {ranges[name][1]:5d}")
    print("\nBridges (intervals, length, end spacings, clearance to block walls):")
    for name, b in (('bridge_slat_upper', up_slat), ('bridge_slat_lower', lo_slat),
                    ('bridge_flap_upper', up_flap), ('bridge_flap_lower', lo_flap)):
        print(f"  {name:20s} n={b['n']:4d}  L={b['L']:.4f}  ds0={b['ds0']:.2e} "
              f"ds1={b['ds1']:.2e}  clearance={clear[name]:.4f}")
    print("\nTurn at each junction (deg, + = concave, - = convex):")
    for k, i1 in junctions.items():
        print(f"  {k:20s} point {i1:5d}  {turn_at(i1):+7.2f}")
    print(f"\nLargest concave turns: " +
          ", ".join(f"{turns[w]:.1f} deg at point {w + 2}" for w in worst))
    print(f"Max neighbour spacing ratio: {ratios.max():.3f} at point "
          f"{int(np.argmax(ratios)) + 2}")

    for w_ in warnings:
        print("WARNING: " + w_)

    if failures:
        print("\nREFUSING to write the curve:")
        for f_ in failures:
            print("  - " + f_)
        raise SystemExit(1)

    # ---- write --------------------------------------------------------------
    write_curve_dat(out, '30P-30N envelope (step 1)', curve)
    meta = {
        'stage': 'step1',
        'curve_path': out,
        'n_total': ntot,
        'source_dir': os.path.relpath(SRC, HERE),
        'params': vars(args),
        'landings': {'main_upper': mu_land, 'flap_upper': fu_land,
                     'main_lower': ml_land, 'flap_lower': fl_land},
        'landing_angles_requested': theta,
        'bridge_intervals': {k: b['n'] for k, b in (
            ('slat_upper', up_slat), ('slat_lower', lo_slat),
            ('flap_upper', up_flap), ('flap_lower', lo_flap))},
        'downstream_lines': {'intervals': nline, 'growth': rline, 'first_step': d0,
                             'far_end': args.far_end, 'turn_length': Lturn,
                             'te_gap_vector': g.tolist(),
                             'far_ends': [lower_line[-1].tolist(),
                                          upper_line[-1].tolist()]},
        'block_walls_source': {
            'slat_cove': ['slat', SLAT_TIP, SLAT_HOOK],
            'main_nose': ['main', ml_land, mu_land],
            'main_cove': ['main', MAIN_TE, MAIN_COVE_BOT],
            'flap_nose': ['flap', fl_land, fu_land],
        },
        'surface_source': {name: {'element': elem, 'source_indices': idx,
                                  'curve_positions': [ranges[name][0] + p
                                                      for p in spos[name]]}
                           for name, _, elem, idx in pieces if elem is not None},
        'refinement': {'y0': y0, 'ypls': args.ypls, 'recd': args.recd,
                       'cfrc': args.cfrc, 'wall_growth': args.wall_growth,
                       'refine': args.refine, 'no_refine': args.no_refine,
                       'spline_corners': {'slat': [SLAT_TIP, SLAT_HOOK],
                                          'main': [MAIN_TE, MAIN_COVE_TOP, MAIN_COVE_BOT],
                                          'flap': []}},
        'junction_turns_deg': {k: turn_at(i1) for k, i1 in junctions.items()},
        'bridge_clearance_mid_half': clear,
        'warnings': warnings,
        'slot_gap': slot_gap,
        # 1-based inclusive ranges, read by visualize_curve.py
        **{f'i_{name}': r for name, r in ranges.items()},
        **{f'i_if_{name}': with_ends(name) for name in
           ('bridge_slat_upper', 'bridge_slat_lower',
            'bridge_flap_upper', 'bridge_flap_lower')},
        'i_hook_tip': junctions['hook_tip'],
        'i_slat_le': ranges['slat_bottom_outer'][0]
                     + spos['slat_bottom_outer'][SLAT_LE - SLAT_HOOK],
        'i_slat_tip': junctions['slat_tip'],
        'i_main_te': junctions['main_te'],
        'i_cove_bottom_corner': junctions['cove_bottom_corner'],
    }
    write_meta(meta_fn, meta)
    print(f"\nWrote {out}\nWrote {meta_fn}")


if __name__ == '__main__':
    main()
