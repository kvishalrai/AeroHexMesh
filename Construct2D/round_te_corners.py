#!/usr/bin/env python3
"""
Locally round the two sharp trailing-edge corners of an airfoil .dat file
that is "closed" via a phantom TE-midpoint (i.e. the file starts and ends
at the same point, e.g. (1.0, 0.0), with the two real TE corners at
indices 2 and n-1 (1-based), or 1 and n-2 (0-based)).

This does NOT touch the rest of the airfoil, does NOT change the overall
chord/extents (the rounding cuts a small notch INTO the existing corner
using a quadratic Bezier blend that is tangent to both adjacent segments
and bounded by the original corner triangle -- it can never poke out past
the original geometry) and does NOT require opening the TE gap, does not
run any re-paneling/rescaling. It just replaces each sharp corner vertex
with a handful of points along a smooth, tangent-continuous local blend.

Usage: round_te_corners.py input.dat output.dat [frac] [npts]
  frac : fraction (0-1) of the shorter adjacent segment length used as
         the rounding "radius" at each corner (default 0.85 -- a good
         skew/growth-ratio balance found for the OAT15 case; try
         0.6-0.95 to trade off skew angle vs. local cell stretching)
  npts : number of points used to represent each rounded corner,
         including its two blend endpoints (default 9)
"""
import sys

def read_dat(path):
    with open(path) as f:
        lines = f.readlines()
    header = lines[0]
    pts = []
    for line in lines[1:]:
        line = line.strip()
        if not line:
            continue
        x, y = line.split()[:2]
        pts.append((float(x), float(y)))
    return header, pts

def write_dat(path, header, pts):
    with open(path, "w") as f:
        f.write(header if header.endswith("\n") else header + "\n")
        for x, y in pts:
            f.write(f"    {x:.10e}     {y:.10e}\n")

def dist(a, b):
    return ((a[0]-b[0])**2 + (a[1]-b[1])**2) ** 0.5

def round_corner(prev_pt, corner, next_pt, frac, npts):
    """Return a list of npts points (including endpoints) replacing `corner`,
    blending tangentially into the prev_pt->corner and corner->next_pt segments."""
    d1 = dist(prev_pt, corner)
    d2 = dist(corner, next_pt)
    r = frac * min(d1, d2)

    def along(a, b, dist_from_a):
        L = dist(a, b)
        t = dist_from_a / L
        return (a[0] + t*(b[0]-a[0]), a[1] + t*(b[1]-a[1]))

    P1 = along(corner, prev_pt, r)   # point on the prev-corner segment, r from corner
    P2 = along(corner, next_pt, r)   # point on the corner-next segment, r from corner

    pts = []
    for i in range(npts):
        t = i / (npts - 1)
        x = (1-t)**2 * P1[0] + 2*(1-t)*t*corner[0] + t**2 * P2[0]
        y = (1-t)**2 * P1[1] + 2*(1-t)*t*corner[1] + t**2 * P2[1]
        pts.append((x, y))
    return pts

def main():
    if len(sys.argv) < 3:
        print(__doc__)
        sys.exit(1)
    infile, outfile = sys.argv[1], sys.argv[2]
    frac = float(sys.argv[3]) if len(sys.argv) > 3 else 0.85
    npts = int(sys.argv[4]) if len(sys.argv) > 4 else 9

    header, pts = read_dat(infile)
    n = len(pts)
    assert pts[0] == pts[-1], "Expected a closed loop (first point == last point)"

    # corner indices (0-based): 1 (just after the phantom tip) and n-2 (just before it)
    lower_corner_idx = 1
    upper_corner_idx = n - 2

    lower_blend = round_corner(pts[0], pts[lower_corner_idx], pts[lower_corner_idx+1], frac, npts)
    upper_blend = round_corner(pts[upper_corner_idx-1], pts[upper_corner_idx], pts[-1], frac, npts)

    new_pts = []
    new_pts.append(pts[0])                     # phantom tip (unchanged)
    new_pts.extend(lower_blend)                # replaces old lower corner point
    new_pts.extend(pts[lower_corner_idx+2:upper_corner_idx-1])  # unchanged middle surface
    new_pts.extend(upper_blend)                # replaces old upper corner point
    new_pts.append(pts[-1])                    # phantom tip again (== pts[0])

    write_dat(outfile, header, new_pts)
    print(f"{infile}: {n} points -> {outfile}: {len(new_pts)} points "
          f"(frac={frac}, npts/corner={npts})")

if __name__ == "__main__":
    main()
