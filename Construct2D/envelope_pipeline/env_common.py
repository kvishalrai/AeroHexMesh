"""
Shared helpers for the 30P-30N envelope pipeline (see PLAN.md next to this
file).

Reuses the xcut pipeline's helpers where they fit (meta I/O, growth-rate
solver, Hermite evaluation, p3d reading, quality stats) by importing
xcut_common from ../xcut_pipeline, which is left untouched. Writers here use
full float64 precision ("%.17g" for curves, "%25.16E" for p3d) so surface
points survive every round trip bit-for-bit (PLAN.md, constraint C5).

All files this pipeline writes are NEW files: every writer refuses to
overwrite an existing path.
"""
import os
import sys
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                '..', 'xcut_pipeline'))
from xcut_common import (geometric_growth_rate, _hermite_points,  # noqa: E402,F401
                         write_meta, read_meta, read_p3d,
                         compute_quality_stats)


# ---------------------------------------------------------------------
# File safety
# ---------------------------------------------------------------------

def refuse_existing(*paths):
    """Stop before writing anything if any output path already exists."""
    existing = [p for p in paths if p and os.path.exists(p)]
    if existing:
        raise SystemExit("Refusing to overwrite existing file(s):\n  "
                         + "\n  ".join(existing))


# ---------------------------------------------------------------------
# Raw element I/O
# ---------------------------------------------------------------------

def read_raw_element(fn, closed):
    """Plain 'x y' per line. For a closed element (slat, main) the file
    repeats its first point at the end; that duplicate is dropped, so
    indices match PLAN.md's geometry table."""
    a = np.loadtxt(fn)
    if closed:
        if not np.array_equal(a[0], a[-1]):
            raise ValueError(f"{fn}: expected a repeated closing point")
        a = a[:-1]
    return a


def write_curve_dat(fn, header, pts):
    """Construct2D (XFoil labeled) format: header line, then 'x y' rows at
    full precision."""
    with open(fn, 'w') as f:
        f.write(header + '\n')
        for x, y in pts:
            f.write(f"{x:.17g} {y:.17g}\n")


def read_curve_dat(fn):
    with open(fn) as f:
        f.readline()
        return np.loadtxt(f)


def write_p3d_full(fn, blocks):
    """Plot3D 2D ASCII, one value per line at full precision. Single block
    is written without the block-count line (Construct2D's own style)."""
    with open(fn, 'w') as f:
        if len(blocks) > 1:
            f.write(f"{len(blocks)}\n")
        for x, _ in blocks:
            f.write(f"{x.shape[0]} {x.shape[1]}\n")
        for x, y in blocks:
            for arr in (x, y):
                for v in arr.flatten(order='F'):
                    f.write(f"{v:25.16E}\n")


# ---------------------------------------------------------------------
# Geometry helpers
# ---------------------------------------------------------------------

def unit(v):
    v = np.asarray(v, dtype=float)
    return v / np.linalg.norm(v)


def heading(v):
    """Direction of a vector in degrees, (-180, 180]."""
    return float(np.degrees(np.arctan2(v[1], v[0])))


def turn_deg(u, v):
    """Signed turn from direction u to direction v in degrees (ccw > 0)."""
    return float(np.degrees(np.arctan2(u[0] * v[1] - u[1] * v[0],
                                       u[0] * v[0] + u[1] * v[1])))


def seg_lengths(pts):
    return np.hypot(*np.diff(pts, axis=0).T)


def vinokur(L, n, ds0, ds1):
    """n intervals on [0, L], first spacing ds0, last ds1 (Vinokur 1983
    two-sided stretching). Returns n+1 arc-length stations."""
    s0, s1 = ds0 / L, ds1 / L
    A = np.sqrt(s1 / s0)
    B = 1.0 / (n * np.sqrt(s0 * s1))
    xi = np.arange(n + 1) / n
    if abs(B - 1.0) < 1e-9:
        u = xi
    elif B > 1.0:
        d = _solve(lambda d: np.sinh(d) / d - B, 1e-8, 50.0)
        u = 0.5 * (1.0 + np.tanh(d * (xi - 0.5)) / np.tanh(0.5 * d))
    else:
        d = _solve(lambda d: np.sin(d) / d - B, 1e-8, np.pi - 1e-9)
        u = 0.5 * (1.0 + np.tan(d * (xi - 0.5)) / np.tan(0.5 * d))
    s = u / (A + (1.0 - A) * u)
    s[0], s[-1] = 0.0, 1.0
    return s * L


def _solve(f, lo, hi, tol=1e-14, maxiter=300):
    flo = f(lo)
    for _ in range(maxiter):
        mid = 0.5 * (lo + hi)
        fm = f(mid)
        if abs(hi - lo) < tol:
            break
        if (fm > 0) == (flo > 0):
            lo, flo = mid, fm
        else:
            hi = mid
    return 0.5 * (lo + hi)


def max_neighbour_ratio(ds):
    r = ds[1:] / ds[:-1]
    return float(np.max(np.maximum(r, 1.0 / r))) if len(r) else 1.0


def choose_bridge_count(L, ds0, ds1, rmax, dsmax, nmin=2, nmax=2000):
    """Smallest interval count whose Vinokur distribution keeps every
    neighbour ratio (including the two joins to the surfaces, whose
    spacings are ds0/ds1) at or below rmax and every spacing at or below
    dsmax (never tighter than the end spacings themselves)."""
    dsmax = max(dsmax, ds0, ds1) * (1.0 + 1e-9)
    for n in range(nmin, nmax):
        ds = np.diff(vinokur(L, n, ds0, ds1))
        full = np.concatenate([[ds0], ds, [ds1]])
        if max_neighbour_ratio(full) <= rmax and ds.max() <= dsmax:
            return n
    raise RuntimeError("no bridge count satisfies the spacing limits")


def sample_by_arclength(dense, stations):
    """Interpolate a densely sampled polyline at the given arc lengths."""
    arc = np.concatenate([[0.0], np.cumsum(seg_lengths(dense))])
    x = np.interp(stations, arc, dense[:, 0])
    y = np.interp(stations, arc, dense[:, 1])
    return np.column_stack([x, y])


def hermite_dense(P0, P1, t0, t1, pull0, pull1, n=4000):
    t = np.linspace(0.0, 1.0, n)
    return _hermite_points(t, np.asarray(P0, float), np.asarray(P1, float),
                           unit(t0) * pull0, unit(t1) * pull1)


def arc_length(dense):
    return float(np.sum(seg_lengths(dense)))


def segments_intersect(p1, p2, q1, q2):
    """Proper intersection of segments p1p2 and q1q2 (shared endpoints and
    collinear touching do not count)."""
    def orient(a, b, c):
        return (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])
    d1, d2 = orient(q1, q2, p1), orient(q1, q2, p2)
    d3, d4 = orient(p1, p2, q1), orient(p1, p2, q2)
    return (d1 * d2 < 0) and (d3 * d4 < 0)


def polyline_crossings(A, B=None):
    """Index pairs of crossing segments between polylines A and B, or of
    A with itself (non-adjacent segments only) when B is None. Uses a
    bounding-box prefilter, so a few thousand points is fine."""
    self_check = B is None
    if self_check:
        B = A
    a0, a1 = A[:-1], A[1:]
    b0, b1 = B[:-1], B[1:]
    amin, amax = np.minimum(a0, a1), np.maximum(a0, a1)
    bmin, bmax = np.minimum(b0, b1), np.maximum(b0, b1)
    hits = []
    for i in range(len(a0)):
        cand = np.where((bmin[:, 0] <= amax[i, 0]) & (bmax[:, 0] >= amin[i, 0]) &
                        (bmin[:, 1] <= amax[i, 1]) & (bmax[:, 1] >= amin[i, 1]))[0]
        for j in cand:
            if self_check and abs(i - j) <= 1:
                continue
            if self_check and j < i:
                continue
            if segments_intersect(a0[i], a1[i], b0[j], b1[j]):
                hits.append((i, int(j)))
    return hits


def point_to_polyline_distance(P, poly):
    """Distance from each point in P to the polyline poly."""
    a, b = poly[:-1], poly[1:]
    ab = b - a
    L2 = np.maximum(np.sum(ab * ab, axis=1), 1e-300)
    out = np.empty(len(P))
    for k, p in enumerate(P):
        t = np.clip(np.sum((p - a) * ab, axis=1) / L2, 0.0, 1.0)
        proj = a + t[:, None] * ab
        out[k] = np.min(np.hypot(*(p - proj).T))
    return out
