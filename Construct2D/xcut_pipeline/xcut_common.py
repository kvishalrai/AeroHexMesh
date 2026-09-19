"""
Shared I/O and geometry helpers for the OAT15-style blunt-trailing-edge
XCUT meshing pipeline (zone (a) construction).

Curve/mesh conventions used throughout this pipeline (established and
validated interactively before being turned into these scripts):

  - Airfoil chord is normalized with the trailing edge at x = 1.0.
  - The assembled "whole curve" fed to Construct2D (CGRD, NWKE=0, BUFF)
    is ordered:
        [close point (Lfar, 0.0)]
      + [bottom wake, Lfar -> 1.0, y = -h, EXCLUDING the shared bottomTE
         point which starts the airfoil segment]
      + [airfoil surface, bottomTE -> (around the LE) -> topTE]
      + [top wake, 1.0 -> Lfar, y = +h, EXCLUDING the shared topTE point
         which ends the airfoil segment]
      + [close point (Lfar, 0.0), identical to the first point]
    The leading/trailing (Lfar, 0.0) duplicate is a deliberate trick: it
    makes x(1)==x(npoints) and y(1)==y(npoints) exactly, which is what
    Construct2D's setup_airfoil_data checks to decide the curve is
    "closed" (tegap=False), which is required to use the BUFF option.
    It also means the elements immediately adjacent to i=1 and i=imax
    are geometrically degenerate (near-zero-area slivers connecting the
    true wake tips at y=+-h to the fake apex at y=0) -- this is expected,
    and is exactly what step 3 (trimming) removes.
  - With NWKE=0, Construct2D treats the WHOLE input curve as "the
    surface" 1:1 -- i.e. the resulting p3d's i=1..imax exactly matches
    the input curve's points 1..npoints, with no insertion/removal.
    This lets us compute bottomTE/topTE column indices in the output
    mesh purely from the point counts used to build the curve (see
    curve_indices() below), with no need to re-detect them from the
    mesh afterwards.

All files this pipeline writes are NEW files (never overwrites an
existing file in place), matching the user's standing preference.
"""
import json
import numpy as np


# ---------------------------------------------------------------------
# Airfoil .dat curve I/O
# ---------------------------------------------------------------------

def read_airfoil_dat(fn):
    """Read a simple (header line + x y pairs) airfoil .dat file."""
    with open(fn) as f:
        lines = [l for l in f if l.strip()]
    header = lines[0].strip()
    pts = np.array([[float(v) for v in l.split()[:2]] for l in lines[1:]])
    return header, pts[:, 0], pts[:, 1]


def write_airfoil_dat(fn, header, x, y):
    # 17 significant digits round-trips any float64 exactly (re-parsing the
    # written text reproduces the identical bit pattern) -- required since
    # unmodified airfoil-surface coordinates must come back byte-identical,
    # not merely "close", and this project's source .dat files carry more
    # precision than a 10-decimal fixed format preserves.
    with open(fn, 'w') as f:
        f.write(f"{header}\n")
        for xi, yi in zip(x, y):
            f.write(f"{xi:24.17g}  {yi:24.17g}\n")


def strip_te_closure_artifact(x, y, tol=1e-9):
    """
    Many airfoil files in this project (oat15.dat, oat15_slc.dat) are
    pre-closed for some earlier meshing attempt by literally repeating
    a single point at index 0 and index -1 (typically the TE-face
    midpoint, e.g. (1.0, 0.0)) so that Construct2D sees a closed curve.

    We don't want that convention here -- this pipeline builds its own
    closure out at x=Lfar instead. So: if the first and last points of
    the input curve coincide (to `tol`), drop both, keeping only the
    genuine airfoil-surface points in between (which may still include
    real TE-face geometry, e.g. a half-gap point, if the source file
    resolves the blunt TE face with more than one point per side --
    that's kept, since it's real surface, not an artifact).

    Returns the (possibly trimmed) x, y and a bool saying whether
    anything was stripped.
    """
    if len(x) >= 2 and abs(x[0] - x[-1]) <= tol and abs(y[0] - y[-1]) <= tol:
        return x[1:-1], y[1:-1], True
    return x, y, False


# ---------------------------------------------------------------------
# PLOT3D I/O (single block and multi-block, matching this project's
# established fixed-width ASCII format: one value per line, "%17.8E")
# ---------------------------------------------------------------------

def read_p3d(fn):
    with open(fn) as f:
        first = f.readline().split()
        if len(first) == 1:
            # multi-block: first line = nblocks, then one "imax jmax" per block
            nblk = int(first[0])
            dims = [tuple(map(int, f.readline().split())) for _ in range(nblk)]
            vals = []
            for line in f:
                vals.extend(line.split())
            vals = np.array(vals, dtype=float)
            blocks = []
            pos = 0
            for imax, jmax in dims:
                n = imax * jmax
                x = vals[pos:pos + n].reshape(jmax, imax).T
                y = vals[pos + n:pos + 2 * n].reshape(jmax, imax).T
                pos += 2 * n
                blocks.append((x, y))
            return blocks
        else:
            imax, jmax = int(first[0]), int(first[1])
            vals = []
            for line in f:
                vals.extend(line.split())
            vals = np.array(vals, dtype=float)
            n = imax * jmax
            x = vals[:n].reshape(jmax, imax).T
            y = vals[n:2 * n].reshape(jmax, imax).T
            return [(x, y)]


def write_p3d(fn, x, y):
    """Single-block p3d."""
    imax, jmax = x.shape
    with open(fn, 'w') as f:
        f.write(f"{imax:12d}{jmax:12d}\n")
        for arr in (x, y):
            for v in arr.T.flatten():
                f.write(f"{v:17.8E}\n")


def write_p3d_multiblock(fn, blocks):
    """blocks: list of (x, y) arrays, each shape (imax_k, jmax_k)."""
    with open(fn, 'w') as f:
        f.write(f"{len(blocks):12d}\n")
        for x, y in blocks:
            imax, jmax = x.shape
            f.write(f"{imax:12d}{jmax:12d}\n")
        for x, y in blocks:
            for arr in (x, y):
                for v in arr.T.flatten():
                    f.write(f"{v:17.8E}\n")


# ---------------------------------------------------------------------
# PLOT3D "function" file I/O -- Construct2D's *_stats.p3d mesh-quality
# files (# comment lines, then "imax jmax kmax nvar", then nvar scalar
# fields in the same per-point layout as the coordinate p3d).
# ---------------------------------------------------------------------

def read_p3d_function(fn):
    with open(fn) as f:
        lines = f.readlines()
    comments = []
    i = 0
    while lines[i].lstrip().startswith('#'):
        comments.append(lines[i].rstrip('\n'))
        i += 1
    imax, jmax, kmax, nvar = map(int, lines[i].split())
    i += 1
    vals = []
    for line in lines[i:]:
        vals.extend(line.split())
    vals = np.array(vals, dtype=float)
    n = imax * jmax * kmax
    arrays = []
    pos = 0
    for _ in range(nvar):
        arrays.append(vals[pos:pos + n].reshape(jmax, imax).T)  # kmax assumed 1
        pos += n
    return comments, imax, jmax, kmax, arrays


def write_p3d_function(fn, comments, imax, jmax, kmax, arrays):
    with open(fn, 'w') as f:
        for c in comments:
            f.write(c + "\n")
        f.write(f"{imax:12d}{jmax:12d}{kmax:12d}{len(arrays):12d}\n")
        for a in arrays:
            for v in a.T.flatten():
                f.write(f"{v:17.8E}\n")


def write_p3d_function_multiblock(fn, comments, dims, blocks_arrays):
    """
    Multi-block counterpart of write_p3d_function, for quality-stats files
    covering more than one grid block (e.g. this pipeline's own 2-block
    XCUT mesh). Mirrors write_p3d_multiblock's own convention -- a block
    count line, then one dims line per block, then all blocks' data
    concatenated block-by-block -- but with an added nvar column on each
    dims line (matching Construct2D's own single-block "imax jmax kmax
    nvar" header) so a reader never has to be told nvar out of band.

    `dims`: list of (imax, jmax, kmax) per block.
    `blocks_arrays`: list (one entry per block) of lists of nvar 2D arrays
    each (same variables, same order, in every block).
    """
    nvar = len(blocks_arrays[0])
    with open(fn, 'w') as f:
        for c in comments:
            f.write(c + "\n")
        f.write(f"{len(dims):12d}\n")
        for imax, jmax, kmax in dims:
            f.write(f"{imax:12d}{jmax:12d}{kmax:12d}{nvar:12d}\n")
        for arrays in blocks_arrays:
            for a in arrays:
                for v in a.T.flatten():
                    f.write(f"{v:17.8E}\n")


def read_p3d_function_multiblock(fn):
    """Inverse of write_p3d_function_multiblock: returns
    (comments, dims, blocks_arrays)."""
    with open(fn) as f:
        lines = f.readlines()
    comments = []
    i = 0
    while lines[i].lstrip().startswith('#'):
        comments.append(lines[i].rstrip('\n'))
        i += 1
    nblk = int(lines[i].split()[0])
    i += 1
    dims = []
    nvar = None
    for _ in range(nblk):
        imax, jmax, kmax, nvar = map(int, lines[i].split())
        dims.append((imax, jmax, kmax))
        i += 1
    vals = []
    for line in lines[i:]:
        vals.extend(line.split())
    vals = np.array(vals, dtype=float)
    blocks_arrays = []
    pos = 0
    for imax, jmax, kmax in dims:
        n = imax * jmax * kmax
        arrays = []
        for _ in range(nvar):
            arrays.append(vals[pos:pos + n].reshape(jmax, imax).T)  # kmax assumed 1
            pos += n
        blocks_arrays.append(arrays)
    return comments, dims, blocks_arrays


# ---------------------------------------------------------------------
# Grid quality metrics -- Python port of Construct2D's own
# compute_quality_stats (src/surface_grid.f90) and its angle()/growth()
# helpers (src/math_deps.f90). Used to generate the same skew-angle /
# xi-growth / eta-growth metrics Construct2D itself writes to a
# *_stats.p3d file, for blocks this pipeline builds directly (e.g. the
# TE gap-filler box) that never go through Construct2D and so never get
# a stats file from it.
# ---------------------------------------------------------------------

def compute_quality_stats(x, y):
    """
    x, y: (imax, jmax) coordinate arrays for one structured block.
    Returns (skewang, growthz, growthn), each (imax, jmax):
      - skewang: max deviation from 90 deg of the 4 angles formed by the
        grid lines crossing at each point (0 = perfectly orthogonal).
      - growthz: normalized cell-size change in the xi (i) direction
        (0 at j=1/jmax edges and at every point immediately along an
        i=1/imax edge, by the same convention Construct2D uses).
      - growthn: normalized cell-size change in the eta (j) direction.

    Verified against a real Construct2D mesh + its own *_stats.p3d
    (naca0012.p3d/_stats.p3d): skewang and growthn match Construct2D's
    own reported values directly (mean abs diff ~1e-4, i.e. within the
    stats file's own 8-sig-fig write precision); growthz matches only
    after negating math_deps.f90's growth() result here (mean abs diff
    ~2e-7 once negated, vs ~0.06 unnegated) -- Construct2D's own
    convention for the sign of xi-direction growth is the opposite of
    what a literal reading of growth(p_{i+1}, p_i, p_{i-1}) gives, for
    reasons not evident from the source alone; this reproduces its
    actual reported values rather than the naive formula.
    """
    imax, jmax = x.shape

    def pt(i, j):
        return np.array([x[i, j], y[i, j]])

    def ang(p1, p2, p0):
        v1, v2 = p1 - p0, p2 - p0
        m1, m2 = np.hypot(v1[0], v1[1]), np.hypot(v2[0], v2[1])
        c = np.dot(v1, v2) / (m1 * m2)
        c = min(1.0, max(-1.0, c))
        return np.degrees(np.arccos(c))

    def grow(p1, p0, pm1):
        len1 = np.hypot(*(p1 - p0))
        len0 = np.hypot(*(p0 - pm1))
        return (len1 - len0) / min(len0, len1)

    skewang = np.zeros((imax, jmax))
    growthz = np.zeros((imax, jmax))
    growthn = np.zeros((imax, jmax))

    for j in range(jmax):
        for i in range(imax):
            left, right = i == 0, i == imax - 1
            bot, top = j == 0, j == jmax - 1
            p0 = pt(i, j)

            if left and bot:
                a1, a2, a3, a4 = ang(pt(i+1,j), pt(i,j+1), p0), 90.0, 90.0, 90.0
                gz = gn = 0.0
            elif left and top:
                a1, a2, a3, a4 = 90.0, 90.0, 90.0, ang(pt(i,j-1), pt(i+1,j), p0)
                gz = gn = 0.0
            elif right and top:
                a1, a2, a3, a4 = 90.0, 90.0, ang(pt(i-1,j), pt(i,j-1), p0), 90.0
                gz = gn = 0.0
            elif right and bot:
                a1, a2, a3, a4 = 90.0, ang(pt(i,j+1), pt(i-1,j), p0), 90.0, 90.0
                gz = gn = 0.0
            elif left:
                a1 = ang(pt(i+1,j), pt(i,j+1), p0)
                a2, a3 = 90.0, 90.0
                a4 = ang(pt(i,j-1), pt(i+1,j), p0)
                gz, gn = 0.0, grow(pt(i,j+1), p0, pt(i,j-1))
            elif right:
                a1 = 90.0
                a2 = ang(pt(i,j+1), pt(i-1,j), p0)
                a3 = ang(pt(i-1,j), pt(i,j-1), p0)
                a4 = 90.0
                gz, gn = 0.0, grow(pt(i,j+1), p0, pt(i,j-1))
            elif bot:
                a1 = ang(pt(i+1,j), pt(i,j+1), p0)
                a2 = ang(pt(i,j+1), pt(i-1,j), p0)
                a3, a4 = 90.0, 90.0
                gz, gn = grow(pt(i+1,j), p0, pt(i-1,j)), 0.0
            elif top:
                a1, a2 = 90.0, 90.0
                a3 = ang(pt(i-1,j), pt(i,j-1), p0)
                a4 = ang(pt(i,j-1), pt(i+1,j), p0)
                gz, gn = grow(pt(i+1,j), p0, pt(i-1,j)), 0.0
            else:
                a1 = ang(pt(i+1,j), pt(i,j+1), p0)
                a2 = ang(pt(i,j+1), pt(i-1,j), p0)
                a3 = ang(pt(i-1,j), pt(i,j-1), p0)
                a4 = ang(pt(i,j-1), pt(i+1,j), p0)
                gz = grow(pt(i+1,j), p0, pt(i-1,j))
                gn = grow(pt(i,j+1), p0, pt(i,j-1))

            skewang[i, j] = max(90.0 - abs(a1), 90.0 - abs(a2),
                                 90.0 - abs(a3), 90.0 - abs(a4))
            growthz[i, j] = -gz  # sign convention -- see docstring
            growthn[i, j] = gn

    return skewang, growthz, growthn


# ---------------------------------------------------------------------
# Metadata sidecar (keeps step-to-step bookkeeping explicit instead of
# re-detecting indices from geometry at every stage)
# ---------------------------------------------------------------------

def write_meta(fn, meta):
    with open(fn, 'w') as f:
        json.dump(meta, f, indent=2)


def read_meta(fn):
    with open(fn) as f:
        return json.load(f)


# ---------------------------------------------------------------------
# Geometric-growth wake line construction
# ---------------------------------------------------------------------

def geometric_growth_rate(L, d0, n, tol=1e-13, maxiter=200):
    """Solve d0*(r**(n-1)-1)/(r-1) = L for r > 1 (n-1 intervals)."""
    if n < 2:
        raise ValueError("need at least 2 points")
    if abs(n - 1) * d0 >= L:
        # uniform spacing already covers/exceeds L
        return 1.0
    lo, hi = 1.0 + 1e-12, 10.0
    def f(r):
        if abs(r - 1.0) < 1e-14:
            return d0 * (n - 1) - L
        return d0 * (r**(n - 1) - 1) / (r - 1) - L
    flo, fhi = f(lo), f(hi)
    while fhi < 0:
        hi *= 2
        fhi = f(hi)
        if hi > 1e8:
            raise RuntimeError("could not bracket growth ratio")
    for _ in range(maxiter):
        mid = 0.5 * (lo + hi)
        fm = f(mid)
        if abs(fm) < tol:
            return mid
        if (fm > 0) == (flo > 0):
            lo, flo = mid, fm
        else:
            hi, fhi = mid, fm
    return 0.5 * (lo + hi)


def build_wake_line(x0, y_const, d0, n, xfar):
    """n points from x0 (inclusive) to xfar (inclusive), constant y,
    geometric growth starting at step d0."""
    L = xfar - x0
    r = geometric_growth_rate(L, d0, n)
    xs = [x0]
    step = d0
    for _ in range(n - 1):
        xs.append(xs[-1] + step)
        step *= r
    xs = np.array(xs)
    xs[-1] = xfar  # kill roundoff drift
    ys = np.full(n, y_const)
    return np.column_stack([xs, ys]), r, L


def _hermite_points(t, P0, P1, m0, m1):
    t = np.asarray(t)
    t2, t3 = t * t, t * t * t
    h00 = (2 * t3 - 3 * t2 + 1)[:, None]
    h10 = (t3 - 2 * t2 + t)[:, None]
    h01 = (-2 * t3 + 3 * t2)[:, None]
    h11 = (t3 - t2)[:, None]
    return h00 * P0 + h10 * m0 + h01 * P1 + h11 * m1


def build_wake_centerline_spline(x0, y0, tangent0, xfar, yfar, tangent1,
                                  d0, n, pull_length=1.0, n_sample=4000):
    """
    A single smooth (cubic Hermite) curve from (x0,y0) with direction
    tangent0, to (xfar,yfar) with direction tangent1, sampled at n
    points with geometric-growth spacing (first step d0) measured
    along the curve's own arc length (not raw x).

    `pull_length` sets how far each tangent's influence reaches, as an
    ABSOLUTE length in the same units as the curve (chords, given this
    project's normalization) -- NOT a fraction of the total P0->P1
    distance. This matters here because both endpoints sit at y=0 (the
    curve returns to the far-field y-level), so tangent magnitude
    alone controls how far it wanders before flattening; scaling it by
    the full (possibly very long) wake distance would way overshoot.
    ~1 chord is a gentle, visually-sensible blend for a typical
    few-degree TE-bisector angle; scale up/down to taste.

    Intended use here: build ONE centerline this way, then offset it
    by +-h in y to get the top/bottom wake lines. Because both sides
    then share literally the same (x, y_centerline) samples shifted by
    a constant, they come out exactly parallel and exactly equidistant
    (gap = 2h everywhere) by construction, and with identical point
    distributions -- rather than being two independently-fit splines
    that could drift apart or wobble unevenly.
    """
    P0 = np.array([x0, y0], dtype=float)
    P1 = np.array([xfar, yfar], dtype=float)
    t0u = np.asarray(tangent0, dtype=float)
    t0u = t0u / np.linalg.norm(t0u)
    t1u = np.asarray(tangent1, dtype=float)
    t1u = t1u / np.linalg.norm(t1u)
    m0 = t0u * pull_length
    m1 = t1u * pull_length

    ts = np.linspace(0.0, 1.0, n_sample)
    pts = _hermite_points(ts, P0, P1, m0, m1)
    seg = np.hypot(np.diff(pts[:, 0]), np.diff(pts[:, 1]))
    arc = np.concatenate([[0.0], np.cumsum(seg)])
    L_arc = arc[-1]

    r = geometric_growth_rate(L_arc, d0, n)
    s = [0.0]
    step = d0
    for _ in range(n - 1):
        s.append(s[-1] + step)
        step *= r
    s = np.array(s)
    s[-1] = L_arc

    x_out = np.interp(s, arc, pts[:, 0])
    y_out = np.interp(s, arc, pts[:, 1])
    x_out[-1], y_out[-1] = xfar, yfar
    return np.column_stack([x_out, y_out]), r, L_arc
