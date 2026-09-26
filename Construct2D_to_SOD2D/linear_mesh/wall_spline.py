"""Fit an open (non-periodic) natural cubic spline through the airfoil wall's
linear Construct2D corner points, and export it as a compact corner +
per-segment-coefficient table consumed by SOD2D's MeshElasticitySolver.

The trailing edge is treated as a genuine geometric corner (a finite wedge
angle, or at best a cusp), never as a smooth wraparound point -- so the
spline is always fit as an open curve with independent (natural) end
conditions at each end, even for the O-grid where both ends of the input
point sequence are numerically the same physical node.

A blunt trailing edge closed by a straight line (Construct2D's straight-
line-closure option -- several grid points sitting on one straight segment
between the upper- and lower-surface TE corners) is kept exactly straight:
that closure run is detected at each end of the wall point sequence (its
points are collinear and it terminates in a sharp corner within a short
arc length) and its segments are written as plain linear coefficients, so
the elasticity solver never bows them out. The cubic spline is then fit
through only the genuine curved surface, between those two TE corners. A
sharp/cusp TE trips nothing here -- no straight run is found, and the fit
is the same single open spline through every wall point as before.

Fortran places each high-order boundary node's target position directly and
parametrically (arc-length interpolation between its two bounding corners,
evaluated via this spline's closed-form per-segment cubic), not via a
nearest-point search against a dense sampling -- see MeshElasticitySolver.f90
for why: a search against a dense curve sampling is a many-to-one map near
high curvature (two distinct nodes straddling the leading edge can share the
same nearest sample), while arc-length is a single monotonic parameter, so
this table only needs to carry the corners themselves and the spline's own
segment coefficients, not a dense sampling.

No scipy dependency: this implements the standard natural-cubic-spline
tridiagonal system directly via a Thomas-algorithm solve.
"""

import numpy as np


def chord_length_param(x, y):
    """Cumulative chord-length parametrization; s[0] = 0."""
    d = np.hypot(np.diff(x), np.diff(y))
    return np.concatenate(([0.0], np.cumsum(d)))


def _turn_metric(x, y, i):
    """Normalized |cross product| of the two edges meeting at corner i:
    0.0 == the three points are collinear, 1.0 == a right-angle bend.
    Scale-free, so a threshold on it is a plain angle test."""
    e1 = np.array([x[i] - x[i - 1], y[i] - y[i - 1]])
    e2 = np.array([x[i + 1] - x[i], y[i + 1] - y[i]])
    n1 = np.hypot(*e1)
    n2 = np.hypot(*e2)
    if n1 == 0.0 or n2 == 0.0:
        return 0.0
    return abs(e1[0] * e2[1] - e1[1] * e2[0]) / (n1 * n2)


def detect_straight_te_closure(x_wall, y_wall, turn_threshold=0.2,
                               max_arclen_frac=0.1, chord_parallel_max=0.3):
    """Find the two trailing-edge corners bounding a straight-line TE closure.

    Walks in from each end of the ordered wall point sequence looking for a
    short run that (a) stays collinear (``_turn_metric`` below the
    threshold), (b) terminates in a sharp corner within
    ``max_arclen_frac`` of the total wall arc length, and (c) runs roughly
    ACROSS the trailing edge -- its direction is nearly perpendicular to
    the chord (|cos angle to chord| <= ``chord_parallel_max``). Condition
    (c) is what tells a blunt straight base apart from the first curved-
    surface segment near a sharp wedge TE: the base spans the TE
    thickness, the surface heads back chordwise. A sharp/cusp TE passes
    neither the collinearity run nor (c), so nothing is trimmed.

    Returns ``(k_start, k_end)``: the indices of the first and last corner
    of the genuine curved surface. ``k_start == 0`` / ``k_end == n-1``
    means no straight closure was found on that side.
    """
    x_wall = np.asarray(x_wall, dtype=float)
    y_wall = np.asarray(y_wall, dtype=float)
    n = len(x_wall)
    s = chord_length_param(x_wall, y_wall)
    total = s[-1]

    # Chord unit vector: trailing edge (corner 0) to the wall corner
    # farthest from it (the leading edge) -- orientation-agnostic.
    te = np.array([x_wall[0], y_wall[0]])
    d_from_te = np.hypot(x_wall - te[0], y_wall - te[1])
    le = np.array([x_wall[d_from_te.argmax()], y_wall[d_from_te.argmax()]])
    chord = te - le
    chord /= np.hypot(*chord) or 1.0

    def is_across_te(i0, i1):
        v = np.array([x_wall[i1] - x_wall[i0], y_wall[i1] - y_wall[i0]])
        nv = np.hypot(*v)
        if nv == 0.0:
            return False
        return abs(np.dot(v / nv, chord)) <= chord_parallel_max

    def corner_from_start():
        for k in range(1, n - 1):
            if s[k] - s[0] > max_arclen_frac * total:
                return 0
            if _turn_metric(x_wall, y_wall, k) >= turn_threshold:
                return k if is_across_te(0, k) else 0
        return 0

    def corner_from_end():
        for k in range(n - 2, 0, -1):
            if s[-1] - s[k] > max_arclen_frac * total:
                return n - 1
            if _turn_metric(x_wall, y_wall, k) >= turn_threshold:
                return k if is_across_te(k, n - 1) else n - 1
        return n - 1

    k_start = corner_from_start()
    k_end = corner_from_end()
    # Need a real surface span left between the two corners.
    if k_end - k_start < 3:
        return 0, n - 1
    return k_start, k_end


def natural_cubic_spline_coeffs(s, y):
    """Natural cubic spline (zero second derivative at both ends).

    Returns per-segment coefficients (a, b, c, d) such that, on segment i
    (s[i] <= t <= s[i+1]):
        y(t) = a[i] + b[i]*dt + c[i]*dt**2 + d[i]*dt**3,  dt = t - s[i]
    """
    n = len(s) - 1  # number of segments
    h = np.diff(s)

    # Tridiagonal system for the second derivatives (natural spline: M[0]=M[n]=0)
    A = np.zeros(n + 1)
    B = np.zeros(n + 1)
    C = np.zeros(n + 1)
    D = np.zeros(n + 1)

    B[0] = 1.0
    B[n] = 1.0
    for i in range(1, n):
        A[i] = h[i - 1]
        B[i] = 2.0 * (h[i - 1] + h[i])
        C[i] = h[i]
        D[i] = 6.0 * ((y[i + 1] - y[i]) / h[i] - (y[i] - y[i - 1]) / h[i - 1])

    # Thomas algorithm (forward sweep + back substitution)
    cp = np.zeros(n + 1)
    dp = np.zeros(n + 1)
    cp[0] = C[0] / B[0]
    dp[0] = D[0] / B[0]
    for i in range(1, n + 1):
        denom = B[i] - A[i] * cp[i - 1]
        cp[i] = C[i] / denom if i < n else 0.0
        dp[i] = (D[i] - A[i] * dp[i - 1]) / denom

    M = np.zeros(n + 1)
    M[n] = dp[n]
    for i in range(n - 1, -1, -1):
        M[i] = dp[i] - cp[i] * M[i + 1]

    a = y[:-1]
    b = (y[1:] - y[:-1]) / h - h * (2.0 * M[:-1] + M[1:]) / 6.0
    c = M[:-1] / 2.0
    d = (M[1:] - M[:-1]) / (6.0 * h)

    return a, b, c, d


def evaluate_spline(s, coeffs, s_query):
    """Evaluate a spline built by natural_cubic_spline_coeffs at s_query
    (assumed sorted, within [s[0], s[-1]])."""
    a, b, c, d = coeffs
    idx = np.searchsorted(s, s_query, side='right') - 1
    idx = np.clip(idx, 0, len(a) - 1)
    dt = s_query - s[idx]
    return a[idx] + dt * (b[idx] + dt * (c[idx] + dt * d[idx]))


def fit_wall_spline(x_wall, y_wall):
    """Fit an open cubic spline through the ordered wall corner points.

    Returns (s, cx, cy) where s is the per-corner cumulative arc-length and
    cx/cy are the (a,b,c,d) per-segment coefficient tuples for x(s)/y(s).
    """
    x_wall = np.asarray(x_wall, dtype=float)
    y_wall = np.asarray(y_wall, dtype=float)

    s = chord_length_param(x_wall, y_wall)
    cx = natural_cubic_spline_coeffs(s, x_wall)
    cy = natural_cubic_spline_coeffs(s, y_wall)
    return s, cx, cy


def check_corner_separation(x_wall, y_wall, min_ratio=0.1):
    """Sanity check: no two *non-adjacent* corners should sit closer together
    than min_ratio times the smallest adjacent-segment length. Catches a
    pathological near-duplicate corner before it ever reaches Fortran's
    corner-matching step (which relies on corners being mutually
    well-separated). Raises ValueError if violated.
    """
    x_wall = np.asarray(x_wall, dtype=float)
    y_wall = np.asarray(y_wall, dtype=float)
    n = len(x_wall)
    seg_len = np.hypot(np.diff(x_wall), np.diff(y_wall))
    min_seg = seg_len.min()
    threshold = min_ratio * min_seg

    for i in range(n):
        for j in range(i + 2, n):
            if i == 0 and j == n - 1:
                continue  # the two open ends may coincide (TE); that's expected
            dist = np.hypot(x_wall[i] - x_wall[j], y_wall[i] - y_wall[j])
            if dist < threshold:
                raise ValueError(
                    f"Wall corners {i} and {j} are only {dist:.3e} apart "
                    f"(< {min_ratio} x min segment length {min_seg:.3e}) -- "
                    f"too close for unambiguous corner matching in Fortran."
                )


def _linear_segment_coeffs(p, s, i):
    """(a, b, c, d) for a straight segment corner i -> i+1 of quantity p,
    parametrized by arc length s: p(t) = p[i] + b*(t - s[i]), c = d = 0."""
    h = s[i + 1] - s[i]
    return p[i], (p[i + 1] - p[i]) / h, 0.0, 0.0


def write_wall_spline_table(path, x_wall, y_wall, preserve_straight_te=True,
                            te_turn_threshold=0.2, te_max_arclen_frac=0.1):
    """Fit the spline and write the compact corner + coefficient table:

      line 1: n_corners
      then n_corners rows: x  y  s
      then (n_corners-1) rows: ax bx cx dx  ay by cy dy
        (segment i's coefficients, matching corner i -> corner i+1)

    With ``preserve_straight_te`` (default), a straight-line trailing-edge
    closure at either end of the wall sequence is detected
    (``detect_straight_te_closure``) and its segments are written as plain
    linear coefficients; the cubic spline is fit through only the curved
    surface between the two TE corners. A sharp/cusp TE finds no such run
    and the result is the single open spline through all corners as before.
    """
    check_corner_separation(x_wall, y_wall)
    x_wall = np.asarray(x_wall, dtype=float)
    y_wall = np.asarray(y_wall, dtype=float)
    n = len(x_wall)
    s = chord_length_param(x_wall, y_wall)

    k_start, k_end = 0, n - 1
    if preserve_straight_te:
        k_start, k_end = detect_straight_te_closure(
            x_wall, y_wall, te_turn_threshold, te_max_arclen_frac
        )

    # Cubic spline through the curved surface span only (k_start..k_end).
    # Its own arc-length sub-range keeps the same values as the global s,
    # so its per-segment coefficients drop straight into global segments
    # k_start .. k_end-1 (Fortran evaluates dt = t - s[global_i]).
    ax_s, bx_s, cx_s, dx_s = natural_cubic_spline_coeffs(s[k_start:k_end + 1], x_wall[k_start:k_end + 1])
    ay_s, by_s, cy_s, dy_s = natural_cubic_spline_coeffs(s[k_start:k_end + 1], y_wall[k_start:k_end + 1])

    ax = np.zeros(n - 1); bx = np.zeros(n - 1); ccx = np.zeros(n - 1); dx = np.zeros(n - 1)
    ay = np.zeros(n - 1); by = np.zeros(n - 1); ccy = np.zeros(n - 1); dy = np.zeros(n - 1)
    for i in range(n - 1):
        if k_start <= i < k_end:
            j = i - k_start
            ax[i], bx[i], ccx[i], dx[i] = ax_s[j], bx_s[j], cx_s[j], dx_s[j]
            ay[i], by[i], ccy[i], dy[i] = ay_s[j], by_s[j], cy_s[j], dy_s[j]
        else:  # straight TE-closure segment
            ax[i], bx[i], ccx[i], dx[i] = _linear_segment_coeffs(x_wall, s, i)
            ay[i], by[i], ccy[i], dy[i] = _linear_segment_coeffs(y_wall, s, i)

    if (k_start, k_end) != (0, n - 1):
        print(f"  wall spline: straight TE closure preserved -- linear segments "
              f"[0..{k_start - 1}] and [{k_end}..{n - 2}], cubic spline through corners {k_start}..{k_end}")

    with open(path, 'w') as f:
        f.write(f"{n}\n")
        for i in range(n):
            f.write(f"{x_wall[i]:.10e} {y_wall[i]:.10e} {s[i]:.10e}\n")
        for i in range(n - 1):
            f.write(
                f"{ax[i]:.10e} {bx[i]:.10e} {ccx[i]:.10e} {dx[i]:.10e} "
                f"{ay[i]:.10e} {by[i]:.10e} {ccy[i]:.10e} {dy[i]:.10e}\n"
            )


def _extract_2d_wall_run(x2d, y2d, face, s1, e1):
    """Extract a boundary's own point run from a raw 2D (idim, jmax) grid,
    generically by its 2D NMF face id -- face3 (j=0, the usual wall row)
    varies I with J fixed at 0, but a TE-closure block's own wall can just
    as well be an i-const face (face1: i=0, or face2: i=idim-1), varying J
    instead. s1/e1 is the corresponding 1-based inclusive range (I range
    for face3/4, J range for face1/2 -- see mesh_extrusion._read_2d_nmf /
    _remap_2d_face_to_3d for the same per-face-parity convention used
    throughout this pipeline)."""
    if face == 3:
        return x2d[s1 - 1: e1, 0], y2d[s1 - 1: e1, 0]
    if face == 4:
        return x2d[s1 - 1: e1, -1], y2d[s1 - 1: e1, -1]
    if face == 1:
        return x2d[0, s1 - 1: e1], y2d[0, s1 - 1: e1]
    if face == 2:
        return x2d[-1, s1 - 1: e1], y2d[-1, s1 - 1: e1]
    raise ValueError(f"Unexpected 2D face id {face} for a VISCOUS boundary")


def _assemble_2block_wall(blocks, viscous):
    """Assemble the full wall polyline for a 2-block C-grid where the
    blunt trailing edge is resolved as its own small block (its wall is
    the straight TE-closure line) instead of a couple of extra rows
    inside one O-grid -- so the two VISCOUS boundaries live in two
    different blocks' arrays instead of one.

    The longer VISCOUS run is the main curved airfoil surface, open-
    ended at the two TE corners; the shorter one is the TE-closure
    block's own short run between the SAME two corners. Concatenating
    them the right way round closes the loop exactly like the O-grid's
    i=0/i=idim-1 wrap (build_wall_spline_table's OGRD branch): the two
    open ends of the combined sequence coincide in coordinates but are
    not adjacent in the fitting order, so the closing segment is a real,
    normal-length wall element.

    Which end of the shorter run needs to come first is decided by
    matching its own endpoint coordinates against the main run's open
    ends (nearest match wins) -- not by assuming a fixed NMF index
    convention, since the ONE_TO_ONE welds already handle connectivity
    for the volume mesh independently of this.
    """
    main, cap = sorted(viscous, key=lambda b: b['e1'] - b['s1'], reverse=True)  # longer run first
    xm, ym = blocks[main['b1'] - 1]
    xc, yc = blocks[cap['b1'] - 1]

    main_x, main_y = _extract_2d_wall_run(xm, ym, main['f1'], main['s1'], main['e1'])
    cap_x, cap_y = _extract_2d_wall_run(xc, yc, cap['f1'], cap['s1'], cap['e1'])

    end = np.array([main_x[-1], main_y[-1]])
    d_start = np.hypot(cap_x[0] - end[0], cap_y[0] - end[1])
    d_end = np.hypot(cap_x[-1] - end[0], cap_y[-1] - end[1])
    if d_end < d_start:
        cap_x, cap_y = cap_x[::-1], cap_y[::-1]
        d_start = d_end

    main_span = np.hypot(main_x[-1] - main_x[0], main_y[-1] - main_y[0]) or 1.0
    if d_start > 1e-6 * main_span:
        raise ValueError(
            f"2-block wall assembly: the TE-closure block's VISCOUS run "
            f"doesn't share an endpoint with the main block's wall "
            f"(closest distance {d_start:.3e})"
        )

    # Drop the cap run's now-leading point -- it coincides with main's
    # last point (both the same physical TE corner); its trailing point
    # coincides with main's first point but, like the O-grid's own
    # closure, is kept as a distinct non-adjacent entry, not dropped.
    wall_x = np.concatenate([main_x, cap_x[1:]])
    wall_y = np.concatenate([main_y, cap_y[1:]])
    return wall_x, wall_y


def build_wall_spline_table(blocks, boundaries, mesh_type, out_path):
    """Extract the wall sub-range(s) from the raw 2D grid, fit the spline,
    and write the compact corner+coefficient table to out_path.

    :blocks: list of (x2d, y2d) raw 2D Plot3D grids, one per block (index
        0 = block 1, 1-based in the NMF), each shape (idim, jmax) with
        j=0 the block's own wall row.
    :boundaries: as returned by mesh_extrusion._read_2d_nmf /
        write_ext_nmf -- covers every block.
    :mesh_type: "OGRD" or "CGRD".
    """
    x2d, y2d = blocks[0]
    if mesh_type.upper() == "OGRD":
        # i=0 and i=idim-1 are the same physical (aliased) node. Fit through
        # ALL idim raw points (not idim-1): for this open, non-periodic
        # spline, a zero-length segment only occurs if the duplicate point
        # is explicitly re-appended and fit adjacent to itself -- point
        # idim-1 here merely coincides with point 0 in coordinates, it is
        # not adjacent to it in the fitting sequence, so the last segment
        # (idim-2 -> idim-1) has completely normal nonzero length and is a
        # real wall element (the one closing the loop back at the TE).
        # Dropping it (as an earlier version of this code did) silently
        # discards that element from the fit entirely.
        wall_x = x2d[:, 0]
        wall_y = y2d[:, 0]
    elif len(blocks) == 1:
        # NMF s1/e1 are 1-based and INCLUSIVE on both ends (standard NMF/Plot3D
        # convention -- confirmed against mesh_extrusion.py's own _remap_2d_face_to_3d,
        # which passes s1/e1 straight through as inclusive block ranges). A Python
        # slice's stop is exclusive, so the endpoint is `e1` (not `e1 - 1`) to include
        # the true last node. Getting this wrong silently drops the VISCOUS range's
        # last point -- for a C-grid, that last point coincides exactly with the first
        # (both are the trailing edge; the wall wraps around the whole airfoil just
        # like the O-grid's i=0/i=idim-1 closure), so dropping it isn't just "one point
        # short", it discards the real wall element on the other side of the TE from
        # the fit entirely, the same class of bug already fixed for the O-grid above.
        visc = next(b for b in boundaries if b['name'] == 'VISCOUS')
        wall_x = x2d[visc['s1'] - 1: visc['e1'], 0]
        wall_y = y2d[visc['s1'] - 1: visc['e1'], 0]
    else:
        viscous = [b for b in boundaries if b['name'] == 'VISCOUS']
        if len(blocks) != 2 or len(viscous) != 2:
            raise ValueError(
                f"build_wall_spline_table: unsupported multiblock configuration "
                f"({len(viscous)} VISCOUS boundaries across {len(blocks)} blocks) "
                f"-- only the 2-block main-C-grid + TE-closure-block CGRD shape "
                f"is supported (2 VISCOUS boundaries across 2 blocks)."
            )
        wall_x, wall_y = _assemble_2block_wall(blocks, viscous)

    write_wall_spline_table(out_path, wall_x, wall_y)
    return out_path
