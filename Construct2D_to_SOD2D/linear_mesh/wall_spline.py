"""Fit an open (non-periodic) natural cubic spline through the airfoil wall's
linear Construct2D corner points, and export it as a compact corner +
per-segment-coefficient table consumed by SOD2D's MeshElasticitySolver.

The trailing edge is treated as a genuine geometric corner (a finite wedge
angle, or at best a cusp), never as a smooth wraparound point -- so the
spline is always fit as an open curve with independent (natural) end
conditions at each end, even for the O-grid where both ends of the input
point sequence are numerically the same physical node.

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


def write_wall_spline_table(path, x_wall, y_wall):
    """Fit the spline and write the compact corner + coefficient table:

      line 1: n_corners
      then n_corners rows: x  y  s
      then (n_corners-1) rows: ax bx cx dx  ay by cy dy
        (segment i's coefficients, matching corner i -> corner i+1)
    """
    check_corner_separation(x_wall, y_wall)
    s, cx, cy = fit_wall_spline(x_wall, y_wall)
    ax, bx, ccx, dx = cx
    ay, by, ccy, dy = cy
    n = len(x_wall)

    with open(path, 'w') as f:
        f.write(f"{n}\n")
        for i in range(n):
            f.write(f"{x_wall[i]:.10e} {y_wall[i]:.10e} {s[i]:.10e}\n")
        for i in range(n - 1):
            f.write(
                f"{ax[i]:.10e} {bx[i]:.10e} {ccx[i]:.10e} {dx[i]:.10e} "
                f"{ay[i]:.10e} {by[i]:.10e} {ccy[i]:.10e} {dy[i]:.10e}\n"
            )


def build_wall_spline_table(x2d, y2d, boundaries, mesh_type, out_path):
    """Extract the wall sub-range from the raw 2D grid, fit the spline, and
    write the compact corner+coefficient table to out_path.

    :x2d, y2d: raw 2D Plot3D grid, shape (idim, jmax), j=0 is the wall row.
    :boundaries: as returned by mesh_extrusion._read_2d_nmf.
    :mesh_type: "OGRD" or "CGRD".
    """
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
    else:
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

    write_wall_spline_table(out_path, wall_x, wall_y)
    return out_path
