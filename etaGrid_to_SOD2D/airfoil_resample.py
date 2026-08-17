"""Resample an airfoil curve to a coarser/finer point distribution with
independently-specified leading-edge and trailing-edge point spacing.

Splits the closed airfoil loop into its two natural open segments (upper
surface: TE->LE, lower surface: LE->TE, sharing the LE and TE points --
the standard way airfoil meshing tools handle end-clustering, since a
single arc-length variable around a closed loop has no notion of "two
independently-tunable ends" the way an open curve does), then applies
Vinokur's (1983) two-sided stretching function to each segment
independently: given a point count and target first/last spacing, it's
the standard, exact technique for CFD surface point distribution -- it
matches the requested end-derivatives (spacing) at BOTH ends
simultaneously, unlike a single geometric-growth sequence (which can
only anchor one end).
"""
import numpy as np

from sweep_mesh import load_airfoil_curve, arc_length_param


def vinokur_stretch(n, ds0, ds1):
    """n points in [0,1] (normalized arc length), first spacing ds0 and
    last spacing ds1 (also normalized, i.e. as a fraction of total
    curve length). Standard two-sided (Vinokur) stretching -- see e.g.
    Vinokur, M. (1983), J. Comput. Phys. 50(2):215-234, or the standard
    formulation in Fletcher's "Computational Techniques for Fluid
    Dynamics" vol 2 section 12.4.

    ds0, ds1: physical first/last spacing DIVIDED by total curve length
    (so e.g. ds0=0.001 means the first interval is 0.1% of the curve).
    """
    if n < 2:
        raise ValueError("need at least 2 points")
    nseg = n - 1
    s0 = ds0 * nseg  # normalized by the UNIFORM spacing 1/nseg, per Vinokur's own convention
    s1 = ds1 * nseg

    A = np.sqrt(s1 / s0)
    B = 1.0 / np.sqrt(s0 * s1)

    if B > 1.0 + 1e-12:
        # Solve sinh(dy)/dy = B for dy via Newton's method.
        if B - 1.0 < 0.36:
            # Series expansion (Vinokur's own small-(B-1) approximation),
            # used as the Newton seed and directly when it's accurate
            # enough on its own.
            bm1 = B - 1.0
            dy = np.sqrt(6.0 * bm1) * (
                1.0 - 0.15 * bm1 + 0.057321429 * bm1**2 - 0.024907295 * bm1**3
                + 0.0077424461 * bm1**4 - 0.0010794123 * bm1**5
            )
        else:
            dy = np.log(2.0 * B) + np.log(np.log(2.0 * B)) if B > 2 else np.pi * (1.0 - 1.0 / B)
        for _ in range(50):
            f = np.sinh(dy) / dy - B
            fp = (np.cosh(dy) * dy - np.sinh(dy)) / dy**2
            step = f / fp
            dy -= step
            if abs(step) < 1e-13:
                break
        xi = np.linspace(0.0, 1.0, n)
        u = 0.5 * (1.0 + np.tanh(dy * (xi - 0.5)) / np.tanh(0.5 * dy))
    elif B < 1.0 - 1e-12:
        # B<1 means requested spacings are LARGER than uniform -- valid
        # but the sin-based branch applies (dy solves sin(dy)/dy = B).
        dy = np.pi * (1.0 - B)
        for _ in range(50):
            f = np.sin(dy) / dy - B
            fp = (np.cos(dy) * dy - np.sin(dy)) / dy**2
            step = f / fp
            dy_new = dy - step
            if dy_new <= 0 or dy_new >= np.pi:
                dy_new = 0.5 * (dy + (0 if step > 0 else np.pi))
            dy = dy_new
            if abs(step) < 1e-13:
                break
        xi = np.linspace(0.0, 1.0, n)
        u = 0.5 * (1.0 + np.tan(dy * (xi - 0.5)) / np.tan(0.5 * dy))
    else:
        xi = np.linspace(0.0, 1.0, n)
        u = xi.copy()

    t = u / (A + (1.0 - A) * u)
    t[0], t[-1] = 0.0, 1.0  # exact endpoints, avoid float drift
    return t


def resample_segment(pts, n, ds0, ds1):
    """pts: (M,2) open curve segment. Returns (n,2) resampled via
    Vinokur stretching, using a cubic-spline fit of the ORIGINAL segment
    for the actual position evaluation (so the resampled points still
    lie exactly on the true curve, just redistributed)."""
    from scipy.interpolate import CubicSpline

    s = arc_length_param(pts)
    total = s[-1]
    spline = CubicSpline(s, pts, axis=0, bc_type="not-a-knot")
    t = vinokur_stretch(n, ds0 / total, ds1 / total)
    return spline(t * total)


def resample_airfoil(pts, n_upper, n_lower, le_spacing, te_spacing):
    """pts: (M,2) closed airfoil loop (pts[0]==pts[-1], sharp TE, as
    loaded by sweep_mesh.load_airfoil_curve). Splits at the point
    closest to the leading edge (max distance from the TE along the
    curve's own chord -- robust to camber/asymmetry, doesn't assume
    y=0 marks the LE), resamples each surface (upper: TE->LE, lower:
    LE->TE) with Vinokur stretching using its OWN point count (n_upper,
    n_lower -- independently choosable, e.g. for a cambered/asymmetric
    section that needs more resolution on one side) but the SAME
    le_spacing/te_spacing target on both, so the two surfaces still meet
    the LE and TE with matching spacing despite differing point counts.
    n_upper/n_lower each include BOTH their own shared endpoint (LE or
    TE) -- output is one closed loop of n_upper + n_lower - 1 points (LE
    shared between the two, TE counted once at each end per the usual
    closed-loop convention).
    """
    te = pts[0]
    dist_from_te = np.linalg.norm(pts - te[None, :], axis=1)
    i_le = int(np.argmax(dist_from_te))

    upper = pts[: i_le + 1]       # TE -> LE
    lower = pts[i_le:]            # LE -> TE (pts[-1] == pts[0] == TE, closed)

    upper_r = resample_segment(upper, n_upper, te_spacing, le_spacing)
    lower_r = resample_segment(lower, n_lower, le_spacing, te_spacing)

    # upper_r's last point == lower_r's first point (both the LE) --
    # drop the duplicate when concatenating.
    return np.concatenate([upper_r, lower_r[1:]], axis=0)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--airfoil-file", required=True)
    parser.add_argument("--output-file", required=True)
    parser.add_argument("--n-upper", type=int, required=True, help="Points on the upper surface (TE to LE, inclusive)")
    parser.add_argument("--n-lower", type=int, required=True, help="Points on the lower surface (LE to TE, inclusive)")
    parser.add_argument("--le-spacing", type=float, required=True)
    parser.add_argument("--te-spacing", type=float, required=True)
    args = parser.parse_args()

    pts = load_airfoil_curve(args.airfoil_file)
    print(f"Input: {pts.shape[0]} points, closed-loop gap = {np.linalg.norm(pts[0]-pts[-1]):.3e}")

    new_pts = resample_airfoil(pts, args.n_upper, args.n_lower, args.le_spacing, args.te_spacing)
    print(f"Output: {new_pts.shape[0]} points, closed-loop gap = {np.linalg.norm(new_pts[0]-new_pts[-1]):.3e}")

    d = np.linalg.norm(np.diff(new_pts, axis=0), axis=1)
    i_le = np.argmax(np.linalg.norm(new_pts - new_pts[0][None, :], axis=1))
    print(f"Actual spacing: first(TE)={d[0]:.6e} (target {args.te_spacing:.3e}), "
          f"last(TE)={d[-1]:.6e} (target {args.te_spacing:.3e}), "
          f"near-LE={d[i_le-1]:.6e}/{d[i_le]:.6e} (target {args.le_spacing:.3e})")

    with open(args.output_file, "w") as f:
        f.write("resampled airfoil\n")
        for x, y in new_pts:
            f.write(f"{x:.7f}  {y:.7E}\n")
    print(f"Wrote {args.output_file}")
