"""Adds spanwise-sinusoidal waviness to an already spline-smoothed wall
surface: for each wall node (already sitting on the true airfoil curve,
per sweep_mesh.AirfoilCurve, after MeshElasticitySolver's own spline
correction), re-locate its arc-length position s on that SAME curve via
a nearest-point Newton search (robust since the node already sits almost
exactly on the curve -- this is the "re-evaluate the normals" step, done
fresh at the smoothed position rather than reusing any earlier value),
then displace it along curve.outward_normal(s) by

    amplitude_fraction * chord * sin(2*pi*wavenumber*z/z_max) * taper(s)

taper(s) = sin(pi*s/s_total)**2 -- zero (with zero slope, i.e. a smooth
ease, not a kink) at BOTH ends of the arc-length parametrization, which
are both the trailing edge (s=0 is the TE approached from the upper
surface, s=s_total is the TE approached from the lower surface -- see
AirfoilCurve's own docstring/construction in sweep_mesh.py), peaking at
the leading edge in between.
"""
import numpy as np

from sweep_mesh import AirfoilCurve


def nearest_s_newton(curve, points, s_guess, iters=25):
    """Newton/Gauss-Newton search for the arc-length s minimizing
    |curve.position(s) - point|^2, for a batch of points at once.
    s_guess: initial guess per point (e.g. from a coarse nearest-sample
    lookup). Converges in a handful of iterations given a decent guess
    and a point already close to the curve.
    """
    s = np.clip(s_guess.copy(), curve.s[0], curve.s[-1])
    for _ in range(iters):
        p = curve.position(s)
        t = curve.tangent(s)
        resid = p - points
        f = np.einsum("ij,ij->i", resid, t)
        fp = np.einsum("ij,ij->i", t, t)
        step = np.where(fp > 1e-300, f / fp, 0.0)
        s_new = np.clip(s - step, curve.s[0], curve.s[-1])
        if np.max(np.abs(s_new - s)) < 1e-13:
            s = s_new
            break
        s = s_new
    return s


def initial_s_guess(curve, points):
    """Coarse nearest-sample-point lookup among the curve's own known
    control points, used purely as a Newton starting guess."""
    d2 = np.sum((curve.pts[None, :, :] - points[:, None, :]) ** 2, axis=2)
    idx = np.argmin(d2, axis=1)
    return curve.s[idx]


def taper(curve, s):
    s_total = curve.s[-1]
    return np.sin(np.pi * s / s_total) ** 2


def wavy_wall_displacement(curve, points, z, amplitude_fraction, wavenumber, chord=1.0):
    """points: (N,3) current (already spline-smoothed) 3D wall positions.
    z: (N,) spanwise coordinate (exact, untouched by any prior step).
    Returns (new_points, s, normal) -- s and normal are the freshly
    re-evaluated arc-length/outward-normal at each point, for inspection.
    """
    xy = points[:, :2]
    s_guess = initial_s_guess(curve, xy)
    s = nearest_s_newton(curve, xy, s_guess)

    normal_2d = curve.outward_normal(s)  # (N,2)
    amp = amplitude_fraction * chord * np.sin(2 * np.pi * wavenumber * z / z.max()) * taper(curve, s)

    normal_3d = np.stack([normal_2d[:, 0], normal_2d[:, 1], np.zeros_like(s)], axis=1)
    new_points = points + amp[:, None] * normal_3d
    return new_points, s, normal_2d


if __name__ == "__main__":
    # Isolated sanity check against the same airfoil curve used by the
    # rest of the pipeline -- no HDF5/mesh involved yet.
    from sweep_mesh import load_airfoil_curve

    pts = load_airfoil_curve("naca0012_demo.dat")
    curve = AirfoilCurve(pts)
    s_total = curve.s[-1]
    print(f"chord (x range): {pts[:,0].min():.4f} to {pts[:,0].max():.4f}")
    print(f"arc-length total s_total={s_total:.6f}, LE roughly at s={s_total/2:.6f}")

    # Sample points EXACTLY on the curve at various s (simulating an
    # already spline-smoothed wall, no need to fabricate mesh noise).
    s_test = np.linspace(0.0, s_total, 9)
    p_test = curve.position(s_test)
    p_test_3d = np.concatenate([p_test, np.zeros((len(s_test), 1))], axis=1)
    z_test = np.linspace(0.0, 0.2, 9)  # pretend these came from 9 different spanwise stations

    new_pts, s_recovered, normals = wavy_wall_displacement(
        curve, p_test_3d, z_test, amplitude_fraction=0.002, wavenumber=5, chord=1.0
    )

    print(f"\n{'s (true)':>10} {'s (recovered)':>14} {'taper':>8} {'z':>6} {'disp_mag':>10}")
    for i in range(len(s_test)):
        disp_mag = np.linalg.norm(new_pts[i] - p_test_3d[i])
        sign = np.sign(np.dot(new_pts[i, :2] - p_test_3d[i, :2], normals[i]))
        print(f"{s_test[i]:10.5f} {s_recovered[i]:14.5f} {taper(curve, s_recovered[i]):8.4f} "
              f"{z_test[i]:6.3f} {sign*disp_mag:10.6f}")

    # Endpoints (both TE) should have ~zero displacement; check.
    assert abs(np.linalg.norm(new_pts[0] - p_test_3d[0])) < 1e-9, "TE (s=0) displacement not zero!"
    assert abs(np.linalg.norm(new_pts[-1] - p_test_3d[-1])) < 1e-9, "TE (s=s_total) displacement not zero!"
    print("\nOK: displacement vanishes at both TE ends.")
