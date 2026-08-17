"""Sweep an unstructured (Y,Z) cross-section mesh around a closed 2D
airfoil curve to build a 3D boundary-layer-style volume mesh -- a
generalization of the wall-normal hyperbolic extrusion the rest of this
repo's pyHyp pipeline does, except here the "wall-normal direction" is a
whole pre-built unstructured quad(+tri) cross-section (encoding both the
near-wall growth AND a periodic spanwise extent) rather than a single 1D
stack of points, and the sweep path follows the airfoil's own curved
surface instead of marching outward from it.

Cross-section convention (see project notes -- confirmed with the user
against etaGrid/airfoil.msh's own node coordinates, all at x=0):
  Y = wall-normal offset (eta), growing away from the airfoil surface
  Z = spanwise/periodic direction, literal and unchanged by the sweep

Sweep: at each of N stations s_0..s_{N-1} along the airfoil's own closed
arc-length (using the airfoil file's OWN point distribution as the
station locations, preserving its LE/TE clustering), each cross-section
point (Y,Z) maps to 3D via
    (X,Y)_3d = airfoil_point(s) + Y * airfoil_outward_normal(s)
    Z_3d = Z
Consecutive stations' cross-sections are connected into one hex (from a
quad face) or prism (from a tri face) layer. This wraps the ENTIRE closed
airfoil loop but does NOT connect the last station back to the first --
station 0 and station N-1 are both physically at the trailing edge but
are kept as distinct mesh entities (the two "lips" of a C-grid wake cut),
ready for --wake-only's own straight downstream extrusion (see that
module).
"""
import numpy as np

try:
    import gmsh
except ImportError as e:
    raise ImportError("pip install gmsh (Python API) into this environment first") from e

from scipy.interpolate import CubicSpline


def load_cross_section(msh_file, recombine_to_quad=True):
    """Returns (yz (N,2), quad_conn (Nq,4) 0-based, tri_conn (Nt,3) 0-based)
    from a GMSH file whose nodes live in the Y-Z plane (x assumed ~0).

    recombine_to_quad: if the mesh has any triangles, they'd sweep into
    prism elements downstream -- which gmsh's own order-elevation can't
    handle above order 1/2 (confirmed empirically: "FaceClosureFull not
    implemented for prisms of order 3"). gmsh.model.mesh.recombine()
    (pairing adjacent triangles into one quad) was tried first but left
    this mesh's own 20 triangles stuck -- one region fails "Perfect
    Match in quadrangulation" for any recombination algorithm, since it
    operates on an already-existing imported mesh rather than generating
    fresh from a CAD surface. So when triangles are present, this
    instead does a full all-quad subdivision (Mesh.SubdivisionAlgorithm=1
    + gmsh.model.mesh.refine()) -- splits EVERY element (quad or tri) via
    edge-midpoints, guaranteeing zero triangles, at the cost of
    quadrupling the WHOLE cross-section's resolution (not just locally
    at the 20 problem elements) -- validated empirically to fully
    eliminate triangles on this mesh (7152 quads + 20 tris -> 28668
    quads, 0 tris)."""
    gmsh.initialize()
    try:
        gmsh.open(msh_file)

        if recombine_to_quad:
            has_tri = any(et == 2 for et in gmsh.model.mesh.getElements(dim=2)[0])
            if has_tri:
                gmsh.option.setNumber("Mesh.SubdivisionAlgorithm", 1)
                gmsh.model.mesh.refine()
                still_has_tri = any(et == 2 for et in gmsh.model.mesh.getElements(dim=2)[0])
                print(f"  Recombined to all-quad (subdivision): "
                      f"{'still has triangles -- unexpected' if still_has_tri else 'now 0 triangles'}")

        node_tags, node_coords_flat, _ = gmsh.model.mesh.getNodes()
        node_coords_flat = node_coords_flat.reshape(-1, 3)
        tag_to_idx = {int(t): i for i, t in enumerate(node_tags)}
        yz = node_coords_flat[:, 1:3].copy()

        etypes, etags, enodes = gmsh.model.mesh.getElements(dim=2)
        quad_blocks, tri_blocks = [], []
        for et, tags, nodes in zip(etypes, etags, enodes):
            nodes = np.asarray(nodes, dtype=np.int64).reshape(len(tags), -1)
            idx = np.vectorize(tag_to_idx.get)(nodes)
            if et == 3:  # 4-node quad
                quad_blocks.append(idx)
            elif et == 2:  # 3-node tri
                tri_blocks.append(idx)
    finally:
        gmsh.finalize()

    quad_conn = np.concatenate(quad_blocks, axis=0) if quad_blocks else np.zeros((0, 4), dtype=np.int64)
    tri_conn = np.concatenate(tri_blocks, axis=0) if tri_blocks else np.zeros((0, 3), dtype=np.int64)
    return yz, quad_conn, tri_conn


def load_airfoil_curve(dat_file):
    """Selig-format airfoil coordinate file (one optional header line,
    then 'x y' pairs, closed loop -- first and last point coincide for a
    sharp trailing edge). Returns (N,2) array."""
    pts = []
    with open(dat_file) as f:
        for line in f:
            parts = line.split()
            if len(parts) != 2:
                continue
            try:
                pts.append([float(parts[0]), float(parts[1])])
            except ValueError:
                continue
    return np.array(pts)


def arc_length_param(pts):
    d = np.linalg.norm(np.diff(pts, axis=0), axis=1)
    return np.concatenate([[0.0], np.cumsum(d)])


class AirfoilCurve:
    """Cubic-spline interpolant of a closed 2D airfoil curve, with
    position and outward-normal evaluable at arbitrary arc-length s."""

    def __init__(self, pts):
        self.pts = pts
        self.s = arc_length_param(pts)
        self.total = self.s[-1]
        self.spline = CubicSpline(self.s, pts, axis=0, bc_type="not-a-knot")

    def position(self, s):
        return self.spline(s)

    def tangent(self, s, eps=1e-6):
        s = np.atleast_1d(np.asarray(s, dtype=float))
        s_hi = np.clip(s + eps, 0.0, self.total)
        s_lo = np.clip(s - eps, 0.0, self.total)
        d = (self.spline(s_hi) - self.spline(s_lo)) / (s_hi - s_lo)[:, None]
        return d / np.linalg.norm(d, axis=1, keepdims=True)

    def outward_normal(self, s):
        """Rotate the tangent -90deg. Sign convention fixed empirically
        in this module's own __main__ check against a known point on the
        upper surface (normal should point toward +y, away from the
        airfoil body) -- see sanity_check_normal_sign()."""
        t = self.tangent(s)
        return np.stack([t[:, 1], -t[:, 0]], axis=1)


def sanity_check_normal_sign(curve):
    """Empirical check: at the arc-length station closest to the airfoil's
    own point of max thickness on the UPPER surface (y>0), the outward
    normal should point predominantly in +y. Returns True/False."""
    upper = curve.pts[curve.pts[:, 1] > 0]
    i_thick = np.argmax(upper[:, 1])
    p_thick = upper[i_thick]
    i_in_full = np.argmin(np.linalg.norm(curve.pts - p_thick, axis=1))
    s_thick = curve.s[i_in_full]
    n = curve.outward_normal(s_thick)[0]
    return n[1] > 0, n, p_thick


def sweep_cross_section(yz, quad_conn, tri_conn, curve, s_values):
    """Wraps the cross-section around the airfoil at each s in s_values
    (open wrap -- station -1 is NOT connected back to station 0). Returns
    (coords (n_stations*N, 3), hex_conn (0-based, 8/layer), prism_conn
    (0-based, 6/layer))."""
    N = yz.shape[0]
    n_stations = len(s_values)
    pos = curve.position(np.asarray(s_values))       # (n_stations,2)
    normal = curve.outward_normal(np.asarray(s_values))  # (n_stations,2)

    coords = np.zeros((n_stations, N, 3))
    coords[:, :, 0] = pos[:, 0:1] + yz[None, :, 0] * normal[:, 0:1]
    coords[:, :, 1] = pos[:, 1:2] + yz[None, :, 0] * normal[:, 1:2]
    coords[:, :, 2] = yz[None, :, 1]

    hex_layers = []
    prism_layers = []
    for k in range(n_stations - 1):
        base0, base1 = k * N, (k + 1) * N
        if quad_conn.shape[0]:
            hex_layers.append(np.concatenate([quad_conn + base0, quad_conn + base1], axis=1))
        if tri_conn.shape[0]:
            prism_layers.append(np.concatenate([tri_conn + base0, tri_conn + base1], axis=1))

    hex_conn = np.concatenate(hex_layers, axis=0) if hex_layers else np.zeros((0, 8), dtype=np.int64)
    prism_conn = np.concatenate(prism_layers, axis=0) if prism_layers else np.zeros((0, 6), dtype=np.int64)
    flat_coords = coords.reshape(-1, 3)
    hex_conn, prism_conn = ensure_positive_orientation(flat_coords, hex_conn, prism_conn)
    return flat_coords, hex_conn, prism_conn


def signed_hex_volume(coords, hex_conn):
    """Signed volume proxy (scalar triple product at corner 0, using
    edges to nodes 1, 3, 4 -- standard hex convention: nodes 0-3 one
    face, 4-7 directly 'above' in the same order) for a batch of hexes.
    Positive for a right-handed/valid element."""
    p0, p1, p3, p4 = coords[hex_conn[:, 0]], coords[hex_conn[:, 1]], coords[hex_conn[:, 3]], coords[hex_conn[:, 4]]
    e1, e2, e3 = p1 - p0, p3 - p0, p4 - p0
    return np.einsum("ij,ij->i", np.cross(e1, e2), e3)


def ensure_positive_orientation(coords, hex_conn, prism_conn):
    """Checks a sample of element volumes; if consistently negative,
    swaps the 'bottom 4'/'top 4' node halves (equivalent to reversing
    the local extrusion direction) to flip them positive.

    Why this is needed at all: confirmed empirically (see project notes)
    that this sweep construction can produce hexes whose "station_k then
    station_k+1" node ordering has EITHER handedness relative to the
    cross-section's own fixed 2D quad winding, depending on how the
    local sweep direction relates to that winding at each point -- the
    body wrap came out uniformly negative, the upper wake uniformly
    positive, the lower wake uniformly negative, each internally
    consistent but disagreeing with the others. So this is checked (and
    fixed if needed) separately per sweep call, not assumed globally."""
    if hex_conn.shape[0]:
        n = min(200, hex_conn.shape[0])
        vols = signed_hex_volume(coords, hex_conn[:n])
        if np.sum(vols < 0) > np.sum(vols > 0):
            hex_conn = hex_conn[:, [4, 5, 6, 7, 0, 1, 2, 3]]
    if prism_conn.shape[0]:
        n = min(200, prism_conn.shape[0])
        p0, p1, p2, p3 = (coords[prism_conn[:n, i]] for i in (0, 1, 2, 3))
        e1, e2, e3 = p1 - p0, p2 - p0, p3 - p0
        vols = np.einsum("ij,ij->i", np.cross(e1, e2), e3)
        if np.sum(vols < 0) > np.sum(vols > 0):
            prism_conn = prism_conn[:, [3, 4, 5, 0, 1, 2]]
    return hex_conn, prism_conn


def write_msh(coords, hex_conn, prism_conn, output_file):
    """Writes a standalone GMSH 4.1 file (one discrete entity per element
    type) via the gmsh API -- coords: (N,3) 0-based node array, hex_conn/
    prism_conn: 0-based (M,8)/(M,6) connectivity."""
    gmsh.initialize()
    try:
        gmsh.model.add("swept")
        tag = gmsh.model.addDiscreteEntity(3)
        node_tags = np.arange(1, coords.shape[0] + 1, dtype=np.int64)
        gmsh.model.mesh.addNodes(3, tag, node_tags, coords.flatten())

        elem_types, elem_tags_list, elem_nodes_list = [], [], []
        next_elem_tag = 1
        if hex_conn.shape[0]:
            elem_types.append(5)  # 8-node hex
            n = hex_conn.shape[0]
            elem_tags_list.append(np.arange(next_elem_tag, next_elem_tag + n, dtype=np.int64))
            elem_nodes_list.append((hex_conn + 1).flatten())
            next_elem_tag += n
        if prism_conn.shape[0]:
            elem_types.append(6)  # 6-node prism
            n = prism_conn.shape[0]
            elem_tags_list.append(np.arange(next_elem_tag, next_elem_tag + n, dtype=np.int64))
            elem_nodes_list.append((prism_conn + 1).flatten())
            next_elem_tag += n
        gmsh.model.mesh.addElements(3, tag, elem_types, elem_tags_list, elem_nodes_list)

        gmsh.write(output_file)
    finally:
        gmsh.finalize()


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--cross-section-file", required=True, help="GMSH .msh file, quad(+tri) mesh in the Y-Z plane")
    parser.add_argument("--airfoil-file", required=True, help="Selig-format airfoil coordinate file (closed loop)")
    parser.add_argument("--output-file", required=True, help="Output .msh (body wrap only, no wake yet)")
    args = parser.parse_args()

    yz, quad_conn, tri_conn = load_cross_section(args.cross_section_file)
    print(f"Cross-section: {yz.shape[0]} nodes, {quad_conn.shape[0]} quads, {tri_conn.shape[0]} tris")

    airfoil_pts = load_airfoil_curve(args.airfoil_file)
    print(f"Airfoil curve: {airfoil_pts.shape[0]} points, closed-loop gap = "
          f"{np.linalg.norm(airfoil_pts[0] - airfoil_pts[-1]):.3e}")
    curve = AirfoilCurve(airfoil_pts)
    print(f"Total arc length: {curve.total:.6f}")

    ok, n, p = sanity_check_normal_sign(curve)
    print(f"Outward-normal sign check at upper-surface max-thickness point {p}: normal={n}  {'OK' if ok else 'FLIPPED -- FIX SIGN'}")

    s_values = curve.s  # reuse the airfoil file's own point distribution as stations
    coords, hex_conn, prism_conn = sweep_cross_section(yz, quad_conn, tri_conn, curve, s_values)
    print(f"Wrap: {len(s_values)} stations, {coords.shape[0]} nodes, {hex_conn.shape[0]} hexes, {prism_conn.shape[0]} prisms")

    write_msh(coords, hex_conn, prism_conn, args.output_file)
    print(f"Wrote {args.output_file}")
