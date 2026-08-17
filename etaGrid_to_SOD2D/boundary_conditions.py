"""Boundary-face extraction and physical-group tagging for the swept C-grid
mesh, matching the SAME convention pyHyp_to_SOD2D/cgns_to_gmsh.py already
uses (so gmsh2sod2d.py downstream needs no changes): WALL=1, INLET=2,
OUTLET=3, VOLUME=109 -- extended here with PERIODIC=4 (that other
pipeline has no periodic boundary and reuses 4 for symmetry instead; this
one genuinely needs it, which is exactly why id 4 was left free).

Cross-section (Y,Z) boundary classification (see sweep_mesh.py's own
docstring for the Y=wall-normal/Z=spanwise-periodic convention):
  Y=0        -> WALL on the body wrap; the merged interior centerline in
                the wake (NOT tagged there -- it's not a physical wall)
  Y=Y_max    -> INLET (outer/far-field), on BOTH body wrap and wake
  Z=Z_min/max -> PERIODIC, on BOTH body wrap and wake
  the final wake station's own 2D cap elements -> OUTLET
"""
import numpy as np
from collections import defaultdict


def extract_boundary_edges(quad_conn, tri_conn):
    """Edges of the 2D cross-section mesh used by exactly one element
    (quad or tri) -- the cross-section's own outer boundary loop(s).
    Returns (M,2) int64 array of 0-based node-index pairs."""
    edge_count = defaultdict(int)

    def add(conn):
        n = conn.shape[1]
        for elem in conn:
            for i in range(n):
                a, b = int(elem[i]), int(elem[(i + 1) % n])
                edge_count[(min(a, b), max(a, b))] += 1

    if quad_conn.shape[0]:
        add(quad_conn)
    if tri_conn.shape[0]:
        add(tri_conn)
    edges = np.array([k for k, v in edge_count.items() if v == 1], dtype=np.int64)
    return edges


def classify_boundary_edges(boundary_edges, yz, tol=1e-9):
    """Splits boundary_edges into (wall, outer, periodic) groups by which
    side of the (Y,Z) domain both endpoints sit on. An edge whose two
    endpoints don't agree on a single side (shouldn't happen for a
    conforming boundary loop, but checked) is dropped with a warning."""
    y, z = yz[:, 0], yz[:, 1]
    y0, y1 = y.min(), y.max()
    z0, z1 = z.min(), z.max()

    def side(mask):
        a, b = boundary_edges[:, 0], boundary_edges[:, 1]
        return mask[a] & mask[b]

    wall_mask = side(np.abs(y) < tol)
    outer_mask = side(np.abs(y - y1) < tol)
    per0_mask = side(np.abs(z - z0) < tol)
    per1_mask = side(np.abs(z - z1) < tol)

    unclassified = ~(wall_mask | outer_mask | per0_mask | per1_mask)
    if unclassified.any():
        print(f"  WARNING: {unclassified.sum()} boundary edges did not classify to any known side")

    return (
        boundary_edges[wall_mask],
        boundary_edges[outer_mask],
        boundary_edges[per0_mask],
        boundary_edges[per1_mask],
    )


def edges_to_quads(edges, global_idx_lo, global_idx_hi):
    """Extrudes boundary edges between two stations (global 0-based node
    index arrays, one per station) into quad boundary faces (a@lo, b@lo,
    b@hi, a@hi)."""
    if edges.shape[0] == 0:
        return np.zeros((0, 4), dtype=np.int64)
    a, b = edges[:, 0], edges[:, 1]
    return np.stack([global_idx_lo[a], global_idx_lo[b], global_idx_hi[b], global_idx_hi[a]], axis=1)


def build_boundary_faces_for_stations(quad_conn, tri_conn, yz, global_idx, include_wall):
    """global_idx: (n_stations, N) 0-based global node indices for a
    sequence of stations (body wrap or one wake side). Returns dict
    {group_name: (n,4) quad faces} for wall/inlet/periodic0/periodic1,
    extruded across every consecutive station pair. include_wall=False
    skips the wall group entirely (used for the wake, whose Y=0 line is
    the merged interior centerline, not a physical wall)."""
    boundary_edges = extract_boundary_edges(quad_conn, tri_conn)
    wall_e, outer_e, per0_e, per1_e = classify_boundary_edges(boundary_edges, yz)

    groups = {"inlet": [], "periodic0": [], "periodic1": []}
    if include_wall:
        groups["wall"] = []

    n_stations = global_idx.shape[0]
    for k in range(n_stations - 1):
        lo, hi = global_idx[k], global_idx[k + 1]
        if include_wall:
            groups["wall"].append(edges_to_quads(wall_e, lo, hi))
        groups["inlet"].append(edges_to_quads(outer_e, lo, hi))
        groups["periodic0"].append(edges_to_quads(per0_e, lo, hi))
        groups["periodic1"].append(edges_to_quads(per1_e, lo, hi))

    return {name: np.concatenate(chunks, axis=0) if chunks else np.zeros((0, 4), dtype=np.int64)
            for name, chunks in groups.items()}
