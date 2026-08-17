"""Extends sweep_mesh.py's body wrap with a straight downstream (+X) wake
block, completing the C-grid topology: the wrap's two trailing-edge-
adjacent stations (station 0, approached from the upper surface, and the
last station, approached from the lower surface -- see sweep_mesh.py's
own docstring) each get extruded straight downstream, using their OWN
frozen (position, outward-normal) evaluated once at the trailing edge --
NOT continuing to follow the airfoil curve. Because Y=0 (the wall) always
maps to just the frozen position regardless of the normal direction, the
upper-wake and lower-wake sides' wall (Y=0) nodes are geometrically
IDENTICAL at every downstream station -- these get merged into shared
nodes (the wake centerline), while Y>0 nodes stay distinct between the
two sides (the genuine, physically-real gap between the upper and lower
wake sheets, widening gradually downstream due to the airfoil's own
small-but-nonzero trailing-edge angle). The wake's own x_offset=0 station
is identical to the body wrap's station 0 / last station, so it directly
reuses those nodes too -- the whole mesh (body + wake) ends up as one
node-sharing-consistent assembly, no explicit stitching pass needed.
"""
import numpy as np

from sweep_mesh import (
    load_cross_section, load_airfoil_curve, AirfoilCurve, sweep_cross_section, write_msh,
    ensure_positive_orientation,
)


def geometric_offsets(s0, total, n_layers, tol=1e-10, max_iter=200):
    """n_layers+1 offsets (0, s0, s0+s0*r, ...) via geometric growth,
    solving for the ratio r that reaches `total` in exactly n_layers
    steps from a first-step size of s0 (bisection on the standard
    geometric-series-sum equation, since it isn't solvable in closed
    form for arbitrary n)."""
    if n_layers * s0 >= total:
        raise ValueError(f"s0={s0} alone over {n_layers} layers already reaches/exceeds total={total}")

    def sum_at(r):
        if abs(r - 1.0) < 1e-12:
            return s0 * n_layers
        return s0 * (r**n_layers - 1.0) / (r - 1.0)

    r_lo, r_hi = 1.0, 2.0
    while sum_at(r_hi) < total:
        r_hi *= 1.5
    for _ in range(max_iter):
        r_mid = 0.5 * (r_lo + r_hi)
        if sum_at(r_mid) < total:
            r_lo = r_mid
        else:
            r_hi = r_mid
        if r_hi - r_lo < tol:
            break
    r = 0.5 * (r_lo + r_hi)

    offsets = [0.0]
    step = s0
    for _ in range(n_layers):
        offsets.append(offsets[-1] + step)
        step *= r
    return np.array(offsets), r


def sweep_wake_side(yz, quad_conn, tri_conn, p, n, x_offsets):
    """Straight +X extrusion of the cross-section from a FROZEN (p,n)
    (evaluated once at the trailing edge) -- same layered hex/prism
    connectivity pattern as sweep_cross_section, just varying x_offset
    instead of arc-length s. Returns (coords, hex_conn, prism_conn) with
    the same conventions as sweep_cross_section."""
    N = yz.shape[0]
    n_stations = len(x_offsets)
    coords = np.zeros((n_stations, N, 3))
    coords[:, :, 0] = p[0] + x_offsets[:, None] + yz[None, :, 0] * n[0]
    coords[:, :, 1] = p[1] + yz[None, :, 0] * n[1]
    coords[:, :, 2] = yz[None, :, 1]

    hex_layers, prism_layers = [], []
    for k in range(n_stations - 1):
        base0, base1 = k * N, (k + 1) * N
        if quad_conn.shape[0]:
            hex_layers.append(np.concatenate([quad_conn + base0, quad_conn + base1], axis=1))
        if tri_conn.shape[0]:
            prism_layers.append(np.concatenate([tri_conn + base0, tri_conn + base1], axis=1))
    hex_conn = np.concatenate(hex_layers, axis=0) if hex_layers else np.zeros((0, 8), dtype=np.int64)
    prism_conn = np.concatenate(prism_layers, axis=0) if prism_layers else np.zeros((0, 6), dtype=np.int64)
    return coords.reshape(-1, 3), hex_conn, prism_conn


def build_full_mesh(cross_section_file, airfoil_file, wake_s0, wake_total, wake_layers):
    yz, quad_conn, tri_conn = load_cross_section(cross_section_file)
    N = yz.shape[0]
    wall_idx = np.where(np.abs(yz[:, 0]) < 1e-12)[0]

    airfoil_pts = load_airfoil_curve(airfoil_file)
    curve = AirfoilCurve(airfoil_pts)
    s_values = curve.s
    n_body_stations = len(s_values)

    body_coords, body_hex, body_prism = sweep_cross_section(yz, quad_conn, tri_conn, curve, s_values)

    x_offsets, ratio = geometric_offsets(wake_s0, wake_total, wake_layers)
    print(f"Wake growth ratio solved: {ratio:.5f} ({wake_layers} layers, s0={wake_s0}, total={wake_total})")

    p0, n0 = curve.position(np.array([s_values[0]]))[0], curve.outward_normal(np.array([s_values[0]]))[0]
    pN, nN = curve.position(np.array([s_values[-1]]))[0], curve.outward_normal(np.array([s_values[-1]]))[0]

    upper_coords, upper_hex, upper_prism = sweep_wake_side(yz, quad_conn, tri_conn, p0, n0, x_offsets)
    lower_coords, lower_hex, lower_prism = sweep_wake_side(yz, quad_conn, tri_conn, pN, nN, x_offsets)

    n_wake_stations = len(x_offsets)

    # --- Node assembly with sharing ---
    # Global layout: [body nodes] + [upper-wake nodes for wake stations 1..end,
    # station 0 reuses body's own station-0 nodes] + [lower-wake NON-WALL
    # nodes for stations 1..end, station 0 reuses body's last-station nodes;
    # wall nodes at every wake station reuse the corresponding upper-wake
    # wall node instead of duplicating].
    n_body = body_coords.shape[0]
    global_coords = [body_coords]
    offset = n_body

    # upper-wake station 0 == body station 0 (indices 0..N-1); stations 1..end are new
    upper_global_idx = np.zeros((n_wake_stations, N), dtype=np.int64)
    upper_global_idx[0, :] = np.arange(0, N)
    for w in range(1, n_wake_stations):
        upper_global_idx[w, :] = np.arange(offset, offset + N)
        global_coords.append(upper_coords.reshape(n_wake_stations, N, 3)[w])
        offset += N

    # lower-wake station 0 == body's LAST station (indices (n_body_stations-1)*N .. +N)
    body_last_base = (n_body_stations - 1) * N
    lower_global_idx = np.zeros((n_wake_stations, N), dtype=np.int64)
    lower_global_idx[0, :] = np.arange(body_last_base, body_last_base + N)
    for w in range(1, n_wake_stations):
        row = np.empty(N, dtype=np.int64)
        row[wall_idx] = upper_global_idx[w, wall_idx]  # merge onto upper-wake's wall nodes
        non_wall = np.setdiff1d(np.arange(N), wall_idx, assume_unique=True)
        row[non_wall] = np.arange(offset, offset + len(non_wall))
        global_coords.append(lower_coords.reshape(n_wake_stations, N, 3)[w][non_wall])
        offset += len(non_wall)
        lower_global_idx[w, :] = row

    all_coords = np.concatenate(global_coords, axis=0)
    print(f"Total nodes after merge: {all_coords.shape[0]} "
          f"(vs {n_body + 2 * n_wake_stations * N} if fully duplicated)")

    # Rebuild wake hex/prism connectivity using the GLOBAL indices per station
    # (sweep_wake_side's own hex/prism arrays used LOCAL 0..2N-1 indices per
    # layer; re-derive per-layer using upper/lower_global_idx directly for
    # clarity and correctness rather than trying to offset the local arrays).
    def build_layers(quad_conn, tri_conn, global_idx):
        hex_layers, prism_layers = [], []
        for k in range(global_idx.shape[0] - 1):
            lo, hi = global_idx[k], global_idx[k + 1]
            if quad_conn.shape[0]:
                hex_layers.append(np.concatenate([lo[quad_conn], hi[quad_conn]], axis=1))
            if tri_conn.shape[0]:
                prism_layers.append(np.concatenate([lo[tri_conn], hi[tri_conn]], axis=1))
        hx = np.concatenate(hex_layers, axis=0) if hex_layers else np.zeros((0, 8), dtype=np.int64)
        pr = np.concatenate(prism_layers, axis=0) if prism_layers else np.zeros((0, 6), dtype=np.int64)
        return hx, pr

    upper_hex_g, upper_prism_g = build_layers(quad_conn, tri_conn, upper_global_idx)
    lower_hex_g, lower_prism_g = build_layers(quad_conn, tri_conn, lower_global_idx)

    # sweep_cross_section already self-corrects body's own orientation;
    # upper/lower wake connectivity was rebuilt here with global (merged)
    # indices using the same "lo-then-hi" pattern, so it needs the same
    # check -- confirmed empirically that upper-wake and lower-wake can
    # have OPPOSITE handedness from each other (and from the body), since
    # each depends on how the local sweep direction there relates to the
    # cross-section's own fixed 2D winding.
    upper_hex_g, upper_prism_g = ensure_positive_orientation(all_coords, upper_hex_g, upper_prism_g)
    lower_hex_g, lower_prism_g = ensure_positive_orientation(all_coords, lower_hex_g, lower_prism_g)

    all_hex = np.concatenate([body_hex, upper_hex_g, lower_hex_g], axis=0)
    all_prism = np.concatenate([body_prism, upper_prism_g, lower_prism_g], axis=0)
    print(f"Total elements: {all_hex.shape[0]} hex, {all_prism.shape[0]} prism")

    body_global_idx = np.arange(n_body_stations * N).reshape(n_body_stations, N)

    return {
        "coords": all_coords,
        "hex_conn": all_hex,
        "prism_conn": all_prism,
        "yz": yz,
        "quad_conn": quad_conn,
        "tri_conn": tri_conn,
        "body_global_idx": body_global_idx,
        "upper_global_idx": upper_global_idx,
        "lower_global_idx": lower_global_idx,
    }


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--cross-section-file", required=True)
    parser.add_argument("--airfoil-file", required=True)
    parser.add_argument("--output-file", required=True)
    parser.add_argument("--wake-s0", type=float, default=0.005, help="First wake layer thickness")
    parser.add_argument("--wake-total", type=float, default=20.0, help="Total wake extent in X")
    parser.add_argument("--wake-layers", type=int, default=50)
    args = parser.parse_args()

    mesh = build_full_mesh(
        args.cross_section_file, args.airfoil_file, args.wake_s0, args.wake_total, args.wake_layers,
    )
    write_msh(mesh["coords"], mesh["hex_conn"], mesh["prism_conn"], args.output_file)
    print(f"Wrote {args.output_file}")
