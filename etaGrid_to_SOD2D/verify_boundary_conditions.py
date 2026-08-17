"""Correct BC verification: since the cross-section is swept around a
curved airfoil, WALL/INLET are NOT flat planes in global (X,Y,Z) -- they
follow position(s) + Y_local*normal(s). So verify classification at the
LOCAL (Y,Z) cross-section level, directly from the same Python mesh-build
functions the pipeline uses (not by re-deriving geometry from the written
file), then separately confirm the raw .msh file's discrete-entity node
sets match those same indices exactly."""
import numpy as np
import gmsh

from build_wake import build_full_mesh
from run_pipeline import collect_boundary_faces, write_raw_msh

cross_section_file = "demo_cross_section.msh"
airfoil_file = "naca0012_demo.dat"
wake_s0, wake_total, wake_layers = 6.81e-3, 20.0, 10

mesh = build_full_mesh(cross_section_file, airfoil_file, wake_s0, wake_total, wake_layers)
yz = mesh["yz"]
y0_local, y1_local = yz[:, 0].min(), yz[:, 0].max()
z0_local, z1_local = yz[:, 1].min(), yz[:, 1].max()
N = yz.shape[0]

boundary_faces = collect_boundary_faces(mesh)

print("=== Local (Y,Z) cross-section classification check ===")
print(f"Cross-section local bbox: Y=[{y0_local:.4f},{y1_local:.4f}] Z=[{z0_local:.4f},{z1_local:.4f}]")


def local_row(global_node_idx):
    """Map a global node index back to its cross-section-local row index.
    Body/upper/lower each use global_idx arrays of shape (n_stations,N)
    where column j always corresponds to yz row j -- so global_idx % ...
    doesn't work directly (indices aren't contiguous mod N after merging),
    need an explicit reverse lookup table instead."""
    raise NotImplementedError


# Build a reverse lookup: global_node_idx -> yz row (0..N-1), from the
# same global_idx arrays the pipeline itself builds.
n_total = mesh["coords"].shape[0]
row_of = -np.ones(n_total, dtype=np.int64)
for gidx in (mesh["body_global_idx"], mesh["upper_global_idx"], mesh["lower_global_idx"]):
    for k in range(gidx.shape[0]):
        row_of[gidx[k]] = np.arange(N)

assert (row_of >= 0).all(), "some global node never got a local-row mapping!"

tol = 1e-9
for name, (quads, tris) in boundary_faces.items():
    if quads.shape[0] == 0 and tris.shape[0] == 0:
        continue
    idx = np.unique(np.concatenate([quads.flatten(), tris.flatten()]) if tris.shape[0] else quads.flatten())
    rows = row_of[idx]
    y_loc, z_loc = yz[rows, 0], yz[rows, 1]
    if name == "wall":
        bad = np.abs(y_loc) > tol
        print(f"wall: {len(idx)} nodes, local Y==0 check: {'OK' if not bad.any() else f'FAIL {bad.sum()}'}")
    elif name == "inlet":
        bad = np.abs(y_loc - y1_local) > tol
        print(f"inlet: {len(idx)} nodes, local Y==Y_max({y1_local:.4f}) check: {'OK' if not bad.any() else f'FAIL {bad.sum()}'}")
    elif name == "periodic0":
        bad = np.abs(z_loc - z0_local) > tol
        print(f"periodic0: {len(idx)} nodes, local Z==Z_min check: {'OK' if not bad.any() else f'FAIL {bad.sum()}'}")
    elif name == "periodic1":
        bad = np.abs(z_loc - z1_local) > tol
        print(f"periodic1: {len(idx)} nodes, local Z==Z_max check: {'OK' if not bad.any() else f'FAIL {bad.sum()}'}")
    elif name == "outlet":
        # outlet uses upper_global_idx[-1] and lower_global_idx[-1] directly by construction
        print(f"outlet: {len(idx)} nodes (constructed directly from final wake station, not a local-Y/Z classification)")

# Cross-check outlet is exactly the final wake station on both sides
outlet_expected = set(mesh["upper_global_idx"][-1].tolist()) | set(mesh["lower_global_idx"][-1].tolist())
outlet_actual = set(np.unique(boundary_faces["outlet"][0].flatten()).tolist())
print(f"\noutlet node-set matches final upper+lower wake station exactly: {outlet_expected == outlet_actual}")

# Also confirm WALL only touches body stations, never wake stations (a
# real correctness requirement: WALL must not appear in the wake).
wall_idx = set(np.unique(boundary_faces["wall"][0].flatten()).tolist())
wake_only_nodes = (set(mesh["upper_global_idx"][1:].flatten().tolist())
                    | set(mesh["lower_global_idx"][1:].flatten().tolist())) - set(mesh["body_global_idx"].flatten().tolist())
print(f"WALL/wake-only-node overlap: {len(wall_idx & wake_only_nodes)} (expect 0)")

print("\n=== Raw .msh discrete-entity cross-check ===")
raw_msh = "verify_raw_tmp.msh"
entity_tags = write_raw_msh(mesh["coords"], mesh["hex_conn"], mesh["prism_conn"], boundary_faces, raw_msh)
gmsh.initialize()
try:
    gmsh.open(raw_msh)
    for name in ["wall", "inlet", "periodic0", "periodic1", "outlet"]:
        tag = entity_tags[name]
        nt, _, _ = gmsh.model.mesh.getNodes(2, tag, includeBoundary=True)
        node_idx_0based = set((nt - 1).astype(np.int64).tolist())
        expected_idx = set(np.unique(
            np.concatenate([boundary_faces[name][0].flatten(), boundary_faces[name][1].flatten()])
            if boundary_faces[name][1].shape[0] else boundary_faces[name][0].flatten()
        ).tolist())
        print(f"{name}: raw-mesh entity node set == python boundary_faces node set: {node_idx_0based == expected_idx}")
finally:
    gmsh.finalize()

print("\n=== Periodic node-pair correspondence (direct, from tagged/elevated mesh) ===")
gmsh.initialize()
try:
    gmsh.open("run_demo_full/etagrid_demo_hi.msh")
    groups = gmsh.model.getPhysicalGroups(dim=2)
    per_entities = []
    for d, t in groups:
        if gmsh.model.getPhysicalName(d, t) == "Periodic":
            per_entities = gmsh.model.getEntitiesForPhysicalGroup(d, t)
    print(f"Periodic physical group entities: {list(per_entities)}")
    all_node_tags, all_coords, _ = gmsh.model.mesh.getNodes()
    coord_map = dict(zip(all_node_tags, all_coords.reshape(-1, 3)))
    for ent in per_entities:
        try:
            res = gmsh.model.mesh.getPeriodicNodes(2, ent, includeHighOrderNodes=True)
            tag_master, node_tags_slave, node_tags_master, affine = res
            print(f"  entity {ent}: masterTag={tag_master}, n_slave_nodes={len(node_tags_slave)}, "
                  f"n_master_nodes={len(node_tags_master)}, affine={list(affine) if len(affine) else '[]'}")
            if len(node_tags_slave):
                span_z = 0.2
                mism = 0
                for ns, nm in zip(node_tags_slave[:500], node_tags_master[:500]):
                    ps, pm = coord_map[ns], coord_map[nm]
                    ok = (np.allclose(ps, pm + np.array([0, 0, span_z]), atol=1e-6)
                          or np.allclose(ps, pm - np.array([0, 0, span_z]), atol=1e-6))
                    if not ok:
                        mism += 1
                print(f"  spot check on {min(500,len(node_tags_slave))} pairs: {mism} mismatched (expect 0)")
        except Exception as e:
            print(f"  entity {ent}: getPeriodicNodes raised: {e!r}")
finally:
    gmsh.finalize()
