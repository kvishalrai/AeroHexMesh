"""Orchestrates the full etaGrid -> SOD2D pipeline, mirroring
pyHyp_to_SOD2D/convert_to_sod2d.py's own structure (kept as a separate,
self-contained framework per project notes -- this pipeline needs a
genuine periodic boundary, which that one doesn't support):

  1. sweep_mesh.py + build_wake.py -- build the swept body wrap + C-grid
     wake (see those modules' own docstrings)
  2. Classify boundary faces (wall/inlet/outlet/periodic) and write a raw
     discrete Gmsh mesh, MSH 2.2, with one entity per boundary group plus
     the volume -- boundary_conditions.py
  3. A companion .geo script assigns Physical Surface/Volume tags
     (WALL=1, INLET=2, OUTLET=3, PERIODIC=4, VOLUME=109 -- the exact same
     convention as Construct2D_to_SOD2D/linear_mesh/mesh_extrusion.py's
     own write_geo_file), sets up native Gmsh periodicity via
     'Periodic Surface ... Translate' (verified this mesh's Z=0/Z=zmax
     boundaries correspond exactly -- see project notes), and order-
     elevates via Mesh.ElementOrder -- then `gmsh file.geo -0` produces
     the final tagged mesh.
  4. gmsh2sod2d.py (vendored, reused directly) -- Gmsh mesh -> SOD2D HDF5
  5. tool_meshConversorPar (same submodule) -- partitions for parallel run
"""
import argparse
import subprocess
import sys
from pathlib import Path

import numpy as np
import gmsh

from build_wake import build_full_mesh
from boundary_conditions import build_boundary_faces_for_stations

SCRIPT_DIR = Path(__file__).resolve().parent
SOD2D_TOOLS = SCRIPT_DIR.parent / "Construct2D_to_SOD2D" / "linear_mesh" / "sod2d_tools"
GMSH2SOD2D = SOD2D_TOOLS / "gmsh2sod2d.py"
TOOL_MESHCONVERSORPAR = SOD2D_TOOLS / "tool_meshConversorPar"

WALL_ID, INLET_ID, OUTLET_ID, PERIODIC_ID, VOLUME_ID = 1, 2, 3, 4, 109

GEO_TEMPLATE = """\
Mesh.MshFileVersion = 2.2;

Merge "{raw_msh}";

Delete Physicals;

Physical Surface("WALL",{wall_id}) = {{{wall_tag}}};
Physical Surface("INLET",{inlet_id}) = {{{inlet_tag}}};
Physical Surface("OUTLET",{outlet_id}) = {{{outlet_tag}}};
Physical Surface("Periodic",{periodic_id}) = {{{per0_tag},{per1_tag}}};
Physical Volume("VolumeCode",{volume_id}) = {{{volume_tag}}};

Mesh.ElementOrder = {porder};

Mesh 3;

Recombine Volume {{{volume_tag}}};

Periodic Surface {{{per1_tag}}} = {{{per0_tag}}} Translate {{0,0,{span_z}}};

Save "{final_msh}";
"""


def collect_boundary_faces(mesh):
    body_faces = build_boundary_faces_for_stations(
        mesh["quad_conn"], mesh["tri_conn"], mesh["yz"], mesh["body_global_idx"], include_wall=True,
    )
    upper_faces = build_boundary_faces_for_stations(
        mesh["quad_conn"], mesh["tri_conn"], mesh["yz"], mesh["upper_global_idx"], include_wall=False,
    )
    lower_faces = build_boundary_faces_for_stations(
        mesh["quad_conn"], mesh["tri_conn"], mesh["yz"], mesh["lower_global_idx"], include_wall=False,
    )

    def cat(key, *groups):
        parts = [g[key] for g in groups if key in g and g[key].shape[0]]
        return np.concatenate(parts, axis=0) if parts else np.zeros((0, 4), dtype=np.int64)

    wall = body_faces["wall"]
    inlet = cat("inlet", body_faces, upper_faces, lower_faces)
    periodic0 = cat("periodic0", body_faces, upper_faces, lower_faces)
    periodic1 = cat("periodic1", body_faces, upper_faces, lower_faces)

    # Outlet: the actual 2D cap elements (quad+tri) at the FINAL wake
    # station of each side -- the true downstream-most boundary, not an
    # edge extrusion.
    upper_last, lower_last = mesh["upper_global_idx"][-1], mesh["lower_global_idx"][-1]
    outlet_quads = np.concatenate([upper_last[mesh["quad_conn"]], lower_last[mesh["quad_conn"]]], axis=0) \
        if mesh["quad_conn"].shape[0] else np.zeros((0, 4), dtype=np.int64)
    outlet_tris = np.concatenate([upper_last[mesh["tri_conn"]], lower_last[mesh["tri_conn"]]], axis=0) \
        if mesh["tri_conn"].shape[0] else np.zeros((0, 3), dtype=np.int64)

    print(f"Boundary faces: wall={wall.shape[0]} inlet={inlet.shape[0]} "
          f"periodic0={periodic0.shape[0]} periodic1={periodic1.shape[0]} "
          f"outlet={outlet_quads.shape[0]}(quad)+{outlet_tris.shape[0]}(tri)")

    return {
        "wall": (wall, np.zeros((0, 3), dtype=np.int64)),
        "inlet": (inlet, np.zeros((0, 3), dtype=np.int64)),
        "periodic0": (periodic0, np.zeros((0, 3), dtype=np.int64)),
        "periodic1": (periodic1, np.zeros((0, 3), dtype=np.int64)),
        "outlet": (outlet_quads, outlet_tris),
    }


def write_raw_msh(coords, hex_conn, prism_conn, boundary_faces, output_file):
    """Writes one discrete entity per boundary group (2D) plus the volume
    (3D), all referencing the SAME global node set -- returns
    {group_name: entity_tag} so the caller can build the .geo script's
    Physical Surface/Volume references."""
    gmsh.initialize()
    try:
        gmsh.model.add("etaGrid")
        node_tags = np.arange(1, coords.shape[0] + 1, dtype=np.int64)

        vol_tag = gmsh.model.addDiscreteEntity(3)
        gmsh.model.mesh.addNodes(3, vol_tag, node_tags, coords.flatten())
        elem_types, elem_tags_list, elem_nodes_list = [], [], []
        next_tag = 1
        if hex_conn.shape[0]:
            n = hex_conn.shape[0]
            elem_types.append(5)
            elem_tags_list.append(np.arange(next_tag, next_tag + n, dtype=np.int64))
            elem_nodes_list.append((hex_conn + 1).flatten())
            next_tag += n
        if prism_conn.shape[0]:
            n = prism_conn.shape[0]
            elem_types.append(6)
            elem_tags_list.append(np.arange(next_tag, next_tag + n, dtype=np.int64))
            elem_nodes_list.append((prism_conn + 1).flatten())
            next_tag += n
        gmsh.model.mesh.addElements(3, vol_tag, elem_types, elem_tags_list, elem_nodes_list)

        entity_tags = {"volume": vol_tag}
        for name, (quads, tris) in boundary_faces.items():
            tag = gmsh.model.addDiscreteEntity(2)
            entity_tags[name] = tag
            etypes, etags_list, enodes_list = [], [], []
            if quads.shape[0]:
                n = quads.shape[0]
                etypes.append(3)
                etags_list.append(np.arange(next_tag, next_tag + n, dtype=np.int64))
                enodes_list.append((quads + 1).flatten())
                next_tag += n
            if tris.shape[0]:
                n = tris.shape[0]
                etypes.append(2)
                etags_list.append(np.arange(next_tag, next_tag + n, dtype=np.int64))
                enodes_list.append((tris + 1).flatten())
                next_tag += n
            if etypes:
                gmsh.model.mesh.addElements(2, tag, etypes, etags_list, enodes_list)

        gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
        gmsh.write(str(output_file))
    finally:
        gmsh.finalize()

    return entity_tags


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--cross-section-file", required=True)
    parser.add_argument("--airfoil-file", required=True)
    parser.add_argument("--work-dir", required=True)
    parser.add_argument("--basename", default="etagrid")
    parser.add_argument("--wake-s0", type=float, default=0.005)
    parser.add_argument("--wake-total", type=float, default=20.0)
    parser.add_argument("--wake-layers", type=int, default=50)
    parser.add_argument("--porder", type=int, default=3)
    parser.add_argument("--num-partitions", type=int, default=1)
    parser.add_argument(
        "--stop-before-conversion", action="store_true",
        help="Stop after writing the tagged, order-elevated .msh (step 3) -- skip gmsh2sod2d.py "
        "and the partitioner, for inspecting the mesh first.",
    )
    args = parser.parse_args()

    work_dir = Path(args.work_dir)
    work_dir.mkdir(parents=True, exist_ok=True)

    print("[1/5] Building swept body + wake mesh...")
    mesh = build_full_mesh(
        args.cross_section_file, args.airfoil_file, args.wake_s0, args.wake_total, args.wake_layers,
    )

    print("[2/5] Classifying boundary faces...")
    boundary_faces = collect_boundary_faces(mesh)

    span_z = float(mesh["yz"][:, 1].max() - mesh["yz"][:, 1].min())

    raw_msh = work_dir / f"{args.basename}_raw.msh"
    print(f"[3/5] Writing raw tagged mesh: {raw_msh}")
    entity_tags = write_raw_msh(mesh["coords"], mesh["hex_conn"], mesh["prism_conn"], boundary_faces, raw_msh)
    print(f"  entity tags: {entity_tags}")

    final_msh = work_dir / f"{args.basename}_hi.msh"
    geo_file = work_dir / f"{args.basename}.geo"
    geo_file.write_text(GEO_TEMPLATE.format(
        raw_msh=raw_msh.name,
        wall_id=WALL_ID, inlet_id=INLET_ID, outlet_id=OUTLET_ID, periodic_id=PERIODIC_ID, volume_id=VOLUME_ID,
        wall_tag=entity_tags["wall"], inlet_tag=entity_tags["inlet"], outlet_tag=entity_tags["outlet"],
        per0_tag=entity_tags["periodic0"], per1_tag=entity_tags["periodic1"], volume_tag=entity_tags["volume"],
        porder=args.porder, span_z=span_z, final_msh=final_msh.name,
    ))
    print(f"[3/5] Wrote {geo_file}; running gmsh to tag + order-elevate to porder={args.porder}...")
    subprocess.run(["gmsh", geo_file.name, "-0"], check=True, cwd=work_dir)

    if args.stop_before_conversion:
        print(f"--stop-before-conversion set: stopping here. Tagged/elevated mesh at {final_msh}")
        return

    print(f"[4/5] Gmsh -> SOD2D HDF5 mesh: {work_dir / (args.basename + '_hi.h5')}")
    subprocess.run(
        [sys.executable, str(GMSH2SOD2D), f"{args.basename}_hi", "-r", str(args.porder), "-s", "500000",
         "-p", str(PERIODIC_ID)],
        check=True, cwd=work_dir,
    )

    print(f"[5/5] Partitioning into {args.num_partitions} piece(s)")
    import json
    input_json = {
        "gmsh_filePath": "", "gmsh_fileName": f"{args.basename}_hi",
        "mesh_h5_filePath": "", "mesh_h5_fileName": f"{args.basename}_hi",
        "num_partitions": args.num_partitions, "eval_mesh_quality": 0,
        "lineal_output": True, "uns_per_links": False,
    }
    (work_dir / "input.json").write_text(json.dumps(input_json, indent=4))
    print("(Run tool_meshConversorPar separately -- needs the SOD2D partitioner's own environment, see "
          "Construct2D_to_SOD2D/linear_mesh/README.md.)")


if __name__ == "__main__":
    main()
