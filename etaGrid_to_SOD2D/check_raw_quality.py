"""Builds the raw (untagged, unelevated) mesh and reports minSICN -- the
standing process check before ever submitting/running order-elevation."""
import sys
from pathlib import Path

import gmsh

from build_wake import build_full_mesh
from run_pipeline import collect_boundary_faces, write_raw_msh

work_dir = Path(sys.argv[1])
cross_section_file = sys.argv[2]
airfoil_file = sys.argv[3]
wake_s0 = float(sys.argv[4])
wake_total = float(sys.argv[5])
wake_layers = int(sys.argv[6])
basename = sys.argv[7] if len(sys.argv) > 7 else "etagrid"

work_dir.mkdir(parents=True, exist_ok=True)

print("[1/3] Building swept body + wake mesh...")
mesh = build_full_mesh(cross_section_file, airfoil_file, wake_s0, wake_total, wake_layers)

print("[2/3] Classifying boundary faces...")
boundary_faces = collect_boundary_faces(mesh)

raw_msh = work_dir / f"{basename}_raw.msh"
print(f"[3/3] Writing raw mesh: {raw_msh}")
entity_tags = write_raw_msh(mesh["coords"], mesh["hex_conn"], mesh["prism_conn"], boundary_faces, raw_msh)
print(f"  entity tags: {entity_tags}")

gmsh.initialize()
try:
    gmsh.open(str(raw_msh))
    etypes, etags, _ = gmsh.model.mesh.getElements(dim=3)
    all_neg = 0
    all_min = 1e9
    for et, tags in zip(etypes, etags):
        q = gmsh.model.mesh.getElementQualities(tags, "minSICN")
        neg = (q < 0).sum()
        all_neg += neg
        all_min = min(all_min, q.min())
        print(f"  type={et} n={len(tags)} minSICN={q.min():.6f} maxSICN={q.max():.6f} negative={neg}")
    print(f"TOTAL negative SICN elements: {all_neg} (overall min={all_min:.6f})")
finally:
    gmsh.finalize()
