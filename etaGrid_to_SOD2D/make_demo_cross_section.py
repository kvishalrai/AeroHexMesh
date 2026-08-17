"""Builds a clean, structured (transfinite) YZ cross-section mesh with
the same general domain shape as etaGrid/airfoil.msh (Y = wall-normal,
growing from 0 with geometric clustering near the wall, out to a far
Y-extent; Z = periodic/spanwise, uniform) -- but guaranteed well-formed
(no degenerate/inconsistently-wound elements), for isolating whether the
sweep/extrusion algorithm itself is correct, independent of the real
mesh's own quality (see project notes: the real mesh's negative-SICN
elements traced to ~40 degenerate quads from its own quad-subdivision
step, NOT the sweep algorithm -- this demo mesh has no such quads by
construction, since transfinite+recombine always yields a clean
structured grid).
"""
import argparse

import gmsh


def build_demo_cross_section(ny, nz, y_max, z_max, y_progression, output_file):
    gmsh.initialize()
    try:
        gmsh.model.add("demo_cross_section")

        p1 = gmsh.model.geo.addPoint(0, 0, 0)
        p2 = gmsh.model.geo.addPoint(0, y_max, 0)
        p3 = gmsh.model.geo.addPoint(0, y_max, z_max)
        p4 = gmsh.model.geo.addPoint(0, 0, z_max)

        l1 = gmsh.model.geo.addLine(p1, p2)  # y: 0 -> y_max, at z=0
        l2 = gmsh.model.geo.addLine(p2, p3)  # z: 0 -> z_max, at y=y_max
        l3 = gmsh.model.geo.addLine(p3, p4)  # y: y_max -> 0, at z=z_max
        l4 = gmsh.model.geo.addLine(p4, p1)  # z: z_max -> 0, at y=0

        loop = gmsh.model.geo.addCurveLoop([l1, l2, l3, l4])
        surf = gmsh.model.geo.addPlaneSurface([loop])

        # Fine spacing near y=0 (the wall) on BOTH z=0 (l1) and z=z_max
        # (l3) edges -- l1 runs 0->y_max so a positive progression puts
        # small steps at its start; l3 runs y_max->0 so it needs the
        # NEGATED progression to also put small steps at ITS end (y=0).
        gmsh.model.geo.mesh.setTransfiniteCurve(l1, ny, "Progression", y_progression)
        gmsh.model.geo.mesh.setTransfiniteCurve(l3, ny, "Progression", -y_progression)
        gmsh.model.geo.mesh.setTransfiniteCurve(l2, nz)
        gmsh.model.geo.mesh.setTransfiniteCurve(l4, nz)

        gmsh.model.geo.mesh.setTransfiniteSurface(surf, "Left", [p1, p2, p3, p4])
        gmsh.model.geo.mesh.setRecombine(2, surf)

        gmsh.model.geo.synchronize()
        gmsh.model.mesh.generate(2)

        node_tags, _, _ = gmsh.model.mesh.getNodes()
        etypes, etags, _ = gmsh.model.mesh.getElements(dim=2)
        print(f"Demo cross-section: {len(node_tags)} nodes, "
              f"{sum(len(t) for t in etags)} elements (types: {list(etypes)})")

        gmsh.option.setNumber("Mesh.MshFileVersion", 4.1)
        gmsh.write(output_file)
    finally:
        gmsh.finalize()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--ny", type=int, default=33)
    parser.add_argument("--nz", type=int, default=17)
    parser.add_argument("--y-max", type=float, default=20.0)
    parser.add_argument("--z-max", type=float, default=0.2)
    parser.add_argument("--y-progression", type=float, default=1.2)
    parser.add_argument("--output-file", required=True)
    args = parser.parse_args()

    build_demo_cross_section(args.ny, args.nz, args.y_max, args.z_max, args.y_progression, args.output_file)
    print(f"Wrote {args.output_file}")
