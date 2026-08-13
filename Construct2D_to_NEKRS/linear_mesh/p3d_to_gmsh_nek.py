#!/usr/bin/env python3
"""
Convert a *native* (non-extruded) 2D Construct2D Plot3D mesh + Neutral Map
File into a linear (straight-sided) 2D Gmsh v2 mesh with wall/inlet/outlet
physical line groups, ready for order elevation (Gmsh's own `-order 2`,
see write_geo_file_2d() below) and then Nek5000's `gmsh2nek` utility.

gmsh2nek only accepts 2nd-order elements (QUAD8/QUAD9, LINE3) and doesn't
use physical *names* for anything -- only the physical group *tag* matters,
which becomes `boundaryID` in the .re2 file; real Nek5000 BC types
('W  ', 'v  ', 'O  ', ...) are assigned afterward in the case's .usr file
keyed off that boundaryID (see gmsh2nek's own README). So unlike the SOD2D
pipeline's p3d_to_gmsh.py, there's no need to bake specific BC-type
semantics into the physical group names here -- 'wall'/'inlet'/'outlet' are
just human-readable labels for the same fixed-id convention.

Node-id aliasing for the O-grid's periodic i-closure and the C-grid's
wake-cut fold, and the angle-of-attack-based inlet/outlet split on the
farfield arc, are carried over from p3d_to_gmsh.py
(Construct2D_to_SOD2D/linear_mesh/) with the same logic, just without the
spanwise/extrusion dimension -- this script works from the *native* 2D
.p3d/.nmf pair directly (Construct2D's own output), not the
extruded/remapped files run_pipeline.py's mesh_extrusion.py produces for
the SOD2D pipeline.

Native 2D Neutral Map File face numbers (KDIM=1): 1=Imin, 2=Imax,
3=Jmin (wall), 4=Jmax (farfield). An O-grid's Imin/Imax are ONE_TO_ONE
(periodic i-closure); a C-grid's are two independent FARFIELD (outlet)
faces, and its Jmin face is split into a VISCOUS (wall) range plus a
ONE_TO_ONE wake-cut fold (two i sub-ranges of the same face).
"""

import argparse
import math
import os.path

import numpy as np

WALL_ID = 1
INLET_ID = 2
OUTLET_ID = 3
DOMAIN_ID = 109


def read_plot3d_2d(filename):
    """Read a native (KDIM=1) Construct2D Plot3D file: 'idim jmax' header,
    then all x values, then all y values, Fortran-ordered."""
    with open(filename, 'r') as f:
        ni, nj = map(int, f.readline().split())

        x_vals = []
        while len(x_vals) < ni * nj:
            x_vals += list(map(float, f.readline().split()))
        x2d = np.array(x_vals).reshape((ni, nj), order='F')

        y_vals = []
        while len(y_vals) < ni * nj:
            y_vals += list(map(float, f.readline().split()))
        y2d = np.array(y_vals).reshape((ni, nj), order='F')

    return x2d, y2d


def read_2d_nmf(nmf_path):
    """Parse a native 2D (KDIM=1) Neutral Map File. Returns (idim, jmax,
    boundaries), boundaries a list of dicts with keys name, f1, s1, e1,
    and, for ONE_TO_ONE lines only, a 'pair' dict with f2, s1, e1."""
    with open(nmf_path, 'r') as fp:
        non_comment = [l for l in fp if not l.lstrip().startswith('#')]
    non_blank = [l for l in non_comment if l.strip() != '']

    nblocks = int(non_blank[0].split()[0])
    dims_tokens = non_blank[1].split()
    idim, jmax = int(dims_tokens[1]), int(dims_tokens[2])
    boundary_lines = non_blank[1 + nblocks:]

    boundaries = []
    for l in boundary_lines:
        tokens = l.replace('\\', ' ').split()
        if not tokens:
            continue
        name = tokens[0].strip('\'"')
        is_one_to_one = name.upper() in ('ONE_TO_ONE', 'ONE-TO-ONE')
        nums_needed = 12 if is_one_to_one else 6
        nums = list(map(int, tokens[1:1 + nums_needed]))
        rec = {'name': name, 'f1': nums[1], 's1': nums[2], 'e1': nums[3]}
        if is_one_to_one:
            rec['pair'] = {'f2': nums[7], 's1': nums[8], 'e1': nums[9]}
        boundaries.append(rec)

    return idim, jmax, boundaries


class GmshFile2DNek(object):
    """Builds a linear (straight-sided) 2D Gmsh v2 mesh: quadrangle (type
    3) domain elements, line (type 1) boundary elements. Order elevation
    to what gmsh2nek needs (QUAD8/9, LINE3) happens separately, via Gmsh
    itself -- see write_geo_file_2d()."""

    def __init__(self):
        self._nodes = []
        self._elements = []
        self._groups = [
            (2, DOMAIN_ID, 'domain'),
            (1, WALL_ID, 'wall'),
            (1, INLET_ID, 'inlet'),
            (1, OUTLET_ID, 'outlet'),
        ]
        self._element_id = 0
        self._mesh_type = None
        self._wake_cut_alias = None

    def _next_element_id(self):
        self._element_id += 1
        return self._element_id

    def _next_group_id(self):
        return max(self._groups, key=lambda g: g[1])[1] + 1

    def _get_or_create_group(self, name, dim=1):
        for g in self._groups:
            if g[2] == name:
                return g[1]
        gid = self._next_group_id()
        self._groups.append((dim, gid, name))
        return gid

    @staticmethod
    def _node_id(jmax, i, j):
        return 1 + j + jmax * i

    def _node_id_closed(self, jmax, i, j):
        """Apply O-grid periodic-i / C-grid wake-cut-fold node aliasing --
        same logic as p3d_to_gmsh.py's GmshFile._p3d_node_id_closed_i, just
        without the spanwise (k) dimension."""
        if self._mesh_type == 'OGRD' and i == self._idim - 1:
            i = 0
        elif (self._mesh_type == 'CGRD' and self._wake_cut_alias is not None
                and j == 0 and self._wake_cut_alias['i_lo_b'] <= i <= self._wake_cut_alias['i_hi_b']):
            i = self._wake_cut_alias['offset'] - i
        return self._node_id(jmax, i, j)

    @staticmethod
    def _build_wake_cut_alias(bdry):
        """bdry['f1']/bdry['pair']['f2'] must both be 3 (Jmin, the wall
        row); bdry['pair']'s (s1,e1) is the I-range of the fold's other
        side, listed in reverse order (s1 > e1), encoding the fold's
        orientation -- see p3d_to_gmsh.py's _build_wake_cut_alias."""
        f1, f2 = bdry['f1'], bdry['pair']['f2']
        if f1 != 3 or f2 != 3:
            raise ValueError(
                f"Expected the CGRD wake-cut ONE_TO_ONE to connect two "
                f"sub-ranges of face3 (J=0, the wall row), got faces {f1},{f2}"
            )
        s1_a, e1_a = bdry['s1'], bdry['e1']
        s1_b, e1_b = bdry['pair']['s1'], bdry['pair']['e1']
        return {
            'i_lo_a': s1_a - 1, 'i_hi_a': e1_a - 1,
            'i_lo_b': min(s1_b, e1_b) - 1, 'i_hi_b': max(s1_b, e1_b) - 1,
            'offset': s1_b - s1_a,
        }

    def _classify_inlet_outlet(self, x2d, y2d, i, angle_of_attack):
        """Same approach as p3d_to_gmsh.py's GmshFile._classify_inlet_outlet:
        approximate the local outward direction as wall->farfield at this
        i, compare against the free-stream direction from angle_of_attack."""
        jmax = x2d.shape[1] - 1
        dx = x2d[i, jmax] - x2d[i, 0]
        dy = y2d[i, jmax] - y2d[i, 0]
        norm = math.hypot(dx, dy)
        if norm > 0:
            dx, dy = dx / norm, dy / norm
        aoa_rad = math.radians(angle_of_attack)
        ux, uy = math.cos(aoa_rad), math.sin(aoa_rad)
        return 'inlet' if (dx * ux + dy * uy) < 0 else 'outlet'

    def consume(self, x2d, y2d, mesh_type, angle_of_attack, nmf_boundaries):
        idim, jmax = x2d.shape
        self._idim = idim
        self._mesh_type = mesh_type
        self._wake_cut_alias = None

        one_to_one = [b for b in nmf_boundaries if b['name'].upper() in ('ONE_TO_ONE', 'ONE-TO-ONE')]
        if mesh_type == 'CGRD' and one_to_one:
            self._wake_cut_alias = self._build_wake_cut_alias(one_to_one[0])

        # Nodes. An OGRD's i=idim-1 plane coincides with i=0 (see
        # _node_id_closed) and is deliberately not stored; a CGRD's
        # wake-cut-fold sub-range of j=0 is likewise aliased, not stored.
        idim_nodes = idim - 1 if mesh_type == 'OGRD' else idim
        for i in range(idim_nodes):
            for j in range(jmax):
                if (mesh_type == 'CGRD' and self._wake_cut_alias is not None
                        and j == 0 and self._wake_cut_alias['i_lo_b'] <= i <= self._wake_cut_alias['i_hi_b']):
                    continue
                nid = self._node_id(jmax, i, j)
                self._nodes.append((nid, x2d[i, j], y2d[i, j], 0.0))

        # Domain elements: one linear quad per (i,j) cell. j (radial) is the
        # OUTER loop and i (circumferential) the inner one deliberately --
        # tools/gmsh2nek preserves relative element order for elements of
        # the same physical surface, and Nek5000 case .usr files written
        # against an O-grid convention (e.g. smooth_geom0 in
        # smooth_2D/naca_set_e448.usr) index elements ring-by-ring
        # (element e+lev_offset, lev_offset=n_ogrid*(elev-1)): the first
        # NUMBER_ELEMENTS_X consecutive Nek elements must be one full
        # circumferential ring (constant j), not a radial spoke (constant
        # i) -- which is what looping i-outer would produce instead.
        for j in range(jmax - 1):
            for i in range(idim - 1):
                el_id = self._next_element_id()
                n1 = self._node_id_closed(jmax, i,     j)
                n2 = self._node_id_closed(jmax, i + 1, j)
                n3 = self._node_id_closed(jmax, i + 1, j + 1)
                n4 = self._node_id_closed(jmax, i,     j + 1)
                self._elements.append([el_id, 3, 2, DOMAIN_ID, DOMAIN_ID, n1, n2, n3, n4])

        # Boundaries.
        for bdry in nmf_boundaries:
            if bdry['name'].upper() in ('ONE_TO_ONE', 'ONE-TO-ONE'):
                continue  # realized by node aliasing above, no boundary element

            f1, s1, e1 = bdry['f1'], bdry['s1'], bdry['e1']

            if f1 == 3:  # Jmin: the wall (VISCOUS)
                gid = self._get_or_create_group('wall')
                for i in range(s1 - 1, e1 - 1):
                    el_id = self._next_element_id()
                    n1 = self._node_id_closed(jmax, i, 0)
                    n2 = self._node_id_closed(jmax, i + 1, 0)
                    self._elements.append([el_id, 1, 2, gid, gid, n1, n2])

            elif f1 == 4:  # Jmax: farfield arc, split inlet/outlet per element
                for i in range(s1 - 1, e1 - 1):
                    group_name = self._classify_inlet_outlet(x2d, y2d, i, angle_of_attack)
                    gid = self._get_or_create_group(group_name)
                    el_id = self._next_element_id()
                    n1 = self._node_id_closed(jmax, i, jmax - 1)
                    n2 = self._node_id_closed(jmax, i + 1, jmax - 1)
                    self._elements.append([el_id, 1, 2, gid, gid, n1, n2])

            elif f1 in (1, 2):  # Imin/Imax: a CGRD's two wake-end faces (always outlet)
                i_fixed = 0 if f1 == 1 else self._idim - 1
                gid = self._get_or_create_group('outlet')
                for j in range(s1 - 1, e1 - 1):
                    el_id = self._next_element_id()
                    n1 = self._node_id_closed(jmax, i_fixed, j)
                    n2 = self._node_id_closed(jmax, i_fixed, j + 1)
                    self._elements.append([el_id, 1, 2, gid, gid, n1, n2])

            else:
                raise ValueError(f"Unknown 2D NMF face identifier: {f1}")

    def save(self, filename):
        with open(filename, 'w') as out:
            out.write('$MeshFormat\n2.2 0 8\n$EndMeshFormat\n')
            out.write('$PhysicalNames\n%d\n' % len(self._groups))
            for grp in self._groups:
                out.write('%d %d "%s"\n' % grp)
            out.write('$EndPhysicalNames\n')
            out.write('$Nodes\n%d\n' % len(self._nodes))
            for node in self._nodes:
                out.write('%d %15.13e %15.13e %15.13e\n' % node)
            out.write('$EndNodes\n')
            out.write('$Elements\n%d\n' % len(self._elements))
            for el in self._elements:
                out.write('%s\n' % ' '.join(map(str, el)))
            out.write('$EndElements\n')

    @property
    def groups(self):
        return self._groups


def p3d2gmsh_nek(p3d_file, angle_of_attack, mesh_type, map_file=None, output_file=None):
    """Convert one native 2D Construct2D Plot3D file to a linear 2D Gmsh
    .msh file. Returns (output_file, groups) -- see GmshFile2DNek.consume
    for why inlet/outlet ids can't be predicted ahead of time."""
    if not os.path.exists(p3d_file):
        raise FileNotFoundError(f"Cannot open {p3d_file}")

    name, _ = os.path.splitext(p3d_file)
    if map_file is None:
        map_file = f"{name}.nmf"
    if output_file is None:
        output_file = f"{name}.msh"

    print(f"Reading P3D file : {p3d_file}")
    print(f"Reading NMF file : {map_file}")
    print(f"Writing MSH file : {output_file}")

    x2d, y2d = read_plot3d_2d(p3d_file)
    _, _, boundaries = read_2d_nmf(map_file)

    gmsh = GmshFile2DNek()
    gmsh.consume(x2d, y2d, mesh_type, angle_of_attack, boundaries)
    gmsh.save(output_file)

    return output_file, gmsh.groups


def write_geo_file_2d(msh_file, geo_file, groups, output_msh, order=2):
    """Write a .geo file that merges msh_file, re-declares physical groups
    at their fixed ids (Gmsh silently keeps a stale group if one is
    re-declared at an id another group already used -- see
    mesh_extrusion.write_geo_file's identical comment), elevates element
    order to `order` -- the QUAD8/9 + LINE3 gmsh2nek needs -- and saves the
    result to output_msh. Run `gmsh <geo_file> -0` to process it (mirrors
    run_pipeline.py's Gmsh invocation for the SOD2D pipeline)."""
    def ids_for(name):
        return sorted(g[1] for g in groups if g[2] == name)

    wall_ids = ids_for('wall')
    inlet_ids = ids_for('inlet')
    outlet_ids = ids_for('outlet')
    if not wall_ids:
        raise ValueError("No 'wall' physical group was generated by p3d2gmsh_nek()")
    if not inlet_ids:
        raise ValueError("No 'inlet' physical group was generated -- check angle_of_attack")
    if not outlet_ids:
        raise ValueError("No 'outlet' physical group was generated -- check angle_of_attack")

    physical_lines = [
        f'Physical Line("WALL",{WALL_ID}) = {{{",".join(map(str, wall_ids))}}};',
        f'Physical Line("INLET",{INLET_ID}) = {{{",".join(map(str, inlet_ids))}}};',
        f'Physical Line("OUTLET",{OUTLET_ID}) = {{{",".join(map(str, outlet_ids))}}};',
    ]
    physical_block = '\n    '.join(physical_lines)

    geo_text = f"""
    Mesh.MshFileVersion = 2.2;

    Merge "{msh_file}";

    Delete Physicals;

    {physical_block}
    Physical Surface("domain",{DOMAIN_ID}) = {{{DOMAIN_ID}}};

    Mesh.ElementOrder = {order};

    Mesh 2;

    Save "{output_msh}";
    """

    with open(geo_file, 'w') as f:
        f.write(geo_text)

    print(f"Wrote {geo_file}")


def main():
    parser = argparse.ArgumentParser(
        description="Convert a native 2D Construct2D Plot3D mesh into a linear Gmsh mesh for Nek5000's gmsh2nek",
        add_help=True,
    )
    parser.add_argument("files", nargs="+", help="Native 2D Plot3D files to convert")
    parser.add_argument("-m", "--map-file", nargs=1, help="Neutral Map File. If omitted, uses <filename>.nmf")
    parser.add_argument("-o", "--output-file", nargs=1, help="Output Gmsh .msh file. If omitted, uses <filename>.msh")
    parser.add_argument("--aoa", type=float, required=True, help="Angle of attack in degrees")
    parser.add_argument(
        "--mesh-type", choices=["OGRD", "CGRD"], required=True,
        help="Mesh topology: OGRD (periodically closed in i) or CGRD (wake-cut)",
    )
    args = parser.parse_args()

    for fn in args.files:
        map_file = args.map_file[0] if args.map_file is not None else None
        output_file = args.output_file[0] if args.output_file is not None else None
        p3d2gmsh_nek(
            p3d_file=fn,
            angle_of_attack=args.aoa,
            mesh_type=args.mesh_type,
            map_file=map_file,
            output_file=output_file,
        )


if __name__ == '__main__':
    main()
