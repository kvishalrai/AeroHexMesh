"""Spanwise extrusion of a 2D Plot3D mesh into 3D, NMF remapping, and Gmsh
.geo file generation.

The 2D Plot3D + Neutral Map File pair consumed here is produced externally
with Construct2D (this module does not run Construct2D itself); see the
project README for the expected inputs.
"""

import numpy as np
import json
from pathlib import Path

from p3d_to_gmsh import WALL_ID, INLET_ID, OUTLET_ID, VOLUME_ID

def read_plot3d_2d(filename):
    with open(filename, 'r') as f:
        ni, nj = map(int, f.readline().split())

        # Read x values
        x_vals = []
        while len(x_vals) < ni * nj:
            x_vals += list(map(float, f.readline().split()))
        x2d = np.array(x_vals).reshape((ni, nj), order='F')  # shape: (ni, nj)

        # Read y values
        y_vals = []
        while len(y_vals) < ni * nj:
            y_vals += list(map(float, f.readline().split()))
        y2d = np.array(y_vals).reshape((ni, nj), order='F')

    return x2d, y2d

def extrude_spanwise(x2d, y2d, z_vals):
    ni, nj = x2d.shape
    nk = len(z_vals)  # number of z (spanwise) layers

    x3d = np.zeros((ni, nj, nk))
    y3d = np.zeros_like(x3d)
    z3d = np.zeros_like(x3d)

    for k, zval in enumerate(z_vals):
        x3d[:, :, k] = x2d
        y3d[:, :, k] = y2d
        z3d[:, :, k] = zval

    return x3d, y3d, z3d  # shape: (ni, nj, nk)

def write_plot3d_ascii(filename, x, y, z):
    ni, nj, nk = x.shape
    with open(filename, 'w') as f:
        f.write("1\n")  # One block
        f.write(f"{ni} {nj} {nk}\n")

        def write_array(arr):
            flat = arr.flatten(order='F')  # Fortran-style (i,j,k) → k fastest
            for i in range(0, len(flat), 5):
                f.write(" ".join(f"{v:.8e}" for v in flat[i:i+5]) + "\n")

        write_array(x)
        write_array(y)
        write_array(z)

def sod2d_mesh(filename,work_dir,z_len,z_plane,MESH_TYPE):

    work_dir = Path(work_dir)

    # 1) Expect the user to have already provided a 2D Plot3D + NMF pair
    #    generated externally with Construct2D (this pipeline no longer
    #    invokes Construct2D itself).
    base = Path(filename).stem
    output_file = work_dir / f"{base}.p3d"
    nmf_file = work_dir / f"{base}.nmf"

    if not output_file.exists():
        raise FileNotFoundError(
            f"Expected a pre-generated 2D Plot3D file at {output_file}. "
            f"This pipeline no longer runs Construct2D; provide the 2D "
            f".p3d/.nmf pair yourself (mesh_type={MESH_TYPE})."
        )
    if not nmf_file.exists():
        raise FileNotFoundError(
            f"Expected a pre-generated 2D Neutral Map File at {nmf_file}. "
            f"This pipeline no longer runs Construct2D; provide the 2D "
            f".p3d/.nmf pair yourself (mesh_type={MESH_TYPE})."
        )

    # 2) Read 2D grid
    x2d, y2d = read_plot3d_2d(output_file)

    # 3) Define extrusion in z (spanwise direction)
    z_plane = int(z_plane)
    z_len = float(z_len)
    z_planes = np.linspace(0.0, z_len, z_plane)  # You can set more layers if needed

    x3d, y3d, z3d = extrude_spanwise(x2d, y2d, z_planes)

    # 4) Write 3D extended Plot3D in work_dir with _ext suffix
    base = Path(filename).stem
    out_p3d = work_dir / f"{base}_ext.p3d"
    write_plot3d_ascii(out_p3d, x3d, y3d, z3d)
    
def _read_2d_nmf(nmf_path):
    """Parse a user-supplied 2D (KDIM=1) Construct2D Neutral Map File.

    Returns (idim, jmax, boundaries) where boundaries is a list of dicts
    with keys name, b1, f1, s1, e1, s2, e2, and, for ONE_TO_ONE/
    ONE-TO-ONE lines only, an additional 'pair' dict (b2, f2, s1, e1, s2,
    e2) describing the matched second range -- the standard NMF ONE_TO_ONE
    line packs two ranges (12 ints) instead of the usual one (6 ints).
    """
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
        rec = {
            'name': name,
            'b1': nums[0], 'f1': nums[1],
            's1': nums[2], 'e1': nums[3], 's2': nums[4], 'e2': nums[5],
        }
        if is_one_to_one:
            rec['pair'] = {
                'b2': nums[6], 'f2': nums[7],
                's1': nums[8], 'e1': nums[9], 's2': nums[10], 'e2': nums[11],
            }
        boundaries.append(rec)

    return idim, jmax, boundaries


def _remap_2d_face_to_3d(face_2d, s1, e1, s2, e2, z_planes):
    """Remap a 2D (KDIM=1) NMF face + range to the corresponding face of
    the spanwise-extruded 3D block.

    Face numbering: 2D face1/2 (i=const) -> 3D face3/4 (i=const); 2D
    face3/4 (j=const) -> 3D face5/6 (j=const). 3D faces 1/2 (k=const) are
    the new spanwise end caps and don't exist in the 2D file.

    Axis order (matching GmshFile._gen_boundary in p3d_to_gmsh.py): i-const
    faces use (S1,E1)=J-range, (S2,E2)=K-range -- unchanged from the 2D
    file except the trivial K range is replaced by the full spanwise
    extent. j-const faces use a reversed (S1,E1)=K-range, (S2,E2)=I-range
    -- so the 2D file's I sub-range moves from (S1,E1) to (S2,E2).
    """
    face_3d = face_2d + 2
    if face_3d in (3, 4):
        return face_3d, s1, e1, 1, z_planes
    elif face_3d in (5, 6):
        return face_3d, 1, z_planes, s1, e1
    raise ValueError(f"Unexpected 2D face id {face_2d}")


def write_ext_nmf(airfoil_file, work_dir, z_planes, mesh_type):
    """Derive the extruded 3D Neutral Map File from the user-supplied 2D
    NMF, remapping face numbers/ranges for the spanwise extrusion. IDIM/JDIM
    are read directly from the 2D NMF's own header, not passed in. Returns
    (idim, jmax, boundaries_2d) so callers needing the 2D NMF's own contents
    don't have to re-parse it.

    mesh_type selects how the 2D file's j=jmin (wall) face is expected to
    be structured:
      - "OGRD": a single full-face closure (ONE_TO_ONE connecting the
        i=imin and i=imax faces), periodically closed in i.
      - "CGRD": the j=jmin face is split into the real airfoil wall
        (VISCOUS) plus a wake-cut ONE_TO_ONE fold matching two i
        sub-ranges of the SAME face. i is not periodically closed for a
        C-grid, but the wake-cut nodes ARE physically coincident, so
        they're merged the same way as OGRD's closure.

    Both cases are kept as a single ONE_TO_ONE line (remapped onto the
    extruded block). p3d_to_gmsh.py's GmshFile never emits a boundary face
    for a ONE_TO_ONE connection -- instead it aliases node ids across the
    fold/closure (GmshFile._p3d_node_id_closed_i /
    GmshFile._build_wake_cut_alias), so no Gmsh-side periodic surface or
    node merge is needed for either case.

    The other (non-ONE_TO_ONE) boundaries -- the airfoil wall and every
    farfield-type face -- are classified as wall/inlet/outlet by
    p3d_to_gmsh.py itself (GmshFile._classify_inlet_outlet), since that
    needs per-element geometry; this function only remaps face numbers
    and ranges, carrying the original boundary name through unchanged.
    """
    work_dir = Path(work_dir)
    mesh_type = mesh_type.upper()
    if mesh_type not in ('OGRD', 'CGRD'):
        raise ValueError(f"Unsupported mesh_type '{mesh_type}', expected 'OGRD' or 'CGRD'")

    base = Path(airfoil_file).stem
    src_nmf = work_dir / f"{base}.nmf"
    nmf_file = work_dir / f"{base}_ext.nmf"

    z_planes = int(z_planes)

    if not src_nmf.exists():
        raise FileNotFoundError(
            f"Expected the user-supplied 2D Neutral Map File at {src_nmf}"
        )

    idim, jmax, boundaries_2d = _read_2d_nmf(src_nmf)

    lines = []

    for rec in boundaries_2d:
        is_one_to_one = rec['name'].upper() in ('ONE_TO_ONE', 'ONE-TO-ONE')
        f3d, s1, e1, s2, e2 = _remap_2d_face_to_3d(
            rec['f1'], rec['s1'], rec['e1'], rec['s2'], rec['e2'], z_planes)

        if not is_one_to_one:
            lines.append(f"{rec['name']:<12s}{rec['b1']:5d}{f3d:5d}    {s1:5d}{e1:5d}    {s2:5d}{e2:5d}\n")
            continue

        pair = rec['pair']
        f3d_b, s1b, e1b, s2b, e2b = _remap_2d_face_to_3d(
            pair['f2'], pair['s1'], pair['e1'], pair['s2'], pair['e2'], z_planes)

        if mesh_type == 'OGRD' and rec['f1'] == pair['f2']:
            raise ValueError(
                f"OGRD expects the ONE_TO_ONE closure to connect two "
                f"different faces (i=imin/i=imax), got face {rec['f1']} "
                f"matched to itself in {src_nmf}"
            )
        if mesh_type == 'CGRD' and rec['f1'] != pair['f2']:
            raise ValueError(
                f"CGRD expects the ONE_TO_ONE wake-cut to connect two "
                f"sub-ranges of the SAME face (the wall face), got "
                f"face {rec['f1']} matched to face {pair['f2']} in {src_nmf}"
            )
        # Either the O-grid's full periodic closure or the C-grid's
        # wake-cut fold: keep as a single ONE_TO_ONE line. p3d_to_gmsh.py's
        # GmshFile aliases node ids across it instead of emitting a
        # boundary face, so it consumes no physical-group id.
        lines.append(
            f"ONE_TO_ONE  {rec['b1']:5d}{f3d:5d}    {s1:5d}{e1:5d}    {s2:5d}{e2:5d}  "
            f"{pair['b2']:5d}{f3d_b:5d}    {s1b:5d}{e1b:5d}    {s2b:5d}{e2b:5d} FALSE\n"
        )

    with open(nmf_file, "w") as f:
        f.write("# ==================== Neutral Map File (extruded, derived from the user-supplied 2D NMF) ====================\n")
        f.write("# Block#   IDIM   JDIM   KDIM\n")
        f.write("       1\n\n")
        f.write(f"       1    {idim:3d}    {jmax:3d}      {z_planes:2d}\n\n")
        f.write("# Type         B1  F1     S1   E1     S2   E2    B2  F2     S1   E1     S2   E2  Swap\n")
        f.write(f"SYMMETRY-Y    1    1      1    {idim}      1   {jmax}\n")
        f.write(f"SYMMETRY-Y    1    2      1    {idim}      1   {jmax}\n")
        f.writelines(lines)

    print(f"Wrote {nmf_file}")
    return idim, jmax, boundaries_2d


# Fixed physical-id convention for the final SOD2D mesh, used for every
# case (O-grid or C-grid): 1=WALL, 2=INLET, 3=OUTLET, 4=Periodic,
# 109=VolumeCode. WALL_ID/INLET_ID/OUTLET_ID/VOLUME_ID come from
# p3d_to_gmsh.py, which assigns them directly when building the raw mesh
# (see the comment there); PERIODIC_ID is only ever created here.
PERIODIC_ID = 4


def write_geo_file(airfoil_file, work_dir, span_z, porder, groups):
    """Write the periodic Gmsh geometry file.

    :groups: the (dim, id, name) list returned alongside the .msh file by
        p3d2gmsh() -- p3d_to_gmsh.py assigns 'wall'/'inlet'/'outlet' groups
        per-element from geometry (see GmshFile._classify_inlet_outlet),
        so the entity ids behind each final name can't be predicted ahead
        of time and must be read back from what was actually generated.
        The two spanwise SYMMETRY-Y faces are similarly looked up by name
        rather than assumed to be ids 2/3, since WALL/INLET/OUTLET/
        VolumeCode now claim the low fixed ids instead.
    """
    work_dir = Path(work_dir)
    base_name = airfoil_file.replace(".dat", "")

    def ids_for(name):
        return sorted(g[1] for g in groups if g[2] == name)

    wall_ids = ids_for('wall')
    inlet_ids = ids_for('inlet')
    outlet_ids = ids_for('outlet')
    # Face1 (k=0) is always consumed before face2 (k=kmax) in
    # GmshFile.consume(), and ids are handed out in that order, so the
    # lower id is always k=0 (the periodic master) and the higher is
    # k=kmax (the translated copy) -- matching the Translate direction
    # below.
    symmetry_ids = sorted(g[1] for g in groups if 'SYMMETRY-Y' in g[2])

    if not wall_ids:
        raise ValueError("No 'wall' physical group was generated by p3d2gmsh()")
    if not inlet_ids:
        raise ValueError(
            "No 'inlet' physical group was generated by p3d2gmsh() -- "
            "check angle_of_attack and mesh orientation"
        )
    if not outlet_ids:
        raise ValueError(
            "No 'outlet' physical group was generated by p3d2gmsh() -- "
            "check angle_of_attack and mesh orientation"
        )
    if len(symmetry_ids) != 2:
        raise ValueError(
            f"Expected exactly 2 SYMMETRY-Y (spanwise) groups, got {symmetry_ids}"
        )

    k0_id, kmax_id = symmetry_ids

    physical_surfaces = [
        f'Physical Surface("WALL",{WALL_ID}) = {{{",".join(map(str, wall_ids))}}};',
        f'Physical Surface("INLET",{INLET_ID}) = {{{",".join(map(str, inlet_ids))}}};',
        f'Physical Surface("OUTLET",{OUTLET_ID}) = {{{",".join(map(str, outlet_ids))}}};',
        f'Physical Surface("Periodic",{PERIODIC_ID}) = {{{k0_id},{kmax_id}}};',
    ]
    physical_block = '\n    '.join(physical_surfaces)

    geo_text = f"""
    span_z = {span_z};

    //+ Set the output mesh file version
    Mesh.MshFileVersion = 2.2;

    Merge "{base_name}.msh";

    // Make all physical groups to Zero
    Delete Physicals;

    // Create physical groups
    {physical_block}

    Physical Volume("VolumeCode",{VOLUME_ID}) = {{{VOLUME_ID}}};

    //+ Options controlling mesh generation/
    Mesh.ElementOrder = {porder};

    Mesh 3;

    Recombine Volume {{{VOLUME_ID}}};

    Periodic Surface {{{kmax_id}}} = {{{k0_id}}} Translate {{0,0,span_z}};

    Save "{base_name}_per.msh";
    """

    with open(work_dir / "airfoil_per.geo", "w") as f:
        f.write(geo_text)

    print("Wrote airfoil_per.geo")
    
def write_partition_input_json(base_name, work_dir, n_partitions):

    input_data = {
        "gmsh_filePath": "",
        "gmsh_fileName": f"{base_name}_per",

        "mesh_h5_filePath": "",
        "mesh_h5_fileName": f"{base_name}_per",

        "num_partitions": int(n_partitions),
        "eval_mesh_quality": 0,
        "lineal_output": True,
        "uns_per_links": False
    }

    with open(work_dir / "input.json", "w") as f:
        json.dump(input_data, f, indent=4)

    print("Wrote input.json")
