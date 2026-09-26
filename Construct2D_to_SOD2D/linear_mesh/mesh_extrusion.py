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

def read_plot3d_2d_multiblock(filename):
    """Read a 2D Plot3D file, single- or multi-block.

    Construct2D's own single-block output has NO leading block-count line
    -- it starts directly with "ni nj". A multi-block file (e.g. a 2-block
    C-grid where the trailing edge is resolved as its own small block
    instead of a couple of extra rows inside one O-grid) starts with a
    block-count line, then one "ni nj" dims line per block, then the
    coordinate data: per block, all of x (flattened Fortran-order, i
    fastest) then all of y -- the same block-major/variable-within-block
    convention as p3d_to_gmsh.P3DfmtFile.load (the 3D equivalent format).

    Returns a list of (x2d, y2d) arrays, one per block (index 0 = block 1
    in 1-based NMF numbering).
    """
    with open(filename, 'r') as f:
        first_line = f.readline().split()
        if len(first_line) == 2:
            # Single-block: no leading count, this line IS the dims.
            dims = [tuple(map(int, first_line))]
        else:
            nblocks = int(first_line[0])
            dims = [tuple(map(int, f.readline().split())) for _ in range(nblocks)]

        def read_values(n):
            vals = []
            while len(vals) < n:
                vals += list(map(float, f.readline().split()))
            return np.array(vals)

        blocks = []
        for ni, nj in dims:
            npts = ni * nj
            x2d = read_values(npts).reshape((ni, nj), order='F')
            y2d = read_values(npts).reshape((ni, nj), order='F')
            blocks.append((x2d, y2d))
    return blocks


def read_plot3d_2d(filename):
    """Backward-compatible single-block reader: unchanged behavior/output
    for the single-block files this pipeline has always taken (identical
    lines consumed in the identical order as before), built on top of
    read_plot3d_2d_multiblock. Raises if the file actually has >1 block --
    use read_plot3d_2d_multiblock directly for those."""
    blocks = read_plot3d_2d_multiblock(filename)
    if len(blocks) != 1:
        raise ValueError(
            f"{filename}: expected a single-block 2D Plot3D file, got "
            f"{len(blocks)} blocks -- use read_plot3d_2d_multiblock instead"
        )
    return blocks[0]

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

def write_plot3d_ascii_multiblock(filename, blocks_3d):
    """blocks_3d: list of (x, y, z), each shape (ni, nj, nk). Always writes
    the leading block-count + per-block dims header expected by
    p3d_to_gmsh.P3DfmtFile.load (its 3D reader already requires this,
    even for one block), then per block all of x, then y, then z
    (flattened Fortran-order, i fastest) -- so a single-block call here
    reproduces write_plot3d_ascii's own output exactly."""
    with open(filename, 'w') as f:
        f.write(f"{len(blocks_3d)}\n")
        for x, y, z in blocks_3d:
            ni, nj, nk = x.shape
            f.write(f"{ni} {nj} {nk}\n")

        def write_array(arr):
            flat = arr.flatten(order='F')  # Fortran-style (i,j,k) → i fastest
            for i in range(0, len(flat), 5):
                f.write(" ".join(f"{v:.8e}" for v in flat[i:i+5]) + "\n")

        for x, y, z in blocks_3d:
            write_array(x)
            write_array(y)
            write_array(z)


def write_plot3d_ascii(filename, x, y, z):
    """Backward-compatible single-block writer: unchanged output, built on
    top of write_plot3d_ascii_multiblock."""
    write_plot3d_ascii_multiblock(filename, [(x, y, z)])

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

    # 2) Read 2D grid (single- or multi-block -- see
    #    read_plot3d_2d_multiblock's own docstring for the format
    #    difference; a single-block file naturally comes back as one
    #    block, so this is a strict generalization of the old single-
    #    block-only reader).
    blocks_2d = read_plot3d_2d_multiblock(output_file)

    # 3) Define extrusion in z (spanwise direction), same z_planes for
    #    every block (a multiblock case's blocks all share one spanwise
    #    extent -- see write_ext_nmf).
    z_plane = int(z_plane)
    z_len = float(z_len)
    z_planes = np.linspace(0.0, z_len, z_plane)  # You can set more layers if needed

    blocks_3d = [extrude_spanwise(x2d, y2d, z_planes) for x2d, y2d in blocks_2d]

    # 4) Write 3D extended Plot3D in work_dir with _ext suffix
    base = Path(filename).stem
    out_p3d = work_dir / f"{base}_ext.p3d"
    write_plot3d_ascii_multiblock(out_p3d, blocks_3d)
    
def _read_2d_nmf(nmf_path):
    """Parse a user-supplied 2D (KDIM=1) Construct2D Neutral Map File,
    single- or multi-block (the NMF format itself always carries a
    leading block-count + one dims line per block, regardless of block
    count -- unlike the .p3d grid file, see read_plot3d_2d_multiblock).

    Returns (block_dims, boundaries) where block_dims is a list of
    (idim, jmax) tuples, one per block (index 0 = block 1), and
    boundaries is a list of dicts with keys name, b1, f1, s1, e1, s2, e2,
    and, for ONE_TO_ONE/ONE-TO-ONE lines only, an additional 'pair' dict
    (b2, f2, s1, e1, s2, e2) describing the matched second range -- the
    standard NMF ONE_TO_ONE line packs two ranges (12 ints) instead of
    the usual one (6 ints). Every boundary from every block is returned
    (each record's own b1/pair.b2 says which); callers that only handle
    a single block can just look at boundaries with b1 == 1.
    """
    with open(nmf_path, 'r') as fp:
        non_comment = [l for l in fp if not l.lstrip().startswith('#')]
    non_blank = [l for l in non_comment if l.strip() != '']

    nblocks = int(non_blank[0].split()[0])
    block_dims = []
    for bn in range(nblocks):
        dims_tokens = non_blank[1 + bn].split()
        block_dims.append((int(dims_tokens[1]), int(dims_tokens[2])))
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

    return block_dims, boundaries


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
        they're merged the same way as OGRD's closure. A second,
        genuinely different CGRD shape is also accepted: a 2-block file
        where the blunt trailing edge is resolved as its own small block
        (its wall is the straight TE-closure line, its far end is
        FARFIELD) instead of a couple of extra rows inside one O-grid --
        there the two ONE_TO_ONE folds connect DIFFERENT blocks (not
        different faces of the same one), so the same-face check below
        only applies when both sides are the same block.

    Every ONE_TO_ONE fold (same-block wake-cut, O-grid closure, or a
    cross-block weld) is kept as a single ONE_TO_ONE line (remapped onto
    the extruded blocks). p3d_to_gmsh.py's GmshFile never emits a
    boundary face for a ONE_TO_ONE connection -- instead it aliases node
    ids across it (GmshFile._p3d_node_id_closed_i /
    GmshFile._build_wake_cut_alias for the same-block case,
    GmshFile._build_multiblock_weld for a cross-block one), so no
    Gmsh-side periodic surface or node merge is needed for any of them.

    The other (non-ONE_TO_ONE) boundaries -- the airfoil wall and every
    farfield-type face, on any block -- are classified as wall/inlet/
    outlet by p3d_to_gmsh.py itself (GmshFile._classify_inlet_outlet for
    the AoA-dependent farfield arc; VISCOUS is always wall, and the flat
    C-grid cut ends are always outlet, regardless of which block they're
    on), since that needs per-element geometry; this function only
    remaps face numbers and ranges, carrying the original boundary name
    through unchanged.

    Returns (idim, jmax, boundaries_2d): idim/jmax describe block 1 only
    (the only block for every case before the 2-block one above, and the
    main C-grid block there too), for callers that only handle a single
    block; boundaries_2d covers every block (each record's own b1/
    pair['b2'] says which).
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

    block_dims, boundaries_2d = _read_2d_nmf(src_nmf)

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

        same_block = rec['b1'] == pair['b2']
        if mesh_type == 'OGRD' and rec['f1'] == pair['f2']:
            raise ValueError(
                f"OGRD expects the ONE_TO_ONE closure to connect two "
                f"different faces (i=imin/i=imax), got face {rec['f1']} "
                f"matched to itself in {src_nmf}"
            )
        if mesh_type == 'CGRD' and same_block and rec['f1'] != pair['f2']:
            raise ValueError(
                f"CGRD expects the ONE_TO_ONE wake-cut to connect two "
                f"sub-ranges of the SAME face (the wall face), got "
                f"face {rec['f1']} matched to face {pair['f2']} in {src_nmf}"
            )
        if mesh_type == 'CGRD' and not same_block and (rec['f1'] not in (3, 4) or pair['f2'] not in (3, 4)):
            raise ValueError(
                f"CGRD expects a cross-block ONE_TO_ONE weld (block "
                f"{rec['b1']} <-> block {pair['b2']}) to connect two "
                f"j-const (2D face3/4) faces, got face {rec['f1']} "
                f"matched to face {pair['f2']} in {src_nmf}"
            )
        # The O-grid's full periodic closure, a same-block C-grid's
        # wake-cut fold, or a cross-block weld: keep as a single
        # ONE_TO_ONE line. p3d_to_gmsh.py's GmshFile aliases node ids
        # across it instead of emitting a boundary face, so it consumes
        # no physical-group id.
        lines.append(
            f"ONE_TO_ONE  {rec['b1']:5d}{f3d:5d}    {s1:5d}{e1:5d}    {s2:5d}{e2:5d}  "
            f"{pair['b2']:5d}{f3d_b:5d}    {s1b:5d}{e1b:5d}    {s2b:5d}{e2b:5d} FALSE\n"
        )

    with open(nmf_file, "w") as f:
        f.write("# ==================== Neutral Map File (extruded, derived from the user-supplied 2D NMF) ====================\n")
        f.write("# Block#   IDIM   JDIM   KDIM\n")
        # NeutralMapFile's own reader (p3d_to_gmsh.py, used for the
        # extruded NMF) skips this whole section by a FIXED line count --
        # one blank line, then exactly one line per block (back-to-back,
        # no blank between them), then one more blank line -- not by
        # recognizing dims lines individually, so the blank lines must
        # only bracket the whole block, never sit between two dims lines
        # (matches the actual Construct2D multiblock NMF format, e.g.
        # Construct2D/oat15_full_2block.nmf).
        f.write(f"       {len(block_dims)}\n\n")
        for bn, (bi, bj) in enumerate(block_dims, start=1):
            f.write(f"       {bn}    {bi:3d}    {bj:3d}      {z_planes:2d}\n")
        f.write("\n")
        f.write("# Type         B1  F1     S1   E1     S2   E2    B2  F2     S1   E1     S2   E2  Swap\n")
        for bn, (bi, bj) in enumerate(block_dims, start=1):
            # Every block's own spanwise end caps need a SYMMETRY-Y entry
            # -- for a single block this reproduces the previous fixed
            # "block 1" lines exactly; a second (or later) block's own
            # pair is additive, needed so its k=0/kmax faces get
            # boundary elements/periodic tagging too (see
            # GmshFile._gen_boundary_multiblock, which merges every
            # block's own k=0 entries into one shared physical group,
            # and likewise for kmax, rather than one group per block).
            f.write(f"SYMMETRY-Y    {bn}    1      1    {bi}      1   {bj}\n")
            f.write(f"SYMMETRY-Y    {bn}    2      1    {bi}      1   {bj}\n")
        f.writelines(lines)

    print(f"Wrote {nmf_file}")
    idim, jmax = block_dims[0]
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
