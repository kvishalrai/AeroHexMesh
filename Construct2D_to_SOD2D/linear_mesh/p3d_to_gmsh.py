#!/usr/bin/env python3
"""
Convert a (spanwise-extruded) Plot3D airfoil mesh + Neutral Map File into
Gmsh's MSH format, ready for the mesh_extrusion.write_geo_file /
sod2d_tools/gmsh2sod2d.py steps of the pipeline.

Handles both O-grid (periodically closed in i) and C-grid (wake-cut fold,
merged via node-id aliasing rather than a boundary face) topologies, and
classifies every farfield-type boundary point as wall/inlet/outlet using
the fixed physical-id convention defined below.

The Plot3D/NMF/Gmsh I/O classes (read_chunk, NeutralMapFile, P3DfmtFile,
GmshFile) are adapted from p3d2gmsh (https://github.com/mrklein/p3d2gmsh):

    The MIT License (MIT)
    Copyright (c) 2015 Alexey Matveichev

    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files (the "Software"), to deal
    in the Software without restriction, including without limitation the rights
    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    copies of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in all
    copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
    SOFTWARE.

The O/C-grid handling, wall/inlet/outlet classification, and periodic
node-id aliasing below are new work on top of that base.
"""

import sys
import argparse
import math
import os.path
import numpy as np

# Fixed physical-id convention for the final SOD2D mesh (shared with
# mesh_extrusion.write_geo_file, which adds "Periodic" = 4 on top of these):
# 1=WALL, 2=INLET, 3=OUTLET, 109=VolumeCode. These are assigned directly
# when the groups are first created (rather than auto-numbered and
# renumbered later), because re-declaring a Gmsh physical group at an id
# some OTHER group already used silently does nothing -- the old group
# lingers instead of being replaced, which is exactly what caused a
# duplicate volume zone ("mesh"=1 alongside "VolumeCode"=109) the first
# time this was tried.
WALL_ID = 1
INLET_ID = 2
OUTLET_ID = 3
VOLUME_ID = 109


def read_chunk(f, t):
    """Read a whitespace-delimited chunk of a file, returns chunk, converted to a given type."""
    res = ""
    while True:
        try:
            c = f.read(1)
            if len(c) == 0:
                return None
            if c.isspace() and len(res) > 0:
                return t(res)
            if not c.isspace():
                res += c
        except EOFError:
            return None

class NeutralMapFile(object):
    """NASA's Neutral Map File representation."""
    @staticmethod
    def skip_comments(fp):
        """Skip lines starting with #."""
        while True:
            pos = fp.tell()
            try:
                l = fp.readline()
                if not l.startswith('#'):
                    fp.seek(pos)
                    break
            except IOError:
                break

    def __init__(self, filename=None):
        """Parse boundaries from :filename:."""
        self.__boundaries = []
        if filename is not None:
            fp = open(filename, 'r')
            # Skip initial comments
            NeutralMapFile.skip_comments(fp)
            # Blocks
            l = fp.readline()
            if l.endswith('\\\n'):
                l = l[0:-2]
            nblocks = int(l)
            fp.readline()
            for _ in range(nblocks):
                fp.readline()
            fp.readline()
            # Middle comments
            NeutralMapFile.skip_comments(fp)
            # Boundaries
            for l in fp:
                if l.endswith('\\'):
                    b = l[0:-2].split()
                else:
                    b = l.split()
                if len(b) > 0:
                    if b[0][0] in ['\'', '"']:
                        b[0] = b[0][1:-1]
                    # A ONE_TO_ONE line packs two ranges (12 ints): the
                    # boundary it connects to (B2,F2,...) in addition to
                    # its own (B1,F1,...). GmshFile.consume() uses both
                    # sides to alias node ids across the connection (an
                    # O-grid's full periodic i-closure, or a C-grid's
                    # wake-cut fold) instead of emitting boundary faces.
                    is_one_to_one = b[0].upper() in ['ONE-TO-ONE', 'ONE_TO_ONE']
                    nums_needed = 12 if is_one_to_one else 6
                    b[1:1 + nums_needed] = list(map(int, b[1:1 + nums_needed]))
                    self.__boundaries.append(tuple(b[0:1 + nums_needed]))
            fp.close()

    @property
    def boundaries(self):
        """Get boundaries list."""
        return self.__boundaries

    def __str__(self):
        """Convert boundary to string."""
        return 'Neutral map file / {0:d} boundaries'.format(
            len(self.__boundaries))


class P3DfmtFile(object):
    """P3Dfmt file representation."""
    def __init__(self, filename=None, **kwargs):
        """Construct from components or load from file."""
        if filename:
            self.load(filename=filename)
        else:
            if kwargs is not None:
                self.__nblocks = kwargs['nblocks'] \
                    if 'nblocks' in kwargs else 0
                self.__coords = kwargs['coords'] \
                    if 'coords' in kwargs else None
            else:
                self.__nblocks = None
                self.__coords = None

    @property
    def nblocks(self):
        """Number of blocks in the file."""
        return self.__nblocks

    def idims(self, nblk=1):
        """Return i-dimensions of the file."""
        return self.__coords[nblk - 1][0].shape[0]

    def jdims(self, nblk=1):
        """Return j-dimensions of the file."""
        return self.__coords[nblk - 1][0].shape[1]

    def kdims(self, nblk=1):
        """Return j-dimensions of the file."""
        return self.__coords[nblk - 1][0].shape[2]

    @property
    def coords(self):
        """Return coordinates stored in the file."""
        return self.__coords

    def load(self, filename):
        """Load mesh blocks from the given file."""
        fp = open(filename)

        # Reading number of blocks
        self.__nblocks = read_chunk(fp, int)

        # Reading dimensions
        idims = np.zeros(self.__nblocks, 'i')
        jdims = np.zeros(self.__nblocks, 'i')
        kdims = np.zeros(self.__nblocks, 'i')

        for i in range(self.__nblocks):
            idims[i] = read_chunk(fp, int)
            jdims[i] = read_chunk(fp, int)
            kdims[i] = read_chunk(fp, int)

        # Reading coordinates
        self.__coords = []

        for b in range(self.__nblocks):
            idim = idims[b]
            jdim = jdims[b]
            kdim = kdims[b]
            x = np.zeros((idim, jdim, kdim), 'f8')
            y = np.zeros((idim, jdim, kdim), 'f8')
            z = np.zeros((idim, jdim, kdim), 'f8')

            coords = [x, y, z]
            for c in coords:
                for k in range(kdim):
                    for j in range(jdim):
                        for i in range(idim):
                            c[i, j, k] = read_chunk(fp, float)

            self.__coords.append((x, y, z))

        fp.close()

    def __str__(self):
        """Convert to string."""
        idims = [x.shape[0] for x, _, _ in self.__coords]
        jdims = [x.shape[1] for x, _, _ in self.__coords]
        kdims = [x.shape[2] for x, _, _ in self.__coords]
        return 'P3Dfmt file (blocks: %d/idims: (%s)/jdims: (%s)/kdims: (%s)' % \
            (self.__nblocks, ' '.join(map(str, idims)),
             ' '.join(map(str, jdims)), ' '.join(map(str, kdims)))

    def save(self, filename=None):
        """Save file, to stdout if no filename is given."""
        raise NotImplementedError

    def dump_coords(self):
        """Dump coordinates of the file as a list."""
        for n in range(self.__nblocks):
            idim = self.__coords[n][0].shape[0]
            jdim = self.__coords[n][0].shape[1]
            kdim = self.__coords[n][0].shape[2]

            x, y, z = self.__coords[n]

            for i in range(idim):
                for j in range(jdim):
                    for k in range(kdim):
                        print('%lf %lf %lf' %
                              (x[i, j, k], y[i, j, k], z[i, j, k]))


class GmshFile(object):
    """Gmsh file representation."""

    # For conversion purposes I need only two types of elements:
    # - quadrangle (3) for the faces
    # - hexagon (5) for the cells
    # So there will be no constants.

    __DEFAULT_GEOMETRY_GROUP = VOLUME_ID

    def __init__(self, nodes=None, elements=None, groups=None, filename=None):
        """Construct from components.

        If filename is provided object is loaded from the file.

        :nodes:
            List of node tuples

        :elements:
            List of element tuples

        :groups:
            List of group tuples

        :filename:
            Name of file to load data from
        """
        self.__element_id = 0
        self.__mesh_type = None
        self.__wake_cut_alias = None
        if filename:
            self.__nodes = []
            self.__elements = []
            self.__groups = []
            self.load(filename)
        else:
            self.__nodes = [] if nodes is None else nodes
            self.__elements = [] if elements is None else elements
            self.__groups = [] if groups is None else groups

    @property
    def nodes(self):
        """Return nodes of the current file."""
        return self.__nodes

    @property
    def elements(self):
        """Return elements of the current file."""
        return self.__elements

    @property
    def groups(self):
        """Return physical groups of the current file."""
        return self.__groups

    def load(self, filename=None):
        """Load nodes, elements, and groups from the given file."""
        raise NotImplementedError

    def __str__(self):
        """Create string representation of the file."""
        return 'GMSH file (nodes: %d, elements: %d, groups: %d)' % \
            (len(self.__nodes), len(self.__elements), len(self.__groups))

    def save(self, filename=None):
        """Save file, to stdout if no filename is given."""
        if filename:
            fp = open(filename, 'w')
        else:
            fp = sys.stdout

        GmshFile._write_header(fp)
        self._write_groups(fp)
        self._write_nodes(fp)
        self._write_elements(fp)

    @staticmethod
    def _write_header(out):
        """Write standard Gmsh file header."""
        out.write('$MeshFormat\n')
        out.write('2.2 0 8\n')
        out.write('$EndMeshFormat\n')

    def _write_groups(self, out):
        """Write Gmsh file physical groups."""
        out.write('$PhysicalNames\n')
        out.write('%d\n' % len(self.__groups))
        for grp in self.__groups:
            out.write('%d %d "%s"\n' % grp)
        out.write('$EndPhysicalNames\n')

    def _write_nodes(self, out):
        """Write Gmsh file nodes."""
        out.write('$Nodes\n')
        out.write('%d\n' % len(self.__nodes))
        for node in self.__nodes:
            out.write('%d %15.13e %15.13e %15.13e\n' % node)
        out.write('$EndNodes\n')

    def _write_elements(self, out):
        """Write Gmsh file elements."""
        out.write('$Elements\n')
        out.write('%d\n' % len(self.__elements))
        for el in self.__elements:
            out.write('%s\n' % ' '.join(map(str, el)))
        out.write('$EndElements\n')

    def consume(self, p3dfmt_file, idim, angle_of_attack,mesh_type,mapfile=None):
        """Convert P3Dfmt file into self.

        :p3dfmt_file:
            P3DfmtFile object to convert.

        :idim:
            Number of grid points in the i-direction.

        :angle_of_attack:
            Angle of attack for the simulation, in degrees. Used to
            classify each farfield-type boundary point as inlet or
            outlet (see GmshFile._classify_inlet_outlet).

        :mesh_type:
            "OGRD" or "CGRD". OGRD blocks are periodically closed in i
            (i=idim-1 is physically the same point as i=0); CGRD blocks are
            not (i=0/i=idim-1 are the two distinct wake-cut/outlet corners),
            so no index wraparound must be applied.

        :mapfile:
            Neutral map file name, for boundary faces
        """
        self.__mesh_type = mesh_type
        self.__wake_cut_alias = None

        if p3dfmt_file.nblocks > 1:
            # A genuinely different shape from the single-block OGRD/CGRD
            # path below (untouched by this branch): e.g. a 2-block
            # C-grid where the blunt trailing edge is resolved as its own
            # small block. Block roles (which face is the wall,
            # farfield, ...) are read from each boundary's own NAME
            # rather than assumed from its face id -- the TE-closure
            # block's own wall sits on a different face than the main
            # block's -- and blocks are welded together by matching
            # coordinates (_build_multiblock_weld) instead of the
            # single-block wake-cut's index-offset alias, which only
            # ever remaps within ONE block.
            if mesh_type != 'CGRD':
                raise ValueError(
                    f"Multiblock input ({p3dfmt_file.nblocks} blocks) is only "
                    f"supported for mesh_type='CGRD', got '{mesh_type}'"
                )
            self._consume_multiblock(p3dfmt_file, angle_of_attack, mapfile)
            return

        one_to_one = [b for b in mapfile.boundaries if b[0].upper() in ('ONE_TO_ONE', 'ONE-TO-ONE')]
        if self.__mesh_type == 'CGRD' and one_to_one:
            # The wake-cut fold: build the node-id alias BEFORE consuming
            # blocks/boundaries, since both the interior volume mesh and
            # the other boundaries reference j=0 nodes across the fold.
            self.__wake_cut_alias = self._build_wake_cut_alias(one_to_one[0])

        # Seed the fixed-id groups directly (rather than auto-numbering
        # and renumbering them later at the .geo stage -- see the
        # WALL_ID/INLET_ID/OUTLET_ID/VOLUME_ID comment above). Face1/2
        # (SYMMETRY-Y) and any other auto-created groups pick up ids
        # above these via _next_group_id()'s max()+1.
        self.__groups.append((3, VOLUME_ID, 'VolumeCode'))
        self.__groups.append((2, WALL_ID, 'wall'))
        self.__groups.append((2, INLET_ID, 'inlet'))
        self.__groups.append((2, OUTLET_ID, 'outlet'))
        for blkn in range(p3dfmt_file.nblocks):
            self._consume_block(p3dfmt_file, blkn)

        for bdry in mapfile.boundaries:
            if bdry[0].upper() in ('ONE_TO_ONE', 'ONE-TO-ONE'):
                # Not a physical boundary: an OGRD's full periodic
                # i-closure and a CGRD's wake-cut fold are both realized
                # by aliasing node ids (see _p3d_node_id_closed_i), so no
                # boundary face is emitted for either.
                continue
            self._gen_boundary(p3dfmt_file, idim, angle_of_attack,bdry)

    @staticmethod
    def _build_wake_cut_alias(bdry):
        """Build the node-id alias for a CGRD wake-cut ONE_TO_ONE boundary.

        bdry = (name, b1,f1,s1,e1,s2,e2, b2,f2,s1,e1,s2,e2), already
        remapped onto the extruded 3D block's face numbering by
        write_ext_nmf. Both sides must be face5 (j=0, the wall/wake-cut
        row); the I-range lives in (s2,e2) for a j=const face (see
        mesh_extrusion._remap_2d_face_to_3d). Side B's range is listed in
        reverse order (s2 > e2), which encodes the fold's orientation.
        """
        f1, s2_a, e2_a = bdry[2], bdry[5], bdry[6]
        f2, s2_b, e2_b = bdry[8], bdry[11], bdry[12]
        if f1 != 5 or f2 != 5:
            raise ValueError(
                f"Expected the CGRD wake-cut ONE_TO_ONE to connect two "
                f"sub-ranges of face5 (j=0, the wall row), got faces {f1},{f2}"
            )
        return {
            'i_lo_a': s2_a - 1, 'i_hi_a': e2_a - 1,
            'i_lo_b': min(s2_b, e2_b) - 1, 'i_hi_b': max(s2_b, e2_b) - 1,
            'offset': s2_b - s2_a,  # i_b = offset - i_a (and vice versa)
        }

    # ---- Multiblock path (>1 block; see consume()'s dispatch) ----------
    #
    # The single-block OGRD/CGRD path above welds at most one fold, always
    # within one block, by an index offset (_build_wake_cut_alias). A
    # multiblock case can have several folds, each possibly crossing two
    # DIFFERENT blocks -- there's no single shared index offset that
    # describes that in general, so instead every node keeps its own
    # global id (_p3d_node_id, already block-aware) and welded folds are
    # resolved by matching physical coordinates: for each ONE_TO_ONE
    # boundary, the corresponding node-id pairs from its two sides are
    # unioned into a small union-find, and every consumer below (node
    # listing, hex elements, boundary faces) looks its ids up through
    # that instead of using an id it computed itself directly.

    @staticmethod
    def _face56_node_ids(p3df, blkn, face, s1, e1, s2, e2):
        """All node ids (not quad corners) along a face5 (j=0) / face6
        (j=jmax) sub-range of one block, in i-outer/k-inner order -- used
        to pair up the two sides of a cross-block ONE_TO_ONE weld
        (_build_multiblock_weld). s1/e1 is the K sub-range (always
        ascending -- spanwise, matches on both sides of any weld by
        construction); s2/e2 is the I sub-range, and its own direction
        (s2>e2 for a reversed side) encodes the fold's orientation --
        same convention as _build_wake_cut_alias's single-block case, see
        mesh_extrusion._remap_2d_face_to_3d.
        """
        x, _, _ = p3df.coords[blkn]
        jdim = x.shape[1]
        j = 0 if face == 5 else jdim - 1
        i_vals = range(s2 - 1, e2 - 2, -1) if s2 > e2 else range(s2 - 1, e2)
        k_vals = range(s1 - 1, e1)
        return [GmshFile._p3d_node_id(p3df, blkn, i, j, k) for i in i_vals for k in k_vals]

    def _build_multiblock_weld(self, p3dfmt_file, mapfile):
        """Build the cross-block node-id union from every ONE_TO_ONE
        boundary, and return a find(raw_id) closure resolving any node's
        raw (_p3d_node_id) id to its weld's canonical (smallest) id --
        callers use find() everywhere a raw id would otherwise be used,
        both for the welded nodes themselves and for every other node
        (find() is a no-op for an id nothing was ever unioned with)."""
        parent = {}

        def find(x):
            root = x
            while parent.get(root, root) != root:
                root = parent[root]
            while parent.get(x, x) != root:
                parent[x], x = root, parent.get(x, x)
            return root

        def union(a, b):
            ra, rb = find(a), find(b)
            if ra != rb:
                lo, hi = (ra, rb) if ra < rb else (rb, ra)
                parent[hi] = lo

        for bdry in mapfile.boundaries:
            if bdry[0].upper() not in ('ONE_TO_ONE', 'ONE-TO-ONE'):
                continue
            b1, f1, s1, e1, s2, e2 = bdry[1:7]
            b2, f2, s1b, e1b, s2b, e2b = bdry[7:13]
            if f1 not in (5, 6) or f2 not in (5, 6):
                raise ValueError(
                    f"Multiblock ONE_TO_ONE welding only supports face5/6 "
                    f"(j=const) connections, got face {f1} (block {b1}) "
                    f"matched to face {f2} (block {b2})"
                )
            ids_a = self._face56_node_ids(p3dfmt_file, b1 - 1, f1, s1, e1, s2, e2)
            ids_b = self._face56_node_ids(p3dfmt_file, b2 - 1, f2, s1b, e1b, s2b, e2b)
            if len(ids_a) != len(ids_b):
                raise ValueError(
                    f"ONE_TO_ONE weld side lengths differ: {len(ids_a)} "
                    f"(block {b1}) vs {len(ids_b)} (block {b2})"
                )
            for a, b in zip(ids_a, ids_b):
                union(a, b)

        return find

    def _consume_multiblock(self, p3dfmt_file, angle_of_attack, mapfile):
        """Convert a multiblock (>1 block) P3Dfmt file into self -- the
        dispatch target from consume() when p3dfmt_file.nblocks > 1; see
        this section's own header comment above."""
        find = self._build_multiblock_weld(p3dfmt_file, mapfile)

        self.__groups.append((3, VOLUME_ID, 'VolumeCode'))
        self.__groups.append((2, WALL_ID, 'wall'))
        self.__groups.append((2, INLET_ID, 'inlet'))
        self.__groups.append((2, OUTLET_ID, 'outlet'))

        for blkn in range(p3dfmt_file.nblocks):
            self._consume_block_multiblock(p3dfmt_file, blkn, find)

        for bdry in mapfile.boundaries:
            if bdry[0].upper() in ('ONE_TO_ONE', 'ONE-TO-ONE'):
                continue
            self._gen_boundary_multiblock(p3dfmt_file, angle_of_attack, bdry, find)

    def _consume_block_multiblock(self, p3dfmt_file, blkn, find):
        """Node/element emission for one block of a multiblock mesh: every
        (i,j,k) still gets a node, but only under its weld's canonical id
        (find() is a no-op for a non-welded node), and every hex corner
        is looked up through find() too -- see _build_multiblock_weld."""
        x, y, z = p3dfmt_file.coords[blkn]
        idim, jdim, kdim = x.shape

        for i in range(idim):
            for j in range(jdim):
                for k in range(kdim):
                    raw_id = self._p3d_node_id(p3dfmt_file, blkn, i, j, k)
                    if find(raw_id) != raw_id:
                        continue
                    self.__nodes.append((raw_id, x[i, j, k], y[i, j, k], z[i, j, k]))

        shifts = [
            [-1, -1, -1],
            [-1, 0, -1],
            [-1, 0, 0],
            [-1, -1, 0],
            [0, -1, -1],
            [0, 0, -1],
            [0, 0, 0],
            [0, -1, 0],
        ]
        for i in range(1, idim):
            for j in range(1, jdim):
                for k in range(1, kdim):
                    el_id = self.get_next_element_id()
                    el = [el_id, 5, 2, VOLUME_ID, VOLUME_ID]
                    for s in shifts:
                        raw = self._p3d_node_id(p3dfmt_file, blkn, i + s[0], j + s[1], k + s[2])
                        el.append(find(raw))
                    self.__elements.append(el)

    def _gen_boundary_multiblock(self, p3df, angle_of_attack, bdry, find):
        """Emit one non-ONE_TO_ONE boundary's quad faces. Unlike the
        single-block _gen_boundary, role (wall vs. farfield-outlet vs.
        farfield-arc) is read from the boundary's own NAME rather than
        assumed from its face id, since a second block's wall doesn't
        have to sit on face5 -- e.g. the TE-closure block's own wall is
        its face3. Node ids go through find() (see
        _build_multiblock_weld) instead of the single-block path's
        _p3d_node_id_closed_i."""
        name = bdry[0].upper()
        blkn = bdry[1] - 1
        face = bdry[2]
        s1, e1, s2, e2 = bdry[3:7]
        x, _, _ = p3df.coords[blkn]
        imax = x.shape[0] - 1
        jmax = x.shape[1] - 1
        kmax = x.shape[2] - 1

        def nid(i, j, k):
            return find(self._p3d_node_id(p3df, blkn, i, j, k))

        # Face 1/2 (k=const, spanwise ends): every block's own k=0/kmax
        # faces merge into ONE shared physical group per role (not one
        # group per block, unlike the single-block path's per-record
        # unique group) -- write_ext_nmf emits one SYMMETRY-Y record per
        # block precisely so this loop sees, and merges, all of them.
        if face in (1, 2):
            role = 'k0-SYMMETRY-Y' if face == 1 else 'kmax-SYMMETRY-Y'
            gid = self._get_or_create_group(role)
            k_fixed = 0 if face == 1 else kmax
            for j in range(s2 - 1, e2 - 1):
                for i in range(s1 - 1, e1 - 1):
                    el_id = self.get_next_element_id()
                    n1, n2 = nid(i, j, k_fixed), nid(i + 1, j, k_fixed)
                    n3, n4 = nid(i + 1, j + 1, k_fixed), nid(i, j + 1, k_fixed)
                    conn = [n1, n4, n3, n2] if face == 1 else [n1, n2, n3, n4]
                    self.__elements.append([el_id, 3, 2, gid, gid] + conn)
            return

        # Face 3/4 (i=const): either a block's own flat C-grid cut end
        # (FARFIELD, always outlet -- same convention as the single-block
        # path) or, for the TE-closure block, its wall (VISCOUS, face3).
        if face in (3, 4):
            i_fixed = 0 if face == 3 else imax
            gid = self._get_or_create_group('wall' if name == 'VISCOUS' else 'outlet')
            for k in range(s2 - 1, e2 - 1):
                for j in range(s1 - 1, e1 - 1):
                    el_id = self.get_next_element_id()
                    n1, n2 = nid(i_fixed, j, k), nid(i_fixed, j + 1, k)
                    n3, n4 = nid(i_fixed, j + 1, k + 1), nid(i_fixed, j, k + 1)
                    conn = [n1, n4, n3, n2] if face == 3 else [n1, n2, n3, n4]
                    self.__elements.append([el_id, 3, 2, gid, gid] + conn)
            return

        # Face 5/6 (j=const): the main block's wall (VISCOUS, face5) or
        # its outer farfield arc (FARFIELD, face6 -- classified inlet/
        # outlet per element exactly as the single-block path does).
        if face in (5, 6):
            j_fixed = 0 if face == 5 else jmax
            for i in range(s2 - 1, e2 - 1):
                for k in range(s1 - 1, e1 - 1):
                    el_id = self.get_next_element_id()
                    n1, n2 = nid(i, j_fixed, k), nid(i + 1, j_fixed, k)
                    n3, n4 = nid(i + 1, j_fixed, k + 1), nid(i, j_fixed, k + 1)
                    if name == 'VISCOUS':
                        gid = self._get_or_create_group('wall')
                        conn = [n1, n2, n3, n4]
                    elif name == 'FARFIELD' and face == 6:
                        group_name = self._classify_inlet_outlet(p3df, blkn, i, k, angle_of_attack)
                        gid = self._get_or_create_group(group_name)
                        conn = [n1, n4, n3, n2]
                    else:
                        raise ValueError(f"Unexpected {name} on face{face} (block {blkn + 1})")
                    self.__elements.append([el_id, 3, 2, gid, gid] + conn)
            return

        raise ValueError(f"Unexpected boundary face id {face} for {name} (block {blkn + 1})")

    @staticmethod
    def __find_smallest_cell(p2dfmt_file):
        dx, dy = 0, 0
        for blk in range(p2dfmt_file.nblocks):
            x, y = p2dfmt_file.coords[blk]
            idim, jdim = x.shape
            dx = x[1, 0] - x[0, 0]
            dy = y[0, 1] - y[0, 0]
            for i in range(1, idim):
                for j in range(1, jdim):
                    dx = min(dx, x[i, j] - x[i - 1, j])
                    dy = min(dy, y[i, j] - y[i, j - 1])
        return min(dx, dy)

    @staticmethod
    def _p3d_node_id(p3dfmt_file, n, i, j, k):
        if n >= p3dfmt_file.nblocks:
            raise IndexError('Block number %d is out of range.' % n)

        basen = 1
        if p3dfmt_file.nblocks > 1:
            for idx in range(n):
                x, _, _ = p3dfmt_file.coords[idx]
                di, dj, dk = x.shape
                basen += di * dj * dk

        x, _, _ = p3dfmt_file.coords[n]
        _, dj, dk = x.shape

        return basen + k + dk * j + dk * dj * i

    def get_next_element_id(self):
        """Generate ID of the next element."""
        self.__element_id += 1
        return self.__element_id
    
    def _p3d_node_id_closed_i(self, p3dfmt_file, n, i, j, k):
        x, _, _ = p3dfmt_file.coords[n]
        idim = x.shape[0]

        # An O-grid is periodically closed in i (i=idim-1 coincides with
        # i=0).
        if self.__mesh_type == 'OGRD' and i == idim - 1:
            i = 0
        # A C-grid's wake-cut fold aliases one sub-range of the j=0 (wall)
        # row onto the matching sub-range on the other side of the fold
        # (see _build_wake_cut_alias) -- everywhere else i=0/i=idim-1 are
        # distinct wake-cut/outlet corners and must never be merged.
        elif (self.__mesh_type == 'CGRD' and self.__wake_cut_alias is not None
                and j == 0 and self.__wake_cut_alias['i_lo_b'] <= i <= self.__wake_cut_alias['i_hi_b']):
            i = self.__wake_cut_alias['offset'] - i

        return GmshFile._p3d_node_id(p3dfmt_file, n, i, j, k)

    def _consume_block(self, p3dfmt_file, blkn):
        x, y, z = p3dfmt_file.coords[blkn]
        idim, jdim, kdim = x.shape

        # For a closed O-grid, the last i-plane coincides with i=0 and is
        # deliberately not stored (see _p3d_node_id_closed_i). A C-grid has
        # no such closure, so every i-plane, including the last, is a real
        # distinct set of nodes -- except the wake-cut-fold side of j=0,
        # which is aliased onto the other side and must likewise not be
        # stored separately.
        idim_nodes = idim - 1 if self.__mesh_type == 'OGRD' else idim

        # Filling nodes list
        for i in range(idim_nodes):
            for j in range(jdim):
                for k in range(kdim):
                    if (self.__mesh_type == 'CGRD' and self.__wake_cut_alias is not None
                            and j == 0 and self.__wake_cut_alias['i_lo_b'] <= i <= self.__wake_cut_alias['i_hi_b']):
                        continue
                    node_id = GmshFile._p3d_node_id(p3dfmt_file, blkn, i, j, k)
                    self.__nodes.append(
                        (node_id, x[i, j, k], y[i, j, k], z[i, j, k]))

        # Generating 3D elements
        shifts = [
            [-1, -1, -1],
            [-1, 0, -1],
            [-1, 0, 0],
            [-1, -1, 0],
            [0, -1, -1],
            [0, 0, -1],
            [0, 0, 0],
            [0, -1, 0],
        ]

        for i in range(1, idim):
            for j in range(1, jdim):
                for k in range(1, kdim):
                    el_id = self.get_next_element_id()
                    el = [el_id, 5, 2, VOLUME_ID, GmshFile.__DEFAULT_GEOMETRY_GROUP]
                    for s in shifts:
                        el.append(
                            self._p3d_node_id_closed_i(p3dfmt_file, blkn, i + s[0],
                                                  j + s[1], k + s[2]))
                    self.__elements.append(el)

    def _next_group_id(self):
        return max(self.__groups, key=lambda n: n[1])[1] + 1

    def _get_or_create_group(self, name, dim=2):
        for group in self.__groups:
            if group[2] == name:
                return group[1]

        gid = self._next_group_id()
        self.__groups.append((dim, gid, name))
        return gid

    def _classify_inlet_outlet(self, p3df, blkn, i, k, angle_of_attack):
        """Classify a farfield-type boundary point as 'inlet' or 'outlet'.

        Approximates the local outward-facing direction as the vector
        from the wall (j=0) to the farfield boundary (j=jmax) at the
        same (i,k), then compares it to the free-stream direction implied
        by angle_of_attack -- U_inf = (cos(AoA), sin(AoA)) in the mesh's
        own x,y frame (the mesh geometry itself is built at zero AoA; AoA
        only rotates the flow direction here). A point is 'inlet' where
        the flow enters (dot < 0), 'outlet' where it leaves (dot >= 0).
        """
        x, y, _ = p3df.coords[blkn]
        jmax = x.shape[1] - 1

        dx = x[i, jmax, k] - x[i, 0, k]
        dy = y[i, jmax, k] - y[i, 0, k]
        norm = math.hypot(dx, dy)
        if norm > 0:
            dx, dy = dx / norm, dy / norm

        aoa_rad = math.radians(angle_of_attack)
        ux, uy = math.cos(aoa_rad), math.sin(aoa_rad)

        return 'inlet' if (dx * ux + dy * uy) < 0 else 'outlet'

    def _gen_boundary(self, p3df, idim, angle_of_attack, bdry):
        blkn = bdry[1] - 1
        x, _, _ = p3df.coords[blkn]

        imax = x.shape[0] - 1
        jmax = x.shape[1] - 1
        kmax = x.shape[2] - 1

        s1, e1, s2, e2 = bdry[3:7]

        # Face 1/2 (k=const, spanwise ends): own group per call, unrelated
        # to the wall/inlet/outlet classification below.
        if bdry[2] in (1, 2):
            gid = self._next_group_id()
            self.__groups.append((2, gid, 'b{0:d}-{1}'.format(gid, bdry[0])))
            k_fixed = 0 if bdry[2] == 1 else kmax
            for j in range(s2 - 1, e2 - 1):
                for i in range(s1 - 1, e1 - 1):
                    el_id = self.get_next_element_id()
                    n1 = self._p3d_node_id_closed_i(p3df, blkn, i,     j,     k_fixed)
                    n2 = self._p3d_node_id_closed_i(p3df, blkn, i + 1, j,     k_fixed)
                    n3 = self._p3d_node_id_closed_i(p3df, blkn, i + 1, j + 1, k_fixed)
                    n4 = self._p3d_node_id_closed_i(p3df, blkn, i,     j + 1, k_fixed)
                    conn = [n1, n4, n3, n2] if bdry[2] == 1 else [n1, n2, n3, n4]
                    self.__elements.append([el_id, 3, 2, gid, gid] + conn)
            return

        # Face 5: j = 0 (the airfoil wall).
        if bdry[2] == 5:
            gid = self._get_or_create_group('wall')
            for i in range(s2 - 1, e2 - 1):
                for k in range(s1 - 1, e1 - 1):
                    el_id = self.get_next_element_id()
                    n1 = self._p3d_node_id_closed_i(p3df, blkn, i,     0, k)
                    n2 = self._p3d_node_id_closed_i(p3df, blkn, i + 1, 0, k)
                    n3 = self._p3d_node_id_closed_i(p3df, blkn, i + 1, 0, k + 1)
                    n4 = self._p3d_node_id_closed_i(p3df, blkn, i,     0, k + 1)
                    self.__elements.append([el_id, 3, 2, gid, gid, n1, n2, n3, n4])
            return

        # Face 3/4: i = 0 / i = imax (a C-grid's two wake-end faces,
        # together forming the flat right-hand boundary of the C, as
        # opposed to face6's curved arc). Always outlet, regardless of
        # angle_of_attack.
        if bdry[2] in (3, 4):
            i_fixed = 0 if bdry[2] == 3 else imax
            gid = self._get_or_create_group('outlet')
            for k in range(s2 - 1, e2 - 1):
                for j in range(s1 - 1, e1 - 1):
                    el_id = self.get_next_element_id()
                    n1 = self._p3d_node_id_closed_i(p3df, blkn, i_fixed, j,     k)
                    n2 = self._p3d_node_id_closed_i(p3df, blkn, i_fixed, j + 1, k)
                    n3 = self._p3d_node_id_closed_i(p3df, blkn, i_fixed, j + 1, k + 1)
                    n4 = self._p3d_node_id_closed_i(p3df, blkn, i_fixed, j,     k + 1)
                    conn = [n1, n4, n3, n2] if bdry[2] == 3 else [n1, n2, n3, n4]
                    self.__elements.append([el_id, 3, 2, gid, gid] + conn)
            return

        # Face 6: j = jmax (the outer farfield arc), split element by
        # element into inlet/outlet based on the local flow direction.
        if bdry[2] == 6:
            for i in range(s2 - 1, e2 - 1):
                for k in range(s1 - 1, e1 - 1):
                    el_id = self.get_next_element_id()
                    n1 = self._p3d_node_id_closed_i(p3df, blkn, i,     jmax, k)
                    n2 = self._p3d_node_id_closed_i(p3df, blkn, i + 1, jmax, k)
                    n3 = self._p3d_node_id_closed_i(p3df, blkn, i + 1, jmax, k + 1)
                    n4 = self._p3d_node_id_closed_i(p3df, blkn, i,     jmax, k + 1)
                    group_name = self._classify_inlet_outlet(p3df, blkn, i, k, angle_of_attack)
                    gid = self._get_or_create_group(group_name)
                    self.__elements.append([el_id, 3, 2, gid, gid, n1, n4, n3, n2])
            return

        raise ValueError('Unknown block face identifier.')

def p3d2gmsh(p3d_file,idim, angle_of_attack, mesh_type, map_file=None, output_file=None):
    """
    Convert one Plot3D file to Gmsh .msh file.
    This allows calling p3d2gmsh() directly from run_pipeline.py.

    Returns (output_file, groups), where groups is the list of
    (dim, id, name) physical groups actually generated -- since inlet vs
    outlet is decided per-element from geometry and angle_of_attack, the
    ids can't be predicted ahead of time, so callers (e.g. write_geo_file)
    need this to know which entity ids ended up tagged 'wall'/'inlet'/
    'outlet'.
    """

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

    p3d = P3DfmtFile()
    p3d.load(p3d_file)

    nmf = NeutralMapFile(map_file)

    gmsh = GmshFile()
    gmsh.consume(p3d,idim=idim, angle_of_attack=angle_of_attack,mesh_type=mesh_type,mapfile=nmf)
    gmsh.save(output_file)

    return output_file, gmsh.groups

def main():
    parser = argparse.ArgumentParser(
        description="Convert P3Dfmt mesh into Gmsh mesh",
        add_help=True
    )

    parser.add_argument(
        "files",
        nargs="+",
        help="Plot3D files to convert"
    )

    parser.add_argument(
        "-m",
        "--map-file",
        nargs=1,
        help="Neutral Map File. If omitted, script uses <filename>.nmf"
    )

    parser.add_argument(
        "-o",
        "--output-file",
        nargs=1,
        help="Output Gmsh .msh file. If omitted, script uses <filename>.msh"
    )
    parser.add_argument(
        "--idim",
        type=int,
        required=True,
        help="Number of airfoil surface points"
    )

    parser.add_argument(
        "--aoa",
        type=float,
        required=True,
        help="Angle of attack in degrees"
    )
    parser.add_argument(
        "--mesh-type",
        choices=["OGRD", "CGRD"],
        required=True,
        help="Mesh topology: OGRD (O-grid, periodically closed in i) or CGRD (C-grid, wake-cut)"
    )
    args = parser.parse_args()

    for fn in args.files:
        map_file = args.map_file[0] if args.map_file is not None else None
        output_file = args.output_file[0] if args.output_file is not None else None

        p3d2gmsh(
            p3d_file=fn,
            idim=args.idim,
            angle_of_attack=args.aoa,
            mesh_type=args.mesh_type,
            map_file=map_file,
            output_file=output_file,
        )


if __name__ == '__main__':
    main()
