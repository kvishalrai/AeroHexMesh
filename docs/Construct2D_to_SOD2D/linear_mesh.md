# `Construct2D_to_SOD2D/linear_mesh/` — math and implementation

Deep-dive companion to
[`Construct2D_to_SOD2D/linear_mesh/README.md`](../../Construct2D_to_SOD2D/linear_mesh/README.md).
Covers the extrusion → NMF remap → Gmsh conversion → SOD2D export
pipeline (`run_pipeline.py`'s four steps).

## 1. Inputs: Plot3D grid + Neutral Map File

Construct2D (run externally) produces a 2D structured grid `x2d, y2d`
(shape `(idim, jdim)`, `j=0` the wall row, `j=jmax` the far boundary) as
ASCII Plot3D, plus a `.nmf` (Neutral Map File) describing what each grid
edge *means*. `_read_2d_nmf()` (`mesh_extrusion.py`) parses it: each
non-`ONE_TO_ONE` line is `name b1 f1 s1 e1 s2 e2` (a named boundary —
`VISCOUS` for the wall, `SYMMETRY-Y`/`XSPLINE` for the farfield, etc.,
over index range `[s1,e1]×[s2,e2]` on face `f1`); a `ONE_TO_ONE` line
packs *two* such ranges (12 ints) — the pair of index ranges that are
physically the same points.

## 2. Spanwise extrusion

`extrude_spanwise()` is the simplest possible 3D construction: pure
translation, no offset/curvature at all —

```
X(i,j,k) = x2d(i,j),   Y(i,j,k) = y2d(i,j),   Z(i,j,k) = z_planes[k]
```

with `z_planes = linspace(0, z_len, z_plane)`. This is the key structural
difference from `etaGrid_to_SOD2D/` (see
[its own math doc](../etaGrid_to_SOD2D/README.md#3-body-wrap-sweep)),
which instead sweeps the cross-section *around* a curved airfoil path —
here the 2D grid (already the full airfoil-to-farfield extent) is simply
copied along a straight spanwise line.

## 3. NMF remapping onto the extruded block

`write_ext_nmf()` derives the 3D NMF from the 2D one, face-by-face:

```
2D face 1/2 (i=const)  →  3D face 3/4 (i=const), same (S1,E1)=J-range,
                          (S2,E2) becomes the FULL spanwise K-range
2D face 3/4 (j=const)  →  3D face 5/6 (j=const), (S1,E1)=K-range (full
                          spanwise extent), (S2,E2)=I-range (the 2D
                          file's own I sub-range)
```

(`_remap_2d_face_to_3d()`). 3D faces 1/2 (`k=const`) are the new
spanwise end caps and have no 2D counterpart — added separately as fixed
`SYMMETRY-Y` entries (these become the Gmsh Periodic pair in step 4). A
`ONE_TO_ONE` line is remapped on **both** sides at once and kept as a
single line — it never becomes a boundary-face physical group (see §3.1).

### 3.1 O-grid closure vs. C-grid wake-cut — the same aliasing mechanism

- **O-grid**: the NMF's `ONE_TO_ONE` connects face `i=0` to face
  `i=imax` (periodic closure all the way around). `GmshFile` never
  stores the `i=idim-1` node plane at all — `_p3d_node_id_closed_i()`
  rewrites any reference to `i=idim-1` as `i=0` before computing the flat
  node id, so the closure is realized purely by **node-id aliasing**, no
  Gmsh-side periodic surface or duplicate node.
- **C-grid**: the NMF's `ONE_TO_ONE` instead connects two *sub-ranges of
  the same face* (`f1 == f2 == 5`, the wall row `j=0`) — the wake-cut
  fold. `_build_wake_cut_alias()` records
  `i_lo_a/i_hi_a` (side A's range), `i_lo_b/i_hi_b` (side B's range,
  listed in *reverse* order in the NMF, encoding the fold's orientation),
  and `offset = s2_b − s2_a`. `_p3d_node_id_closed_i()` then remaps any
  `i` on side B's range to `offset − i` on side A — mirroring the fold —
  and `_consume_block()` skips storing side B's nodes entirely (`continue`
  in the node-generation loop). Mechanically identical to the O-grid case
  (rewrite-before-lookup + skip-on-store), just with a different index
  map.

## 4. Element/node numbering and hex generation

Flat (1-based) node id for structured index `(i,j,k)` in block `n`:

```
id(n,i,j,k) = base_n + k + dk·j + dk·dj·i
```

(`_p3d_node_id()`, `dj,dk` = that block's own J/K extent, `base_n` =
running total of all earlier blocks' node counts — Fortran-style,
K fastest-varying). Hexahedra are built by walking `i∈[1,idim), j∈[1,jdim),
k∈[1,kdim)` and, for each cell, gathering its 8 corners via a fixed
`shifts` table of `(−1,0)` offsets in each axis — the standard
"low/high corner in each of 3 dimensions" enumeration, in Gmsh's own
8-node hex corner order. Every hex gets `VOLUME_ID` as its physical tag.

## 5. Inlet/outlet classification (`_classify_inlet_outlet`)

At each farfield boundary point `(i, k)`:

```
d = normalize( (x[i,jmax,k], y[i,jmax,k]) − (x[i,0,k], y[i,0,k]) )   # wall → farfield, local outward direction
U_inf = (cos(AoA), sin(AoA))                                          # free-stream direction, mesh built at zero AoA
classify = 'inlet' if dot(d, U_inf) < 0 else 'outlet'
```

i.e. inlet where the local outward direction opposes the free stream
(flow entering), outlet where it agrees (flow leaving). **Caveat carried
from the README**: the sign convention (flow tilts toward `+y` for
positive AoA) is a reasonable default, verified to behave consistently
across O-grid and C-grid, but not checked against an independently-known
result — worth confirming before trusting nonzero-AoA results. For a
C-grid, the two wake-end faces (the open ends of the "C") are always
OUTLET regardless of AoA, handled as a separate case in `_gen_boundary()`
rather than through this dot-product test (their own outward direction
isn't meaningfully "wall → farfield").

## 6. Wall-spline fit (`wall_spline.py`)

An **open, natural cubic spline** through the wall's own linear corner
points — used later by SOD2D's `MeshElasticitySolver` to correct the
straight-line-faceted high-order wall nodes (see
[`high_order_mesh.md`](high_order_mesh.md)). Implemented from scratch (no
scipy dependency, so it can be reused directly by SOD2D-adjacent tooling
without adding a dependency) via the standard tridiagonal system for a
natural cubic spline's second derivatives.

Given corner points `(x_i, y_i)` at chord-length arc-lengths `s_i`
(`chord_length_param()`, same cumulative-distance construction as
`etaGrid_to_SOD2D/sweep_mesh.py`'s `arc_length_param()`), fit `y(s)`
piecewise-cubically: on segment `i`, `s_i ≤ t ≤ s_{i+1}`,

```
y(t) = a_i + b_i·dt + c_i·dt² + d_i·dt³,     dt = t − s_i
```

The interior second-derivative values `M_i` solve the classic
tridiagonal system (natural boundary conditions `M_0 = M_n = 0`):

```
h_i = s_{i+1} − s_i
A_i = h_{i−1},   B_i = 2(h_{i−1} + h_i),   C_i = h_i
D_i = 6·[ (y_{i+1} − y_i)/h_i − (y_i − y_{i−1})/h_{i−1} ]
```

solved via the **Thomas algorithm** (forward elimination + back
substitution — `natural_cubic_spline_coeffs()`), then the per-segment
coefficients follow directly:

```
a_i = y_i
b_i = (y_{i+1} − y_i)/h_i − h_i·(2M_i + M_{i+1})/6
c_i = M_i / 2
d_i = (M_{i+1} − M_i) / (6·h_i)
```

Fit independently for `x(s)` and `y(s)`. `check_corner_separation()`
guards against a pathological near-duplicate corner (any two
*non-adjacent* corners closer than `0.1×` the smallest adjacent-segment
length raises an error) before writing the table — corners must be
mutually well-separated for the Fortran side's own nearest-neighbor
corner-matching to be unambiguous (see
[`high_order_mesh.md`](high_order_mesh.md#1-fortran-side-placement)).

For an O-grid, the fit runs through *all* `idim` raw points, even though
`i=0`/`i=idim-1` are the same physical (aliased) node — the two are not
*adjacent* in the fitting sequence, so the final segment
(`idim-2 → idim-1`) is a real, normal-length wall element (the one
closing the loop at the TE), not a degenerate zero-length one. For a
C-grid, the VISCOUS range is sliced as `x2d[s1-1:e1, 0]` — Python's
exclusive slice stop naturally reaching the true inclusive last point,
since NMF ranges are 1-based inclusive on both ends.

## 7. Gmsh tagging, periodicity, order elevation

`write_geo_file()` looks up the physical-group ids `p3d_to_gmsh.py`
actually assigned at runtime (wall/inlet/outlet groups are created
per-element from geometry, so their ids can't be predicted ahead of
time — see `_get_or_create_group()`/`_next_group_id()`), re-tags them to
the fixed convention `WALL=1, INLET=2, OUTLET=3, Periodic=4,
VolumeCode=109`, and emits `Periodic Surface {kmax}={k0} Translate
{0,0,span_z}`. Face1 (`k=0`) is always consumed before face2 (`k=kmax`)
in `GmshFile.consume()`, so the lower id is always the periodic master —
matching the translate direction. `Mesh.ElementOrder = porder; Mesh 3;`
then order-elevates (straight-line placement of new nodes — the defect
[`high_order_mesh.md`](high_order_mesh.md) exists to fix).

## 8. SOD2D export and partitioning

`gmsh2sod2d.py` (vendored from SOD2D's own `utils/gmsh2sod2d/`) reads the
tagged, order-elevated `.msh`, separates interior/boundary/periodic
elements by their physical id, and — for periodic faces — parses the
native Gmsh `$Periodic` section (written by the `Periodic Surface ...`
directive above) directly as text to build the `periodicLinks` node-pair
dataset SOD2D's own HDF5 format needs. `tool_meshConversorPar` then
partitions the resulting HDF5 mesh across MPI ranks
(`write_partition_input_json()` writes its `input.json` driver file).
