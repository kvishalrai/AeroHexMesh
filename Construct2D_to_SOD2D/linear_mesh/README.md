# p3d-to-sod2d

Converts a 2D airfoil mesh (O-grid or C-grid) into a partitioned, 3D
spanwise-extruded SOD2D mesh, via Gmsh.

This is the **middle step** of a 3-step pipeline:

1. **Construct2D** (separate, external) — generates the 2D Plot3D mesh
   (`.p3d`) and its Neutral Map File (`.nmf`) for a given airfoil and mesh
   topology (O-grid or C-grid). Not included here; run it yourself.
2. **This repo** — extrudes the 2D mesh spanwise into 3D, remaps the NMF
   onto the extruded block, converts it to Gmsh format with wall/inlet/
   outlet boundary classification, exports to SOD2D's HDF5 format, and
   partitions it for parallel runs.
3. **SOD2D** (separate, external) — provides the `tool_meshConversorPar`
   partitioning binary this pipeline calls, and (separately) a tool to
   smooth the airfoil boundary into a high-order-conforming curved mesh
   for production runs.

## Pipeline stages

```
2D .p3d/.nmf  →  mesh_extrusion.py   →  extruded 3D .p3d + remapped .nmf
              →  p3d_to_gmsh.py      →  .msh (wall/inlet/outlet/periodic tagged)
              →  Gmsh (.geo)         →  order-elevated, periodic .msh
              →  sod2d_tools/gmsh2sod2d.py  →  SOD2D .h5
              →  sod2d_tools/tool_meshConversorPar  →  partitioned .hdf (per rank)
```

`run_pipeline.py` orchestrates all of the above for one airfoil.

## Prerequisites

- Python 3 with `numpy` (and `h5py` for the `gmsh2sod2d.py` step).
- [Gmsh](https://gmsh.info/), invoked as a CLI (`gmsh airfoil_per.geo -0`).
- MPI (`mpirun`) and HDF5, for the partitioning step.
- `sod2d_tools/tool_meshConversorPar` — a compiled SOD2D binary,
  **not included in this repo** (platform-specific). Build it from the
  SOD2D repo, or copy an existing build, into `sod2d_tools/`.

## Usage

```
python3 run_pipeline.py \
    --work-dir runs/env_0_cgrd \
    --config flow_config_cgrd.json \
    --airfoil-file airfoil
```

`--work-dir` must already contain the 2D mesh pair produced by Construct2D,
named `{airfoil-file}.p3d` / `{airfoil-file}.nmf` (see `examples/` for
sample O-grid and C-grid pairs). `--config` is required and has no
default — pick whichever `flow_config_*.json` matches the mesh topology
you're running (`mesh_type` must agree with the actual `.nmf` content).

### flow_config*.json fields used by this pipeline

| Field | Meaning |
|---|---|
| `mesh_type` | `"OGRD"` (periodically closed in i) or `"CGRD"` (wake-cut, not closed) — must match the actual 2D `.nmf` topology |
| `idim`, `jmax` | 2D mesh dimensions, must match the `.nmf`'s declared `IDIM`/`JDIM` |
| `angle_of_attack` | degrees; drives the inlet/outlet split on the farfield arc (see below) |
| `z_spanwise_len`, `z_spanwise_planes` | spanwise extrusion length and number of planes (>= 2) |
| `porder` | Gmsh element order for the final mesh |
| `num_partitions` | number of MPI ranks to partition the mesh for |

(The config files also carry additional fields consumed by the downstream
SOD2D solver run itself, not by this pipeline.)

## Physical-id convention

Every mesh produced by this pipeline (O-grid or C-grid) uses the same
fixed physical-group ids, defined once in `p3d_to_gmsh.py`
(`WALL_ID`/`INLET_ID`/`OUTLET_ID`/`VOLUME_ID`) and `mesh_extrusion.py`
(`PERIODIC_ID`):

| id | name |
|---|---|
| 1 | WALL |
| 2 | INLET |
| 3 | OUTLET |
| 4 | Periodic (spanwise) |
| 109 | VolumeCode |

Wall/inlet/outlet are classified per boundary element from geometry: at
each farfield-type point, the local outward direction (wall→farfield
vector) is compared against the free-stream direction implied by
`angle_of_attack` (`U_inf = (cos(AoA), sin(AoA))` in the mesh's own,
zero-AoA frame) — inlet where flow enters, outlet where it leaves. For a
C-grid, the two wake-end faces (the flat "right-hand" boundary of the C,
as opposed to the curved farfield arc) are always OUTLET regardless of
AoA.

**Caveat:** the AoA sign convention above (flow tilts toward +y) was
chosen as a reasonable default and works consistently across O-grid and
C-grid, but wasn't checked against an independently-known-good case.
Verify it before trusting results at nonzero AoA, especially if flipping
it turns out to matter for your solver setup.

## Mesh topology notes

- **O-grid**: the periodic i-closure (i=0 / i=idim-1 are the same physical
  point) is realized purely by node-id aliasing in `p3d_to_gmsh.py` — no
  boundary face or Gmsh-side periodic surface is needed for it.
- **C-grid**: the wake-cut fold (the two i-sub-ranges the `.nmf` marks
  `ONE_TO_ONE` on the wall face) is merged the same way — by node-id
  aliasing, not a Gmsh `Periodic Surface`. An earlier attempt to expose it
  as two separate boundary surfaces and relate them with a Gmsh periodic
  directive turned out to be unreliable (see git history if curious); the
  node-alias approach mirrors the O-grid closure and has been verified
  end-to-end through Gmsh, `gmsh2sod2d.py`, and `tool_meshConversorPar`.

## Directory layout

```
run_pipeline.py       # orchestrates the full pipeline for one airfoil
mesh_extrusion.py      # 2D->3D extrusion, NMF remapping, .geo file writer
p3d_to_gmsh.py         # Plot3D+NMF -> Gmsh .msh, wall/inlet/outlet classification
sod2d_tools/
  gmsh2sod2d.py        # vendored SOD2D tool: Gmsh .msh -> SOD2D .h5
  tool_meshConversorPar # vendored SOD2D binary (not included, see Prerequisites)
examples/
  naca0012.p3d/.nmf         # O-grid example
  naca0012_sharp.p3d/.nmf   # C-grid example
flow_config_ogrd.json
flow_config_cgrd.json
```
