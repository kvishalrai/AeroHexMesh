# p3d-to-sod2d

This stage takes the flat 2D grid Construct2D produced and turns it into a
partitioned 3D mesh in SOD2D's own format, ready for
`../high_order_mesh/` to smooth its wall boundary. If you haven't read it
yet, the [top-level README's "big picture" section](../README.md#the-big-picture)
explains *why* each of these steps exists before you dive into the *how*
below.

## What's in a `.p3d`/`.nmf` pair?

Construct2D writes two files per mesh, and this whole pipeline starts
from them:

- **`.p3d`** (Plot3D format) — just a big grid of `(x, y)` coordinates,
  arranged in rows and columns like a spreadsheet. One "direction" of the
  grid runs around the airfoil (circumferentially); the other runs away
  from it, from the wall out to the far boundary.
- **`.nmf`** (Neutral Map File) — a small text file that says what each
  edge of that grid of points *means*: which edge is the solid wall,
  which is the inflow, which is the outflow, and (for a C-grid) which two
  strips of points are actually the same physical wake line, folded onto
  itself.

## The algorithm

This whole pipeline (this directory *and* `../high_order_mesh/`) is four
steps, each one taking the previous step's output as input. This
directory does the first three; `../high_order_mesh/` does the fourth:

```
INPUT: airfoil.p3d, airfoil.nmf           (from Construct2D)

STEP 1  Stretch the 2D grid into a 3D grid by repeating it along a
        straight line (the wing's span), and remap the .nmf's
        boundary info onto that extruded 3D grid.
        →  an extruded 3D Plot3D grid + a remapped .nmf

STEP 2  Convert to a mesh format Gmsh/SOD2D understand, tagging
        which boundary is the wall/inlet/outlet/spanwise-periodic,
        and add the extra nodes a high-order mesh needs (still
        placed by straight-line interpolation -- the wall is still
        jagged at this point).
        →  an order-elevated Gmsh mesh (.msh)

STEP 3  Convert to SOD2D's own mesh format, then split ("partition")
        the mesh into one piece per MPI rank you plan to run on.
        →  a partitioned SOD2D mesh (.hdf, one file set per rank)

STEP 4  (in ../high_order_mesh/ -- see its own README) Fit a smooth
        curve through the wall's corner points, snap every wall
        mesh node onto it, and elastically relax the rest of the
        mesh to match.
```

`run_pipeline.py` runs steps 1–3 above for you, in one command — unlike
the NekRS pipeline, which runs each step as a separate command (see
[`Construct2D_to_NEKRS/linear_mesh/README.md`](../../Construct2D_to_NEKRS/linear_mesh/README.md)
if you're curious why: it comes down to whether the underlying tools are
callable as a Python library, which SOD2D's are and Nek5000's aren't).

## Prerequisites

- Python 3 with `numpy` (and `h5py` for the `gmsh2sod2d.py` step).
- [Gmsh](https://gmsh.info/), invoked as a CLI (`gmsh airfoil_per.geo -0`).
- MPI (`mpirun`) and HDF5, for the partitioning step.
- `sod2d_tools/gmsh2sod2d.py` and `sod2d_tools/tool_meshConversorPar` both
  come from the SOD2D repo itself (`../CFD_code/sod2d_gitlab/utils/gmsh2sod2d/`
  and `.../tool_meshConversorPar/`), **not included in this repo** as
  build/binary artifacts:
  - `gmsh2sod2d.py` is a plain script — copy it (or symlink it) in from the
    submodule.
  - `tool_meshConversorPar` must be **compiled** — it's a CPU-only tool
    (`-DTOOL_MESHPART=ON` when configuring CMake, e.g.
    `utils/buildCPU.sh <threads> <isMN> <setTPP> 1`; not built by the GPU
    build). Copy the resulting binary into `sod2d_tools/`.

### Cluster-specific module/environment setup

`run_pipeline.py` shells out to Gmsh, `python3`, and MPI for its Gmsh,
`gmsh2sod2d.py`, and partitioning steps, and needs each to be on `PATH`
when it runs. The defaults (`AEROHEXMESH_MODULE_SETUP` /
`AEROHEXMESH_PYTHON_MODULE_SETUP`, both env vars) are BSC MareNostrum 5
specific (`module getdefault sod2d` is an MN5-only alias) — **override
them on any other system** rather than editing `run_pipeline.py`:

```bash
export AEROHEXMESH_MODULE_SETUP="module load gmsh openmpi hdf5 python3"
export AEROHEXMESH_PYTHON_MODULE_SETUP=""   # only needed if your base
                                             # module set's python3 lacks
                                             # numpy/h5py
python3 run_pipeline.py --work-dir ... --config ... --airfoil-file ...
```

Leave both unset to keep the MN5 defaults, or set either to `""` for a
no-op (e.g. if you've already activated everything yourself, such as via
a virtualenv, before running the pipeline).

## Usage

```
python3 run_pipeline.py \
    --work-dir runs/env_0_cgrd \
    --config flow_config_cgrd.json \
    --airfoil-file airfoil
```

`--work-dir` must already contain the 2D mesh pair produced by Construct2D,
named `{airfoil-file}.p3d` / `{airfoil-file}.nmf` (see `examples/` for
sample O-grid and C-grid pairs, and `generate_ogrd_mesh.sh`/
`generate_cgrd_mesh.sh` for ready-to-run examples using them).
`--config` is required and has no default — pick whichever
`flow_config_*.json` matches the mesh topology you're running
(`mesh_type` must agree with the actual `.nmf` content).

### flow_config*.json fields

Every field is read by this pipeline — there's nothing extra in these files
for a downstream solver run to consume. The 2D mesh dimensions (`IDIM`/`JDIM`)
are read directly from the `.nmf` file's own header, not from config.

| Field | Meaning |
|---|---|
| `mesh_type` | `"OGRD"` (periodically closed in i) or `"CGRD"` (wake-cut, not closed) — must match the actual 2D `.nmf` topology |
| `angle_of_attack` | degrees; drives the inlet/outlet split on the farfield arc (see below) |
| `z_spanwise_len`, `z_spanwise_planes` | spanwise extrusion length and number of planes (>= 2) |
| `porder` | Gmsh element order for the final mesh |
| `num_partitions` | number of MPI ranks to partition the mesh for |

## How the wall gets tagged (physical-id convention)

Every boundary element in the mesh gets labeled with one of these fixed
numeric ids, used consistently across the pipeline (`p3d_to_gmsh.py`'s
`WALL_ID`/`INLET_ID`/`OUTLET_ID`/`VOLUME_ID`, `mesh_extrusion.py`'s
`PERIODIC_ID`):

| id | name |
|---|---|
| 1 | WALL |
| 2 | INLET |
| 3 | OUTLET |
| 4 | Periodic (spanwise) |
| 109 | VolumeCode |

Inlet vs. outlet isn't hardcoded — it's worked out geometrically: at each
point on the outer/far boundary, the local outward direction (the
wall→farfield vector) is compared against the free-stream direction
implied by `angle_of_attack` (`U_inf = (cos(AoA), sin(AoA))` in the
mesh's own, zero-AoA frame) — inlet where flow enters, outlet where it
leaves. For a C-grid, the two wake-end faces (the flat "open" ends of the
C, as opposed to the curved farfield arc) are always OUTLET regardless of
AoA.

**Caveat:** the AoA sign convention above (flow tilts toward +y) was
chosen as a reasonable default and works consistently across O-grid and
C-grid, but wasn't checked against an independently-known-good case.
Verify it before trusting results at nonzero AoA, especially if flipping
it turns out to matter for your solver setup.

## O-grid vs. C-grid, in this code

Both grid shapes are handled by the *same* scripts; the difference is
entirely in the shape of the input `.nmf` file. See the
[top-level README](../README.md#the-big-picture) for the plain-language
picture; here's the implementation detail:

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
run_pipeline.py       # STEPS 1-3: orchestrates the full pipeline for one airfoil
mesh_extrusion.py      # step 1: 2D->3D extrusion, NMF remapping, .geo file writer
p3d_to_gmsh.py         # step 2: Plot3D+NMF -> Gmsh .msh, wall/inlet/outlet classification
sod2d_tools/
  gmsh2sod2d.py        # step 3: vendored SOD2D tool, Gmsh .msh -> SOD2D .h5
  tool_meshConversorPar # step 3: vendored SOD2D binary (not included, see Prerequisites)
examples/
  naca0012.p3d/.nmf         # O-grid example
  naca0012_sharp.p3d/.nmf   # C-grid example
generate_ogrd_mesh.sh  # runs the pipeline end-to-end on the O-grid example
generate_cgrd_mesh.sh  # runs the pipeline end-to-end on the C-grid example
flow_config_ogrd.json
flow_config_cgrd.json
```
