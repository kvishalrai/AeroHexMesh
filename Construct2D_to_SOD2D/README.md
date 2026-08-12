# Construct2D → SOD2D

This is the **middle and final stage** of the pipeline described in the
[top-level README](../README.md): it takes the 2D `.p3d`/`.nmf` grid that
Construct2D produces and turns it into a smoothed, partitioned, high-order
3D mesh that SOD2D can run production simulations on.

```
2D .p3d/.nmf  ──►  linear_mesh/      ──►  high_order_mesh/   ──►  copy the smoothed
(Construct2D)      extrude, Gmsh,         wall-boundary            mesh out and run
                    SOD2D export,          smoothing                 production
                    partition               (MeshElasticitySolver)
```

## Stages

1. **[`linear_mesh/`](linear_mesh/README.md)** — extrudes the 2D mesh
   spanwise into 3D, remaps the Neutral Map File onto the extruded block,
   converts it to Gmsh format with wall/inlet/outlet boundary classification,
   exports to SOD2D's HDF5 format, and partitions it for parallel runs.
   `run_pipeline.py` orchestrates all of this for one airfoil (or other
   closed-curve geometry) end to end.
2. **[`high_order_mesh/`](high_order_mesh/README.md)** — the mesh out of
   `linear_mesh/` places its extra high-order boundary nodes by straight-line
   interpolation, so the wall is faceted at high order, not truly curved.
   This step uses SOD2D's `MeshElasticitySolver` to snap every wall node
   onto a cubic spline fit through the wall's own corner points and
   elastically relax the interior mesh to match, writing out the smoothed
   mesh. There's no separate production-run stage in this repo — copy that
   smoothed mesh wherever you want to run it.

## `CFD_code/`

`CFD_code/sod2d_gitlab/` is [SOD2D](https://gitlab.com/bsc_sod2d/sod2d_gitlab)
itself, included as a git submodule (see the top-level README's
[Credits](../README.md#credits)). It's currently pinned to the
`277-witness-points-using-wrong-connectivity` branch, which carries a fix for
high-order wall boundary smoothing (parametric arc-length placement,
replacing an earlier nearest-point-search approach that could collide near
regions of high curvature). Both stages above depend on binaries built
from it:

| Binary | Used by | Build notes |
|---|---|---|
| `sod2d_gitlab/utils/gmsh2sod2d/gmsh2sod2d.py` | `linear_mesh/` (SOD2D export step) | plain script, copy into `linear_mesh/sod2d_tools/` |
| `sod2d_gitlab/tool_meshConversorPar` | `linear_mesh/` (partitioning step) | CPU-only, `-DTOOL_MESHPART=ON` at CMake configure time |
| `sod2d_gitlab`'s `sod2d` app (`MeshElasticitySolver` case) | `high_order_mesh/` | GPU build (`build_gpu`), see `high_order_mesh/airfoil0.sh` |

See `sod2d_gitlab`'s own [README](CFD_code/sod2d_gitlab/README.md) for full
build instructions.

## Prerequisites

- Python 3 with `numpy` (and `h5py` for the SOD2D export step).
- [Gmsh](https://gmsh.info/), invoked as a CLI.
- MPI and HDF5, for mesh partitioning and running SOD2D.
- SOD2D itself, built from the `sod2d_gitlab` submodule — see [`CFD_code/`](#cfd_code)
  above and its own [README](CFD_code/sod2d_gitlab/README.md) for build
  instructions.

Generated mesh, results, and log files (`*.hdf`, `*.h5`, `*.msh`, `*.log`,
etc.) are gitignored — regenerate them by running the pipeline rather than
expecting them to be present after a clone.

## Running it end to end

```bash
# 1. Generate the 2D grid (external, not part of this repo)
#    -> {airfoil}.p3d, {airfoil}.nmf

# 2. linear_mesh/: extrude, convert, partition
cd linear_mesh
python3 run_pipeline.py \
    --work-dir runs/env_0 \
    --config flow_config_ogrd.json \
    --airfoil-file airfoil

# 3. high_order_mesh/: copy the mesh + wall spline table from step 2,
#    then smooth the wall boundary
cp runs/env_0/*.hdf runs/env_0/*_wall_spline.dat ../high_order_mesh/
cd ../high_order_mesh
sbatch airfoil0.sh

# 4. Copy the resulting smoothed mesh (mesh_h5_file_newname in
#    MeshElasticitySolver.json) wherever you want to run production.
```
