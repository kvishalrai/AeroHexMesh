# Construct2D → NekRS

This is the **middle and final stage** of the pipeline described in the
[top-level README](../README.md): it takes the 2D `.p3d`/`.nmf` grid that
Construct2D produces and turns it into a spanwise-periodic, wall-spline-
smoothed 3D mesh that NekRS can run on, via Nek5000's own mesh tools.

```
2D .p3d/.nmf  ──►  linear_mesh/      ──►  high_order_mesh/   ──►  run
(Construct2D)      convert, smooth        apply the 2D            NekRS
                    the wall, extrude 3D   correction to 3D
```

Currently set up and verified for a **C-grid** mesh; the tooling also
supports O-grid (see `linear_mesh/README.md`'s
[Computing case parameters](linear_mesh/README.md#computing-case-parameters)),
it just isn't the case checked in here.

## Stages

1. **[`linear_mesh/`](linear_mesh/README.md)** — converts the native 2D
   Construct2D grid into a linear Gmsh mesh with wall/inlet/outlet
   boundary classification, fits and applies a wall-spline correction on
   the 2D mesh (Paul Fischer's `smooth_geom0`, adapted for this repo — see
   the top-level README's [Credits](../README.md#credits)), and extrudes
   the mesh spanwise into a periodic-in-Z 3D mesh.
2. **[`high_order_mesh/`](high_order_mesh/README.md)** — a NekRS case
   (`example_nekrs/`) whose `usrdat2`/`get_smooth_data` applies the 2D
   wall correction from step 1 onto the real 3D mesh (repeating across
   every spanwise level), then runs the actual flow solve.

## `CFD_code/`

`CFD_code/Nek5000` and `CFD_code/nekRS` are included as git submodules
(see the top-level README's [Credits](../README.md#credits)). Both stages
above depend on tools/binaries built from them:

| Tool | Used by | Build notes |
|---|---|---|
| `Nek5000/bin/gmsh2nek` | `linear_mesh/` (2D Gmsh → `.re2`) | `tools/maketools gmsh2nek` |
| `Nek5000/bin/genmap` | `linear_mesh/smooth_2D/` | `tools/maketools genmap` |
| `Nek5000/bin/{re2torea,reatore2,n2to3}` | `linear_mesh/arglist.sh` (3D extrusion) | `tools/maketools re2torea reatore2 n2to3` (`re2torea`/`reatore2` share one makefile) |
| NekRS itself | `high_order_mesh/example_nekrs/` | large OCCA/GPU build — see nekRS's own [build docs](CFD_code/nekRS/README.md), or use a site-provided install as `example_nekrs/airfoil0.sh` does |

## Prerequisites

- Python 3 with `numpy`.
- [Gmsh](https://gmsh.info/), invoked as a CLI.
- MPI, for Nek5000/NekRS.
- Nek5000 tools and a NekRS install, built as above — see
  [`CFD_code/`](#cfd_code).

Generated mesh, geometry-correction, and run output files (`*.re2`,
`geom_3x3*.dat`, `*.ma2`, compiled binaries, `.cache/`, logs, restart
fields) are gitignored — regenerate them by running the pipeline (see
each stage's own README) rather than expecting them to be present after a
clone.

## Running it end to end

```bash
cd linear_mesh
# ... see linear_mesh/README.md for the full 2D-convert/smooth/extrude
# sequence, ending with arglist.sh writing the 3D mesh into
# ../high_order_mesh/example_nekrs/ and patching its case parameters.

cd ../high_order_mesh/example_nekrs
# ... see high_order_mesh/README.md for the spline-correction copy step
# and running the case.
sbatch airfoil0.sh
```
