# AeroHexMesh

A pipeline for generating high-order, spanwise-extruded 3D meshes (O-grid or
C-grid) for [SOD2D](https://gitlab.com/bsc_sod2d/sod2d_gitlab), a
spectral-element CFD solver, starting from a 2D
[Construct2D](https://sourceforge.net/projects/construct2d/) grid. Although
built and tested around airfoils, nothing in the pipeline is airfoil-specific
— it works from any closed 2D curve Construct2D can mesh.

```
Construct2D            linear_mesh/                    high_order_mesh/         prod_run/
2D closed    ───────►  extrude, remap NMF,   ───────►   snap wall nodes    ────► partitioned
curve .p3d/.nmf         Gmsh .msh, partition             onto fitted curve        production run
                        (SOD2D .hdf, tool_               (MeshElasticitySolver)
                        meshConversorPar)
```

## Repository layout

| Path | What it is |
|---|---|
| `Construct2D/` | Vendored 2D grid generator, designed for airfoils but usable for other closed-curve geometries too (source + Windows binary). See [Credits](#credits). |
| `Construct2D_to_SOD2D/linear_mesh/` | 2D → 3D extrusion, NMF remap, Plot3D → Gmsh conversion, SOD2D export, partitioning. See its own [README](Construct2D_to_SOD2D/linear_mesh/README.md). |
| `Construct2D_to_SOD2D/high_order_mesh/` | Snaps the faceted high-order wall boundary onto a cubic spline fit through the mesh's own wall corner points via SOD2D's `MeshElasticitySolver`, elastically relaxing the interior. See its own [README](Construct2D_to_SOD2D/high_order_mesh/README.md). |
| `Construct2D_to_SOD2D/prod_run/` | Production SOD2D run configs for the smoothed meshes. |
| `Construct2D_to_SOD2D/CFD_code/sod2d_gitlab/` | SOD2D itself, as a git submodule (see [Credits](#credits)). |

## Getting the code

```bash
git clone --recurse-submodules https://github.com/kvishalrai/AeroHexMesh.git
```

If you already cloned without `--recurse-submodules`:

```bash
git submodule update --init --recursive
```

The `sod2d_gitlab` submodule is currently pinned to the
`277-witness-points-using-wrong-connectivity` branch, which carries a fix for
high-order wall boundary smoothing (parametric arc-length placement,
replacing an earlier nearest-point-search approach that could collide near
regions of high curvature).

## Prerequisites

- Python 3 with `numpy` (and `h5py` for the SOD2D export step).
- [Gmsh](https://gmsh.info/), invoked as a CLI.
- MPI and HDF5, for mesh partitioning and running SOD2D.
- SOD2D itself, built from the `sod2d_gitlab` submodule (see its own
  [README](Construct2D_to_SOD2D/CFD_code/sod2d_gitlab/README.md) for build
  instructions) — needed both for `tool_meshConversorPar` (partitioning) and
  to run the `MeshElasticitySolver` / production solves.

Generated mesh, results, and log files (`*.hdf`, `*.h5`, `*.msh`, `*.log`,
etc.) are gitignored — regenerate them by running the pipeline rather than
expecting them to be present after a clone.

## Credits

This repository builds on several external projects. Full credit to their
authors:

- **[Construct2D](https://sourceforge.net/projects/construct2d/)** —
  Copyright © 2013–2018 Daniel Prosser, GPLv3. Vendored in full under
  `Construct2D/` (source, Makefiles, license, and docs as distributed
  upstream). Used unmodified as an external tool invoked by the pipeline,
  not linked into anything here.
- **[SOD2D](https://gitlab.com/bsc_sod2d/sod2d_gitlab)** — Copyright ©
  2022 Lucas Gasparino, Jordi Muela and Oriol Lehmkuhl (Barcelona
  Supercomputing Center), MIT License. Included as the
  `Construct2D_to_SOD2D/CFD_code/sod2d_gitlab` git submodule.
  - SOD2D itself uses **[GeMPa](https://gitlab.com/rickbp/gempa)**
    (Ricard Borrell, BSC) for mesh partitioning, and
    **[json-fortran](https://github.com/jacobwilliams/json-fortran)**
    (Jacob Williams) for JSON parsing, both as its own submodules.
- **`Construct2D_to_SOD2D/linear_mesh/p3d_to_gmsh.py`** — its Plot3D/Neutral
  Map File/Gmsh I/O layer (`read_chunk`, `NeutralMapFile`, `P3DfmtFile`,
  `GmshFile`) is adapted from
  **[p3d2gmsh](https://github.com/mrklein/p3d2gmsh)**, Copyright © 2015
  Alexey Matveichev, MIT License (full notice retained in the file header).
  The O/C-grid handling, wall/inlet/outlet classification, and periodic
  node-id aliasing on top of that base are original to this repo.
- **`Construct2D_to_SOD2D/linear_mesh/sod2d_tools/gmsh2sod2d.py`** —
  vendored from the SOD2D project's own tooling (same BSC credit as above).
- **[Gmsh](https://gmsh.info/)** (Christophe Geuzaine and Jean-François
  Remacle) — external dependency, invoked as a CLI; not vendored or
  redistributed here.
- The parametric arc-length wall-boundary placement in
  `MeshElasticitySolver.f90`'s `imposedDisplacement_elasticitySolverBufferSplineWall`
  was reverse-engineered from a working **[Nek5000](https://nek5000.mcs.anl.gov/)**
  reference case (`smooth_geom0` in a `naca_set_e448.usr` case file,
  Argonne National Laboratory / Nek5000 contributors) — an algorithmic
  reference, not vendored code.

## Development

Substantial parts of this repository — the mesh-generation pipeline scripts,
the wall-boundary elasticity fix in `MeshElasticitySolver.f90`, and this
documentation — were developed with assistance from
[Claude Code](https://claude.com/claude-code) (Anthropic), working under
Vishal Kumar's direction and review.

## License

Code original to this repository is released under the MIT License (see
[`LICENSE`](LICENSE)). Vendored components keep their own upstream licenses
as noted above and in their respective directories: `Construct2D/` is
GPLv3 (see `Construct2D/license/gpl.txt`), and the `sod2d_gitlab` submodule
is MIT (see its own `LICENSE`). These are separate programs used as
pipeline stages, not statically combined into one binary.
