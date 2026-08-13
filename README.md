# AeroHexMesh

*Aerodynamics with Hexahedral Meshes*

A framework for generating high-order, spanwise-extruded 3D meshes (O-grid or
C-grid) for spectral-element CFD solvers, starting from a 2D
[Construct2D](https://sourceforge.net/projects/construct2d/) grid. Although
built and tested around airfoils, nothing in the pipeline is airfoil-specific
— it works from any closed 2D curve Construct2D can mesh.

```
Construct2D             Construct2D_to_<SOLVER>/
2D closed    ───────►   extrude, convert, partition, smooth the
curve .p3d/.nmf          wall boundary for the target solver
```

Each `Construct2D_to_<SOLVER>/` directory is a self-contained pipeline from
the shared 2D Construct2D grid to a smoothed, solver-ready 3D mesh — see its
own README for the details of that pipeline's stages.

## References & Acknowledgements

- This work started as a collaboration between Argonne National Laboratory
  (ANL) and the Barcelona Supercomputing Center (BSC).
- V. Kumar acknowledges his AI4S fellowship within the Generación D
  initiative by Red.es, Ministerio para la Transformación Digital y de la
  Función Pública, for talent attraction (C005/24-ED CV1), funded by
  NextGenerationEU through PRTR.
- Part of this work was published as: V. Kumar, A. Tomboulides, P. Fischer,
  and M. Min, "Delayed detached-eddy simulations of NACA wing sections using
  spectral elements," *Journal of Turbulence*, vol. 26, 2025.
  https://doi.org/10.1080/14685248.2025.2608679
- Contributions from several BSC and ANL personnel are acknowledged:
  YuHsiang Lan (ANL), Bedri Yagez (BSC), Jose Maria (BSC).

## Repository layout

| Path | What it is |
|---|---|
| `Construct2D/` | Shared, vendored 2D grid generator, designed for airfoils but usable for other closed-curve geometries too (source + Windows binary). See its own [README](Construct2D/README.md) and [Credits](#credits). |
| `Construct2D_to_SOD2D/` | 2D grid → smoothed high-order mesh for [SOD2D](https://gitlab.com/bsc_sod2d/sod2d_gitlab). See its own [README](Construct2D_to_SOD2D/README.md). |
| `Construct2D_to_NEKRS/` | 2D grid → smoothed, spanwise-periodic 3D mesh for [NekRS](https://github.com/Nek5000/nekRS). See its own [README](Construct2D_to_NEKRS/README.md). |

Additional `Construct2D_to_<SOLVER>/` pipelines may be added following the
same pattern.

## Getting the code

```bash
git clone --recurse-submodules https://github.com/kvishalrai/AeroHexMesh.git
```

If you already cloned without `--recurse-submodules`:

```bash
git submodule update --init --recursive
```

## Prerequisites

Each `Construct2D_to_<SOLVER>/` pipeline has its own prerequisites and build
instructions — see its README (e.g.
[`Construct2D_to_SOD2D/README.md`](Construct2D_to_SOD2D/README.md)).

Generated mesh, results, and log files (`*.hdf`, `*.h5`, `*.msh`, `*.log`,
etc.) are gitignored — regenerate them by running the relevant pipeline
rather than expecting them to be present after a clone.

## Credits

This repository builds on several external projects. Full credit to their
authors:

### Shared across pipelines

- **[Construct2D](https://sourceforge.net/projects/construct2d/)** —
  Copyright © 2013–2018 Daniel Prosser, GPLv3. Vendored in full under
  `Construct2D/` (source, Makefiles, license, and docs as distributed
  upstream). Used unmodified as an external tool invoked by the pipeline,
  not linked into anything here.
- **[Gmsh](https://gmsh.info/)** (Christophe Geuzaine and Jean-François
  Remacle) — external dependency, invoked as a CLI; not vendored or
  redistributed here.

### `Construct2D_to_SOD2D/` pipeline

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
  vendored (with minor local tweaks) from `sod2d_gitlab`'s own
  `utils/gmsh2sod2d/gmsh2sod2d.py` (same BSC credit as above). Likewise
  `sod2d_tools/tool_meshConversorPar` is a compiled copy of
  `sod2d_gitlab/tool_meshConversorPar` — see
  [`Construct2D_to_SOD2D/README.md`](Construct2D_to_SOD2D/README.md#prerequisites).
- The parametric arc-length wall-boundary placement in
  `MeshElasticitySolver.f90`'s `imposedDisplacement_elasticitySolverBufferSplineWall`
  was reverse-engineered from a working **[Nek5000](https://nek5000.mcs.anl.gov/)**
  reference case (`smooth_geom0` in a `naca_set_e448.usr` case file, by
  **Paul Fischer** (UIUC / Argonne National Laboratory)) — an algorithmic
  reference, not vendored code at the time. That same case file is now
  directly used (not just referenced, and further generalized to handle a
  C-grid's open wall as well as the original O-grid) in the
  `Construct2D_to_NEKRS/` pipeline's `linear_mesh/smooth_2D/` — see below.

### `Construct2D_to_NEKRS/` pipeline

- **[Nek5000](https://github.com/Nek5000/Nek5000)** and
  **[NekRS](https://github.com/Nek5000/nekRS)** — Copyright © UChicago
  Argonne, LLC, BSD-3-Clause-style license. Included as the
  `Construct2D_to_NEKRS/CFD_code/{Nek5000,nekRS}` git submodules; their own
  mesh tools (`gmsh2nek`, `genmap`, `re2torea`, `reatore2`, `n2to3`) drive
  most of `linear_mesh/`'s pipeline.
- **`linear_mesh/smooth_2D/naca_gen_spline_info.usr`**'s `smooth_geom0`/
  `smooth_geom` — Paul Fischer's (UIUC / Argonne National Laboratory)
  original wall-smoothing algorithm (see above), renamed from
  `naca_set_e448.usr` and adapted to run against this repo's own meshes
  (originally O-grid only; generalized to also handle a C-grid's wall as a
  sub-range of the ring rather than the whole closed loop).
- **[p3d2nek](https://github.com/yslan/p3d2nek)** (YuHsiang Lan, Argonne
  National Laboratory) — a MATLAB-based Plot3D-to-Nek5000 mesh converter,
  an alternative to this pipeline's own `p3d_to_gmsh_nek.py` + `gmsh2nek`
  route. Explored as a possible route but not used in the pipeline
  documented here (needs MATLAB); referenced for anyone wanting a
  Gmsh-free path instead.

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
GPLv3 (see `Construct2D/license/gpl.txt`), the `sod2d_gitlab` submodule is
MIT (see its own `LICENSE`), and the `Nek5000`/`nekRS` submodules are
BSD-3-Clause-style (Copyright © UChicago Argonne, LLC; see their own
`LICENSE` files). These are separate programs used as pipeline stages, not
statically combined into one binary.
