# AeroHexMesh

*Computational Aerodynamics with Hexahedral Meshes*

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

## Repository layout

| Path | What it is |
|---|---|
| `Construct2D/` | Shared, vendored 2D grid generator, designed for airfoils but usable for other closed-curve geometries too (source + Windows binary). See its own [README](Construct2D/README.md) and [Credits](#credits). |
| `Construct2D_to_SOD2D/` | 2D grid → smoothed high-order mesh for [SOD2D](https://gitlab.com/bsc_sod2d/sod2d_gitlab). See its own [README](Construct2D_to_SOD2D/README.md). |
| `Construct2D_to_NEKRS/` | 2D grid → smoothed, spanwise-periodic 3D mesh for [NekRS](https://github.com/Nek5000/nekRS). See its own [README](Construct2D_to_NEKRS/README.md). |
| `pyHyp/` | 3D surface mesh → volume mesh by hyperbolic extrusion, for real (tapered/swept) 3D shapes rather than a spanwise-extruded 2D section. See its own [README](pyHyp/README.md). |

Additional `Construct2D_to_<SOLVER>/` pipelines may be added following the
same pattern.

## CFD Solvers

This framework builds meshes for these spectral-element CFD codes:

| Code | What it is | Docs | Paper |
|---|---|---|---|
| [SOD2D](https://gitlab.com/bsc_sod2d/sod2d_gitlab) | Incompressible & compressible (high-Mach); GPU-accelerated (BSC) | [Wiki](https://gitlab.com/bsc_sod2d/sod2d_gitlab/-/wikis/home) | [Gasparino et al. 2024](https://doi.org/10.1016/j.cpc.2023.109067) |
| [Nek5000](https://github.com/Nek5000/Nek5000) | Incompressible & low-Mach; CPU-based (ANL) | [Docs](https://nek5000.github.io/NekDoc/) | [Fischer 1997](https://doi.org/10.1006/jcph.1997.5651) |
| [NekRS](https://github.com/Nek5000/nekRS) | Incompressible & low-Mach; GPU-oriented successor to Nek5000 (ANL/UIUC/PSU) | [Docs](https://nekrsdoc.readthedocs.io/en/latest/) | [Fischer et al. 2021](https://arxiv.org/abs/2104.05829) |

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
  is based on the `smooth_geom0` algorithm by **Paul Fischer** (UIUC /
  Argonne National Laboratory) — the same algorithm used directly in the
  `Construct2D_to_NEKRS/` pipeline's `linear_mesh/smooth_2D/`, see below.

### `Construct2D_to_NEKRS/` pipeline

- **[Nek5000](https://github.com/Nek5000/Nek5000)** and
  **[NekRS](https://github.com/Nek5000/nekRS)** — Copyright © UChicago
  Argonne, LLC, BSD-3-Clause-style license. Included as the
  `Construct2D_to_NEKRS/CFD_code/{Nek5000,nekRS}` git submodules; their own
  mesh tools (`gmsh2nek`, `genmap`, `re2torea`, `reatore2`, `n2to3`) drive
  most of `linear_mesh/`'s pipeline.
- **`linear_mesh/smooth_2D/naca_gen_spline_info.usr`**'s `smooth_geom0`/
  `smooth_geom` — Paul Fischer's (UIUC / Argonne National Laboratory)
  wall-smoothing algorithm (see above), adapted to run against this
  repo's own meshes and to support both O-grid and C-grid.
- **[p3d2nek](https://github.com/yslan/p3d2nek)** (YuHsiang Lan, Argonne
  National Laboratory) — a MATLAB-based alternative to this pipeline's
  Python/Gmsh-based `p3d_to_gmsh_nek.py` + `gmsh2nek` route.

## References & Acknowledgements

- This work started as a collaboration between Argonne National Laboratory
  (ANL) and the Barcelona Supercomputing Center (BSC).
- Part of this work was published as: V. Kumar, A. Tomboulides, P. Fischer,
  and M. Min, "Delayed detached-eddy simulations of NACA wing sections using
  spectral elements," *Journal of Turbulence*, vol. 26, 2025.
  https://doi.org/10.1080/14685248.2025.2608679
- Contributions from several BSC and ANL personnel are acknowledged:
  YuHsiang Lan (ANL), Bedri Yagez (BSC), Jose Maria (BSC).
- V. Kumar acknowledges his AI4S fellowship within the Generación D
  initiative by Red.es, Ministerio para la Transformación Digital y de la
  Función Pública, for talent attraction (C005/24-ED CV1), funded by
  NextGenerationEU through PRTR.

## License

Code original to this repository is released under the MIT License (see
[`LICENSE`](LICENSE)). Vendored components keep their own upstream licenses
as noted above and in their respective directories: `Construct2D/` is
GPLv3 (see `Construct2D/license/gpl.txt`), the `sod2d_gitlab` submodule is
MIT (see its own `LICENSE`), and the `Nek5000`/`nekRS` submodules are
BSD-3-Clause-style (Copyright © UChicago Argonne, LLC; see their own
`LICENSE` files). These are separate programs used as pipeline stages, not
statically combined into one binary.

## Contact

Vishal Kumar — kumar14.rai@gmail.com — [LinkedIn](https://www.linkedin.com/in/vishal-kumar-69a32b4a/)
