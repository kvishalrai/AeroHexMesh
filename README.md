# AeroHexMesh

*Computational Aerodynamics with Hexahedral Meshes*

A framework for generating high-order hexahedral 3D meshes for
spectral-element CFD solvers, built around three independent front ends —
pick whichever matches the geometry and resolution you actually have:

```
Construct2D             Construct2D_to_<SOLVER>/
2D closed    ───────►   extrude, convert, partition, smooth the
curve .p3d/.nmf          wall boundary for the target solver

3D surface   ───────►   pyHyp/
mesh (wing,              hyperbolic-march a volume mesh outward
wall-to-wall)            from the surface

YZ cross-    ───────►   etaGrid_to_SOD2D/
section mesh             sweep it around an airfoil curve + a
(any quad grid)          C-grid wake, convert, partition, smooth
```

Each is a self-contained pipeline from its own starting geometry to a
smoothed, solver-ready 3D mesh — see its own README for the details of
that pipeline's stages. Although built and tested around airfoils,
nothing in the `Construct2D_to_<SOLVER>/` pipelines is airfoil-specific —
they work from any closed 2D curve Construct2D can mesh. For the math
behind every algorithm and how it's implemented, see
[`docs/`](docs/README.md).

## Repository layout

| Path | What it is |
|---|---|
| `Construct2D/` | Shared, vendored 2D grid generator, designed for airfoils but usable for other closed-curve geometries too (source + Windows binary). See its own [README](Construct2D/README.md) and [Credits](#credits). |
| `Construct2D_to_SOD2D/` | 2D grid → smoothed high-order mesh for [SOD2D](https://gitlab.com/bsc_sod2d/sod2d_gitlab). See its own [README](Construct2D_to_SOD2D/README.md). |
| `Construct2D_to_NEKRS/` | 2D grid → smoothed, spanwise-periodic 3D mesh for [NekRS](https://github.com/Nek5000/nekRS). See its own [README](Construct2D_to_NEKRS/README.md). |
| `pyHyp/` | 3D surface mesh → volume mesh by hyperbolic extrusion, for real (tapered/swept) 3D shapes rather than a spanwise-extruded 2D section. See its own [README](pyHyp/README.md). |
| `etaGrid_to_SOD2D/` | Arbitrary YZ cross-section mesh, swept around an airfoil curve (not just spanwise-translated) and into a C-grid wake, for [SOD2D](https://gitlab.com/bsc_sod2d/sod2d_gitlab). See its own [README](etaGrid_to_SOD2D/README.md). |
| `CFD_solvers/` | The CFD solver git submodules (SOD2D, Nek5000, NekRS) themselves, shared across every pipeline above that targets that solver — not a pipeline of its own. See [Credits](#credits) for what's included and where. |

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

Each pipeline directory has its own prerequisites and build instructions —
see its README (e.g.
[`Construct2D_to_SOD2D/README.md`](Construct2D_to_SOD2D/README.md),
[`pyHyp/README.md`](pyHyp/README.md),
[`etaGrid_to_SOD2D/README.md`](etaGrid_to_SOD2D/README.md)).

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
  [`CFD_solvers/sod2d_gitlab`](CFD_solvers) git submodule — shared at the
  repo root, since `etaGrid_to_SOD2D/` and `pyHyp_to_SOD2D/` build
  against this same copy too.
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
  [`CFD_solvers/{Nek5000,nekRS}`](CFD_solvers) git submodules — shared at
  the repo root, for any future NekRS-consuming pipeline; their own mesh
  tools (`gmsh2nek`, `genmap`, `re2torea`, `reatore2`, `n2to3`) drive
  most of `linear_mesh/`'s pipeline.
- **`linear_mesh/smooth_2D/naca_gen_spline_info.usr`**'s `smooth_geom0`/
  `smooth_geom` — Paul Fischer's (UIUC / Argonne National Laboratory)
  wall-smoothing algorithm (see above), adapted to run against this
  repo's own meshes and to support both O-grid and C-grid.
- **[p3d2nek](https://github.com/yslan/p3d2nek)** (YuHsiang Lan, Argonne
  National Laboratory) — a MATLAB-based alternative to this pipeline's
  Python/Gmsh-based `p3d_to_gmsh_nek.py` + `gmsh2nek` route.

### `pyHyp/` pipeline

- **[pyHyp](https://github.com/mdolab/pyhyp)** and
  **[cgnsutilities](https://github.com/mdolab/cgnsutilities)** — MDO Lab
  (University of Michigan), external dependencies built from source (not
  vendored); see [`pyHyp/README.md`](pyHyp/README.md#prerequisites) for
  the build recipe. `generate_volume_mesh.py` is a thin CLI wrapper
  around `pyhyp`'s own Python API.
- The bundled `example_m6_wing/m6_small.fmt` surface mesh is taken from
  pyHyp's own test suite.

### `etaGrid_to_SOD2D/` pipeline

- Reuses `Construct2D_to_SOD2D/linear_mesh/wall_spline.py` and
  `sod2d_tools/gmsh2sod2d.py`/`tool_meshConversorPar` unmodified (see
  above) — the sweep/wake/boundary-classification code around them is
  original to this repo.
- Adds a new wall-perturbation mode
  (`imposedDisplacement_elasticitySolverBufferWavyWall`, dispatched via
  `wavy_wall_amplitude_fraction`/`wavy_wall_wavenumber`) directly to
  SOD2D's own `MeshElasticitySolver.f90`, alongside its existing
  Fischer-derived spline mechanism (same BSC credit as above).

## References & Acknowledgements

- This work started as a collaboration between Argonne National Laboratory
  (ANL) and the Barcelona Supercomputing Center (BSC).
- Part of this work was published as: V. Kumar, A. Tomboulides, P. Fischer,
  and M. Min, "Delayed detached-eddy simulations of NACA wing sections using
  spectral elements," *Journal of Turbulence*, vol. 26, 2025.
  https://doi.org/10.1080/14685248.2025.2608679
- `etaGrid_to_SOD2D/`'s underlying concept (an unstructured cross-section
  swept with wall-normal/spanwise grid sizes proportional to the local
  Kolmogorov scale η) comes from: A. Rouhi, V. Kumar, W. Wu, M. Kozul, and
  O. Lehmkuhl, "Leveraging unstructured grids for direct numerical
  simulations of wall turbulence," under review, *Journal of Fluid
  Mechanics*. https://arxiv.org/abs/2605.01015
- Contributions from several BSC and ANL personnel are acknowledged:
  Yu-Hsiang Lan (ANL) -- see his GitHub, https://github.com/yslan/, for
  many Nek5000/NekRS-related tools -- Bedri Yagez (BSC), Jose Maria (BSC).
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
