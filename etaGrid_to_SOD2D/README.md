# etaGrid → SOD2D

A self-contained pipeline, independent of the `Construct2D_to_<SOLVER>/`
family: instead of spanwise-extruding a 2D Construct2D grid, it sweeps an
arbitrary **YZ cross-section mesh** (unstructured but quad, `η`/wall-normal
in Y, spanwise/periodic in Z) around an airfoil curve to build the 3D
volume, then extends it downstream in X for a C-grid wake. This lets the
wall-normal + spanwise resolution come from any YZ grid you can build —
not just a uniform spanwise translation — while still reusing SOD2D's own
wall-spline smoothing mechanism for the airfoil surface.

The underlying η-grid concept -- a cross-section with wall-normal (Y) and
spanwise (Z) grid sizes proportional to the local Kolmogorov scale η, thin
viscous-scaled near the wall and coarsening above it -- comes from Rouhi,
Kumar, Wu, Kozul & Lehmkuhl, "Leveraging unstructured grids for direct
numerical simulations of wall turbulence," under review, *Journal of Fluid
Mechanics* (<https://arxiv.org/abs/2605.01015>); this pipeline builds the
3D swept/wake mesh and SOD2D import around that concept.

## The big picture

1. **Start from a YZ cross-section.** A quad mesh in the (Y, Z) plane:
   Y = wall-normal distance (0 at the wall, growing outward to a far
   boundary), Z = spanwise/periodic. `make_demo_cross_section.py` builds a
   clean, transfinite example; any other quad YZ mesh works too, as long
   as its Y=0 boundary is the wall.
2. **Sweep it around the airfoil (body wrap).** `sweep_mesh.py` fits a
   cubic spline through the airfoil's own points (`airfoil_resample.py`
   can first re-resample that curve to a different point count / LE-TE
   spacing) and sweeps the YZ cross-section along it: each station's
   Y=0 row lands exactly on the true airfoil curve, and the cross-section
   is offset outward along the curve's own local normal — so the wall
   follows the airfoil chordwise while every Z-column keeps its own
   wall-normal resolution.
3. **Extend into a C-grid wake.** `build_wake.py` takes the body wrap's
   first/last stations (frozen at the trailing edge) and extrudes them
   straight downstream in X, merging the two wake sheets' wall (Y=0)
   nodes into a single centerline — the open, wake-cut C-grid topology.
4. **Classify and tag boundaries, convert, partition.**
   `boundary_conditions.py` classifies each 2D boundary edge (wall /
   outer / periodic) purely from the cross-section's own local (Y, Z)
   coordinates, then `run_pipeline.py` writes the raw mesh, tags it via a
   generated Gmsh `.geo` script (WALL=1, INLET=2, OUTLET=3, Periodic=4,
   VolumeCode=109 — the same convention as
   `Construct2D_to_SOD2D/linear_mesh/`), order-elevates it, converts it to
   SOD2D's HDF5 format (`gmsh2sod2d.py`), and partitions it
   (`tool_meshConversorPar`) — all in one command.
5. **Smooth the wall onto the true airfoil spline.** Order-elevation
   places the new high-order wall nodes by straight-line interpolation,
   faceting a smoothly curved airfoil. `elasticity_run/` uses SOD2D's own
   `MeshElasticitySolver`, driven by a 1D cubic spline fit through the
   wall's corner points (`Construct2D_to_SOD2D/linear_mesh/wall_spline.py`,
   reused unmodified), to snap every high-order wall node onto the true
   curve and elastically relax the rest of the mesh to match.
6. **Optional: a spanwise wavy-wall perturbation.** `elasticity_run2/`
   runs `MeshElasticitySolver` a second time on the already spline-smoothed
   mesh, adding a triangular (tent-shaped) spanwise displacement along
   the wall's local outward normal — amplitude tapering smoothly to zero
   at both trailing-edge ends of the arc-length parametrization and at
   both spanwise-periodic boundaries, never displacing inward. This is a
   genuine addition to `CFD_code/sod2d_gitlab`'s
   `MeshElasticitySolver.f90` (`imposedDisplacement_elasticitySolverBufferWavyWall`,
   dispatched via the new `wavy_wall_amplitude_fraction`/
   `wavy_wall_wavenumber` JSON fields) — the submodule needs rebuilding
   (`build_p2`) to pick it up. `wavy_wall.py` is a Python-side companion
   for independently recovering each wall node's arc-length position (via
   a nearest-point Newton search against the same airfoil curve) and
   verifying the resulting displacement field.
7. **Run production.** `sod2d_run/` runs SOD2D's `BluffBodySolverIncomp`
   on the final mesh — wall now fixed (`bc_type_non_slip_adiabatic`, not
   moving), inlet/outlet on the far/wake boundaries.

```
YZ cross-section   ──►  sweep over airfoil (Y)  ──►  extend into wake (X)  ──►  tag + convert + partition
  (any quad mesh)        body wrap, sweep_mesh.py      build_wake.py            run_pipeline.py
                                                                                        │
                                                                                        ▼
                                            spline-smooth wall  ──►  (optional) wavy-wall  ──►  run production
                                              elasticity_run/         elasticity_run2/          sod2d_run/
```

## Example input

`examples/airfoil.msh` + `examples/naca0012_sharp.dat` is the original YZ
cross-section + airfoil curve this framework was built against (as opposed
to `make_demo_cross_section.py`'s clean, synthetic transfinite grid).
**Known issue, not yet fixed:** `airfoil.msh` has ~40 degenerate
(near-zero-area) quads from the quad-subdivision step that removed its
original stray triangles — traced but not yet repaired (see
`make_demo_cross_section.py`'s own docstring for how this was diagnosed).
To be replaced with a repaired version.

## Prerequisites

- Python 3 with `numpy`, `scipy` (`CubicSpline`), and the `gmsh` Python
  package (also used as a CLI for order elevation).
- MPI and HDF5, for mesh partitioning and running SOD2D.
- SOD2D itself, built from `Construct2D_to_SOD2D/CFD_code/sod2d_gitlab`
  (this pipeline reuses that same submodule and its `build_p2` build —
  see its own README) — `MeshElasticitySolver.f90`'s wavy-wall addition
  needs that build rebuilt to pick it up (`build_p2/rebuild_p2.job`).
- `gmsh2sod2d.py` and `tool_meshConversorPar` from
  `Construct2D_to_SOD2D/linear_mesh/sod2d_tools/` (reused directly, not
  duplicated here).

**No mesh, results, or log files are stored in this repo** (`*.hdf`,
`*.h5`, `*.msh`, `*.log`, etc.) — regenerate them by running the pipeline
below rather than expecting them to already be there after a clone.

## Running it end to end

```bash
# 1. Build a YZ cross-section (or supply your own quad mesh)
python3 make_demo_cross_section.py --output-file demo_cross_section.msh

# 2. Resample the airfoil curve to the LE/TE spacing and point count you want
python3 airfoil_resample.py \
    --airfoil-file <naca0012.dat> --output-file naca0012_demo.dat \
    --n-upper 17 --n-lower 17 --le-spacing 0.005 --te-spacing 0.005

# 3. Sweep, tag, convert, partition -- see run_pipeline.py's own --help
python3 run_pipeline.py \
    --cross-section-file demo_cross_section.msh \
    --airfoil-file naca0012_demo.dat \
    --work-dir run_full --wake-s0 6.81e-3 --wake-total 20.0 --wake-layers 10 \
    --porder 2 --num-partitions 4

# 4. Smooth the wall onto the true airfoil spline -- copy the partitioned
#    mesh + wall spline table into elasticity_run/, then:
sbatch elasticity_run/run_elasticity.job

# 5. (optional) Spanwise wavy-wall perturbation -- copy the smoothed mesh
#    from step 4 into elasticity_run2/, then:
sbatch elasticity_run2/run_wavy_wall.job

# 6. Run production -- copy the final mesh into sod2d_run/, then:
sbatch sod2d_run/run_sod2d.job
```

## Directory layout

```
sweep_mesh.py               # cross-section loading, airfoil curve fit, body-wrap sweep
build_wake.py                # C-grid wake extension, full-mesh assembly + node merging
boundary_conditions.py       # wall/inlet/outlet/periodic boundary-face classification
airfoil_resample.py          # airfoil curve re-resampling (Vinokur stretching)
make_demo_cross_section.py   # example transfinite YZ cross-section generator
run_pipeline.py              # orchestrates sweep -> tag -> convert -> partition
check_raw_quality.py         # SICN quality check on the raw (pre-elevation) mesh
verify_boundary_conditions.py# independent geometric verification of BC tagging
wavy_wall.py                 # Python-side verification of the wavy-wall displacement
elasticity_run/               # step 5: wall spline smoothing (MeshElasticitySolver)
elasticity_run2/              # step 6: wavy-wall perturbation (MeshElasticitySolver)
sod2d_run/                    # step 7: production run (BluffBodySolverIncomp)
```
