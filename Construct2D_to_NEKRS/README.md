# Construct2D → NekRS

This is the **middle and final stage** of the pipeline described in the
[top-level README](../README.md). Construct2D gives you a 2D grid — a
cross-section of an airfoil, described by a grid of points. This pipeline
turns that 2D grid into the kind of mesh NekRS actually needs: a real 3D
mesh, with a smooth (not jagged) airfoil surface, ready to simulate airflow
over.

## The big picture

If you're new to this, here's the problem in plain terms before any tool
names come up:

1. **You have a 2D shape.** Construct2D's output (a `.p3d` file) is just a
   grid of points tracing the airfoil's cross-section and the empty space
   around it, out to some far boundary. Think of it as a fine mesh
   wrapped around a wing's cross-section, like a slice through a wing.
2. **NekRS needs a 3D mesh.** A real wing has span (length), so the 3D
   mesh is built by copying that 2D cross-section along a straight line —
   "extruding" it — and telling the solver the flow repeats identically at
   every cross-section (this is called **spanwise-periodic**; it's a
   standard simplification for studying a 2D airfoil shape in a 3D solver,
   valid as long as you're not modeling wingtip effects).
3. **The airfoil surface starts out jagged, not smooth.** The grid points
   from Construct2D only approximate the true, smooth airfoil curve — and
   the extra points the mesh format adds between them (for higher accuracy)
   get placed by simple straight-line interpolation by default, which
   makes a *curved* surface look like a series of tiny flat facets. Before
   running a real simulation, this pipeline bends those points back onto
   the true smooth curve — a **spline fit** through the grid's own corner
   points, evaluated at exactly the positions the solver needs. Skipping
   this step doesn't crash anything; it just makes the simulation less
   accurate near the wall, where accuracy matters most.
4. **Two grid shapes are supported: O-grid and C-grid.** Picture the
   airfoil from above. An **O-grid** wraps the mesh all the way around the
   airfoil in complete rings, like concentric circles — the mesh closes
   perfectly on itself. A **C-grid** instead follows the airfoil's surface
   and then peels off downstream in a wake region, so from above it looks
   like the letter C — better suited to capturing the wake behind the
   airfoil, but its "wall" (the airfoil surface itself) is now only part
   of one ring, not the whole thing, which the smoothing step has to
   handle differently (see [`linear_mesh/README.md`](linear_mesh/README.md)).
   **This repo currently ships a working C-grid example**; the tooling
   supports O-grid too, it's just not the example checked in.

Once you have that picture, the two stages below are just: *build and
smooth the mesh*, then *use it*.

```
2D grid from Construct2D  ──►  linear_mesh/       ──►  high_order_mesh/  ──►  run
   (.p3d + .nmf files)         build the 3D mesh,       apply the smooth      NekRS
                                fit the smooth wall      wall curve to the
                                curve on the 2D mesh     real 3D mesh
```

## Stages

1. **[`linear_mesh/`](linear_mesh/README.md)** — Step 1: convert the 2D
   grid to a mesh format NekRS's tools understand, tag its boundaries
   (wall / inlet / outlet), fit the smooth wall curve, and build the final
   3D mesh (still with the *old*, jagged wall — the fix from step 1 hasn't
   been applied to it yet).
2. **[`high_order_mesh/`](high_order_mesh/README.md)** — Step 2: a NekRS
   case (`example_nekrs/`) that, right before running, stamps the smooth
   wall curve computed in step 1 onto the real 3D mesh, then runs the
   actual airflow simulation.

The wall-smoothing algorithm itself (`smooth_geom0`) was originally
written by **Paul Fischer** (UIUC / Argonne National Laboratory) — see the
top-level README's [Credits](../README.md#credits) for the full
attribution.

## `CFD_solvers/Nek5000` and `CFD_solvers/nekRS`

Included as git submodules at the repo root's
[`CFD_solvers/`](../CFD_solvers) (shared with any future NekRS-consuming
pipeline) — external codebases this pipeline depends on but doesn't own
(see the top-level README's [Credits](../README.md#credits)). Both
stages above call small command-line tools built from Nek5000's source,
plus NekRS itself:

| Tool | Used by | Build notes |
|---|---|---|
| `Nek5000/bin/gmsh2nek` | `linear_mesh/` (2D mesh → Nek5000's own `.re2` format) | `tools/maketools gmsh2nek` |
| `Nek5000/bin/genmap` | `linear_mesh/smooth_2D/` | `tools/maketools genmap` |
| `Nek5000/bin/{re2torea,reatore2,n2to3}` | `linear_mesh/arglist.sh` (builds the 3D mesh) | `tools/maketools re2torea reatore2 n2to3` (`re2torea`/`reatore2` share one makefile) |
| NekRS itself | `high_order_mesh/example_nekrs/` | large GPU build — see nekRS's own [build docs](../CFD_solvers/nekRS/README.md), or use an existing site install as `example_nekrs/airfoil0.sh` does |

## Prerequisites

- Python 3 with `numpy`.
- [Gmsh](https://gmsh.info/) (an open-source mesh-generation tool), invoked
  as a command-line program.
- MPI, for running Nek5000/NekRS.
- Nek5000's tools and a NekRS install, built as described above.

**No mesh files, generated data, or compiled binaries are stored in this
repo.** Everything you'd need to regenerate (`*.re2` mesh files,
`geom_3x3*.dat` correction data, `*.ma2` partition maps, compiled
binaries, build caches, logs, restart files) is excluded on purpose — you
build them by running the pipeline yourself, following each stage's
README. This keeps the repo small and makes sure the instructions below
actually work, rather than silently going stale next to a cached file.

## Running it end to end

```bash
cd linear_mesh
# Step 1 -- see linear_mesh/README.md for the full walkthrough: convert
# the 2D grid, fit the smooth wall curve, and build the 3D mesh. Ends by
# writing the 3D mesh straight into high_order_mesh/example_nekrs/ and
# filling in its matching settings.

cd ../high_order_mesh/example_nekrs
# Step 2 -- see high_order_mesh/README.md: copy in the smooth-wall
# correction data from step 1, then submit the actual simulation.
sbatch airfoil0.sh
```
