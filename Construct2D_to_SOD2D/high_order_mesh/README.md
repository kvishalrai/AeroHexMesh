# Step 4: wall boundary smoothing

This is the last step of the pipeline described in
[`../linear_mesh/README.md`](../linear_mesh/README.md#the-algorithm).
The mesh coming out of step 3 has the right shape, boundary tags, and
partitioning — but its wall is still the jagged, faceted approximation of
the true airfoil curve. This step fixes exactly that.

## What this step actually does

`p3d_to_gmsh.py` + `gmsh2sod2d.py` (steps 1–3) place a high-order mesh's
*extra* boundary nodes (the ones beyond the coarse corner points) by
simple straight-line interpolation — so a smoothly curved wall ends up
looking like a series of short flat facets, especially visible with a
coarse mesh or near sharp curvature like the leading edge.

This step uses SOD2D's own `MeshElasticitySolver` to:

1. Fit a smooth curve (a cubic spline) through the wall's own corner
   points — the same coarse points Construct2D produced.
2. Snap every high-order wall node exactly onto that curve.
3. Elastically relax the rest of the mesh to match, so moving the wall
   nodes doesn't tangle or collapse the interior elements — think of the
   mesh as a stretchy material that gets gently reshaped around the
   corrected wall, rather than the wall nodes being moved in isolation.

Implementation lives in
`../CFD_code/sod2d_gitlab/src/lib_mainBaseClass/sources/MeshElasticitySolver.f90`,
in `imposedDisplacement_elasticitySolverBufferSplineWall` — this is the
**default** binding, so no rebuilding SOD2D is needed to use it.

Node placement on the curve is **parametric (arc-length)**, not a
nearest-point search: each high-order boundary node's target position is
computed directly from its own reference coordinate and the wall spline's
matched corner arc-lengths. This matters because a naive "find the
nearest point on the curve" approach can send two different mesh nodes to
the *same* point on the curve near high curvature (like the leading
edge) — the arc-length approach can't do that, by construction.

## Where the input files come from

None of the files in this folder besides `MeshElasticitySolver.json`,
`airfoil0.sh`, and `mn5_bind.sh` are meant to be hand-written — they're
produced upstream and copied in:

| File | Produced by |
|---|---|
| `naca_re50k-4.hdf` (mesh; gitignored, not tracked) | `linear_mesh/run_pipeline.py`'s partitioning step (`tool_meshConversorPar`) — copy from `linear_mesh/runs/<work-dir>/` |
| `naca_re50k_wall_spline.dat` | `linear_mesh/run_pipeline.py`, via `wall_spline.py`'s `build_wall_spline_table` — written automatically as `{airfoil}_wall_spline.dat` alongside the mesh, then copied here |
| `MeshElasticitySolver.json` | hand-written config for this step (see below) |
| `airfoil0.sh` | hand-written SLURM job script — runs `sod2d MeshElasticitySolver` from the `build_gpu` build |
| `mn5_bind.sh` | hand-written MN5 process/GPU binding wrapper, invoked by `airfoil0.sh` |

Files this folder should **not** keep tracked: anything `sod2d` itself
writes when you run it (`surf_code_*.dat`, `paramQual.txt`, `roughness.dat`,
logs, restart/results HDF5, SLURM `out.o`/`error.e`) — those are gitignored;
rerun the solver to regenerate them.

## JSON config

`MeshElasticitySolver.json` in this folder is a working example. Fields
specific to this step:

```json
"elasticity_problemType": "elasticity_fromALE",
"bouCodes": [
   {"id":1, "bc_type":"bc_type_non_slip_adiabatic_moving"},
   ...
],
"wall_spline_table_file": "naca_re50k_wall_spline.dat"
```

- `wall_spline_table_file` — required. Path to the compact corner +
  per-segment cubic-coefficient table written by `wall_spline.py` (see
  above). Format: line 1 is `n_corners`; then `n_corners` rows of `x y s`;
  then `n_corners-1` rows of `ax bx cx dx ay by cy dy` (segment `i`:
  corner `i` → corner `i+1`).
- `wall_spline_correction_scale` — optional, default `1.0` (no damping).
  Damps the per-call correction; not needed with the current parametric
  placement (a leftover safety net from an earlier, search-based approach
  that could need damping to avoid mesh collapse), but left
  JSON-configurable.
- `wall_spline_substeps` — optional, default `1` (single pass). Re-imposes
  and re-solves from the already-updated coordinates each iteration; same
  "not needed anymore, but available" status as the scale above.
- The wall's `bouCodes` id **must** map to `bc_type_non_slip_adiabatic_moving`,
  not plain `bc_type_non_slip_adiabatic` — the elasticity solver's
  Dirichlet-BC routine
  (`temporary_bc_routine_dirichlet_prim_meshElasticity` in
  `mod_bc_routines_meshElasticity.f90`) only applies a nonzero imposed
  displacement for that specific code. Revert it to
  `bc_type_non_slip_adiabatic` afterward for the actual flow run.

## Running it

`sbatch airfoil0.sh` from this folder (MN5, GPU partition — see the
script for module/queue setup). At `num_partitions>=3` on a large mesh,
two unrelated MPI crashes have been seen and are worked around directly
in `airfoil0.sh` (see its comments for the full explanation of each):
a UCX rendezvous-protocol segfault in halo-exchange, and a separate
segfault in Open MPI's `ompio` component reading a large mesh HDF5 file.

## How to tell if it worked

- The run prints a **max corner-match error** each time it runs; if a wall
  boundary corner doesn't coincide with any point in the spline table within
  `corner_match_tol` (mesh/table mismatch), the run stops with an error
  rather than silently producing a bad mesh.
- Quality (`minQ`/`maxQ`) is logged before and after the elasticity solve —
  it should stay essentially unchanged (the correction is sub-percent-of-chord).
- Load the resulting mesh in ParaView and confirm the wall boundary nodes sit
  on the smooth curve (most visible on a coarse/low-`porder` boundary where
  the faceting was previously obvious).
