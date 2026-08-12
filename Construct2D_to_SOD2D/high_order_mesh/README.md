# NACA 4-digit boundary smoothing

`p3d_to_gmsh.py` + `gmsh2sod2d.py` place extra (non-vertex) high-order boundary
nodes on the airfoil wall by straight-line/GLL interpolation between the linear
Construct2D corner points — the curved wall boundary is faceted at high order,
not truly analytic. This step uses SOD2D's `MeshElasticitySolver` to snap every
wall node exactly onto the analytic NACA 4-digit curve and elastically relax the
interior mesh to match.

Implementation lives in
`../CFD_code/sod2d_gitlab/src/lib_mainBaseClass/sources/MeshElasticitySolver.f90`:
- `MeshElasticitySolver_readJSONNaca` reads the `"naca_airfoil"` JSON block.
- `imposedDisplacement_elasticitySolverBufferNACA` does the per-node projection:
  a damped/clamped Newton search (inline in the parallel loop, not a separate
  `!$acc routine seq` procedure call — that pattern isn't otherwise exercised
  anywhere in this file and caused a real `CUDA_ERROR_ILLEGAL_ADDRESS` crash)
  for the closest point on the analytic curve.

## JSON config

`MeshElasticitySolver.json` in this folder is a working example. Two things
matter beyond the usual `MeshElasticitySolver` keys:

```json
"elasticity_problemType": "elasticity_fromALE",
"bouCodes": [
   {"id":1, "bc_type":"bc_type_non_slip_adiabatic_moving"},
   ...
],
"naca_airfoil": {"digits":"0012", "chord":1.0, "x0":0.0, "y0":0.0, "sharp_te":true}
```

- The airfoil wall's `bouCodes` id **must** map to
  `bc_type_non_slip_adiabatic_moving`, not plain `bc_type_non_slip_adiabatic` —
  the elasticity solver's Dirichlet-BC routine
  (`temporary_bc_routine_dirichlet_prim_meshElasticity` in
  `mod_bc_routines_meshElasticity.f90`) only ever applies a nonzero imposed
  displacement for that specific code; every other wall code is pinned at zero
  displacement. Revert it to `bc_type_non_slip_adiabatic` afterward for the
  actual flow run.
- `"naca_airfoil"` is optional; omitting it (or any of its keys) defaults to
  NACA0012, chord=1, leading edge at the origin, sharp trailing edge —
  `digits` is a plain "MPXX" string (e.g. `"4412"`), `sharp_te:false` uses the
  standard finite/open trailing-edge coefficient instead.
- **`sharp_te` must match the mesh's trailing-edge *topology*, not the
  airfoil's nominal "digits" convention.** An O-grid closes its i-loop into a
  single shared TE node (`i=0` and `i=idim-1` are literally the same aliased
  node — see `p3d_to_gmsh.py`'s node-id aliasing) — i.e. the mesh itself is
  always geometrically sharp at the TE, even for a nominally "open" NACA0012.
  Using `sharp_te:false` there asks the shared TE node to move toward the
  *open*-TE target (a real, if small, ~0.0013·chord half-gap) while its two
  immediate neighbors barely move — a sign-flipping fold right at the TE that
  no choice of `E`/`nu` can un-invert (this is exactly what produced
  `minQ: 0.000000` for every trial in `assessBestElasticityParameters`).
  Always use `sharp_te:true` for an O-grid mesh. A C-grid's wake-cut keeps the
  two TE-approaching boundary nodes physically distinct (not aliased), so
  there `sharp_te` can genuinely match the airfoil's own TE type.

## Activating the subroutine

`imposedDisplacement_elasticitySolverBufferNACA` is **not** wired in by
default — it sits next to the file's other case-specific
`imposedDisplacement_elasticitySolverBuffer*` variants (same pattern as the
existing `...BufferCarlos`), so it doesn't change behavior for any other run
that uses `elasticity_fromALE`. To use it:

1. In `MeshElasticitySolver.f90`, find the `contains` block of the
   `MeshElasticitySolver` type (~line 89) and change
   ```fortran
   procedure, public :: initialBuffer => imposedDisplacement_elasticitySolverBuffer
   ```
   to
   ```fortran
   procedure, public :: initialBuffer => imposedDisplacement_elasticitySolverBufferNACA
   ```
2. Rebuild (`build_p4`, per `airfoil0.sh`).
3. Run once (`sbatch airfoil0.sh` from this folder).
4. Revert the binding (and the wall's `bc_type` in the JSON) before using this
   same solver for anything else.

## Verifying the result

Load the resulting mesh/results in ParaView and confirm the wall boundary nodes
sit on the smooth analytic curve (most visible on a coarse/low-`porder`
boundary where the faceting was previously obvious), and that element quality
stayed acceptable through the elasticity smoothing pass.
