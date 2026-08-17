# Boundary-condition convention (unified reference)

This convention is duplicated near-verbatim across three pipelines'
READMEs (`Construct2D_to_SOD2D/linear_mesh/`, `Construct2D_to_NEKRS/linear_mesh/`,
`etaGrid_to_SOD2D/`) — this page is the one place to check it, with the
per-pipeline variations called out explicitly.

## Physical-id convention (Gmsh-based pipelines)

| id | Name | Used by |
|---|---|---|
| 1 | WALL | All Gmsh-based pipelines |
| 2 | INLET | All Gmsh-based pipelines |
| 3 | OUTLET | All Gmsh-based pipelines |
| 4 | **Periodic** (`Construct2D_to_SOD2D`, `etaGrid_to_SOD2D`) **or Symmetry** (`pyHyp_to_SOD2D`) | See below |
| 109 | VolumeCode | All Gmsh-based pipelines |

**Why id 4 means two different things**: `Construct2D_to_SOD2D` and
`etaGrid_to_SOD2D` both have a genuine spanwise-periodic boundary and use
id 4 for it (`PERIODIC_ID` in `mesh_extrusion.py`/`boundary_conditions.py`).
`pyHyp_to_SOD2D` (a real 3D wing mesh, no spanwise periodicity at all)
reuses the same id for a symmetry plane instead (`SYMMETRY_ID` in
`cgns_to_gmsh.py`) — id 4 was simply the next free slot in a pipeline
that never needed genuine periodicity, not a deliberate overload of
meaning. **Check which one applies before reusing this convention in a
new pipeline.**

`Construct2D_to_NEKRS` uses the *same* id numbers for WALL/INLET/OUTLET
(`p3d_to_gmsh_nek.py` reuses the SOD2D pipeline's own classification
logic directly) but Nek5000's own `.rea`/`.par` format doesn't use a
"Periodic" physical-id the same way — spanwise periodicity there is
realized structurally, by `n2to3`'s own extrusion (see
[`Construct2D_to_NEKRS/linear_mesh.md` step 3](../Construct2D_to_NEKRS/linear_mesh.md#step-3-argosh--2d--3d-via-nek5000s-own-format-converters)),
not a Gmsh `Periodic Surface` directive.

## Inlet/outlet classification (shared math)

Every Gmsh-based pipeline that needs to split a farfield boundary into
inlet/outlet uses the **same** dot-product test, at each farfield
boundary point:

```
d = normalize( position(farfield) − position(wall) )    # local outward direction, same (i,k)
U_inf = (cos(AoA), sin(AoA))                              # free-stream direction; mesh itself is built at zero AoA
classify = 'inlet' if dot(d, U_inf) < 0 else 'outlet'      # entering vs. leaving
```

Implemented independently in `Construct2D_to_SOD2D/linear_mesh/p3d_to_gmsh.py`'s
`_classify_inlet_outlet()` and reused directly by
`Construct2D_to_NEKRS/linear_mesh/p3d_to_gmsh_nek.py`. `etaGrid_to_SOD2D`
doesn't need this test at all — its whole outer boundary is uniformly
INLET regardless of AoA, since it classifies purely from the
cross-section's own local `(Y,Z)` coordinates (see
[`etaGrid_to_SOD2D/README.md` §5](../etaGrid_to_SOD2D/README.md#5-boundary-classification)),
not from flow direction.

**Unverified caveat** (carried from the `Construct2D_to_SOD2D/linear_mesh/`
README): the sign convention above (flow tilts toward `+y` for positive
AoA) is a reasonable default, and behaves consistently across O-grid and
C-grid, but has **not** been checked against an independently-known-good
result. Worth confirming before trusting nonzero-AoA results.

For a C-grid, the two wake-end faces (the open ends of the "C") are
always OUTLET regardless of AoA — a separate special case, not run
through the dot-product test (their own "outward direction" isn't
meaningfully wall→farfield).

## `bc_type` strings (SOD2D JSON configs)

Distinct from the physical-id integers above — these are strings in
`MeshElasticitySolver.json`/`BluffBodySolverIncomp.json`/etc.'s own
`bouCodes` array, mapping a physical id to how SOD2D's Fortran should
treat that boundary for *this particular solver*.

| `bc_type` | Meaning | Used for |
|---|---|---|
| `bc_type_non_slip_adiabatic` | Fixed (non-moving) viscous wall | Production flow solves (wall shape already finalized) |
| `bc_type_non_slip_adiabatic_moving` | Viscous wall, **with an imposed displacement** | `MeshElasticitySolver`'s WALL boundary specifically — the *only* code its Dirichlet routine treats as "moving"; every other code (including plain `bc_type_non_slip_adiabatic`) is held fixed at zero displacement in that same routine |
| `bc_type_far_field` | Farfield/inlet-type boundary | Both elasticity and production configs |
| `bc_type_outlet_incomp` | Incompressible outlet | Both elasticity and production configs |
| `bc_type_symmetry` | Symmetry plane | `pyHyp_to_SOD2D`'s farfield id 4 (see above); also the fallback treatment for any unrecognized code in the elasticity Dirichlet routine |
| `bc_type_slip_wall_model` | Wall-modeled (not fully resolved) viscous wall | Compressible high-Mach production configs (`BluffBody3DSolver`'s own default template) |

**Important, easy to get backwards**: in a `MeshElasticitySolver.json`,
the wall's `bouCodes` entry *must* be `bc_type_non_slip_adiabatic_moving`,
not plain `bc_type_non_slip_adiabatic` — the elasticity solver's own
Dirichlet routine (`temporary_bc_routine_dirichlet_prim_meshElasticity`
in `mod_bc_routines_meshElasticity.f90`) only imposes a nonzero
displacement for that exact code; every other code (that plain one
included) is held fixed. Revert to plain `bc_type_non_slip_adiabatic`
for the subsequent production run, once the wall shape is finalized.
