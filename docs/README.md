# AeroHexMesh — technical documentation

This is the deep-dive companion to the repo's own per-directory
READMEs. Every directory in AeroHexMesh already has a README covering
**what it does and how to run it** (algorithm overview, prerequisites,
usage, file provenance) — this `docs/` tree instead covers **the actual
math behind every algorithm and how it's implemented**, organized to
mirror the repo's own top-level folder structure. Read the relevant
top-level README first if you haven't; these pages assume that context
and go deeper rather than repeating it.

## By folder

| Folder | Docs |
|---|---|
| `Construct2D/` | [`Construct2D/README.md`](Construct2D/README.md) — batch-mode driving, auxiliary scripts |
| `Construct2D_to_SOD2D/` | [`linear_mesh.md`](Construct2D_to_SOD2D/linear_mesh.md) — extrusion, NMF remap, node aliasing, AoA classification, wall-spline fit, Gmsh/SOD2D export<br>[`high_order_mesh.md`](Construct2D_to_SOD2D/high_order_mesh.md) — arc-length wall placement, linear elasticity PDE |
| `Construct2D_to_NEKRS/` | [`linear_mesh.md`](Construct2D_to_NEKRS/linear_mesh.md) — Fischer's `smooth_geom0` algorithm, weighted-Laplace boundary blending, Nek5000 format-conversion chain<br>[`high_order_mesh.md`](Construct2D_to_NEKRS/high_order_mesh.md) — interpolating the 2D correction onto the real 3D mesh |
| `pyHyp/` | [`pyHyp/README.md`](pyHyp/README.md) — surface coarsening, curvature-adaptive elliptic surface smoothing |
| `etaGrid_to_SOD2D/` | [`etaGrid_to_SOD2D/README.md`](etaGrid_to_SOD2D/README.md) — sweep geometry, Vinokur stretching, wake growth series, boundary classification, wavy-wall tent function |

## Cross-cutting reference

| Page | Covers |
|---|---|
| [`reference/boundary-conditions.md`](reference/boundary-conditions.md) | The physical-id convention and `bc_type` strings — duplicated near-verbatim across 3 pipeline READMEs, unified here with the per-pipeline variations called out |
| [`reference/cluster-setup.md`](reference/cluster-setup.md) | The `AEROHEXMESH_*_MODULE_SETUP` env-var pattern, SLURM job-script conventions, MPI-crash workarounds, login-node limits |
| [`reference/troubleshooting.md`](reference/troubleshooting.md) | Every known caveat/gotcha in the repo, in one place |

## A recurring theme worth knowing up front

Three genuinely different wall-smoothing mechanisms exist across this
repo, solving the same underlying problem (order elevation places new
high-order boundary nodes by straight-line interpolation, faceting a
curved wall) three different ways:

1. **SOD2D's `MeshElasticitySolver`** (`Construct2D_to_SOD2D/`,
   `etaGrid_to_SOD2D/`, `pyHyp_to_SOD2D/`): snap onto a fitted spline
   parametrically (by arc length, not nearest-point), then relax the
   rest of the mesh via a genuine linear-elasticity PDE (Lamé parameters
   from `E`/`nu`, spectral-element discretization, conjugate gradient).
2. **Nek5000/NekRS's `smooth_geom0`/`blend_surface`**
   (`Construct2D_to_NEKRS/`): the same parametric arc-length placement
   idea, independently implemented — but relaxes the interior via a
   **weighted Laplace** (not full elasticity) solve, with the weight
   boosted near the wall to protect the boundary layer.
3. **Curvature-adaptive elliptic surface smoothing** (`pyHyp/`): a
   different problem entirely — smooths a *surface mesh's own point
   distribution* before extrusion (not a volume mesh's boundary after
   order elevation), via curvature-weighted equidistribution and a
   monitor-weighted Laplacian with Newton reprojection.

Worth reading [`Construct2D_to_SOD2D/high_order_mesh.md`](Construct2D_to_SOD2D/high_order_mesh.md)
and [`Construct2D_to_NEKRS/linear_mesh.md`](Construct2D_to_NEKRS/linear_mesh.md#step-2-smooth_2dnaca_gen_spline_infousr--fischers-wall-smoothing-algorithm)
side by side if you're curious how two independently-arrived-at
solutions to the same parametric-placement problem compare.
