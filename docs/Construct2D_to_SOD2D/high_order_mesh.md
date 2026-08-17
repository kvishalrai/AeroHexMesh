# `Construct2D_to_SOD2D/high_order_mesh/` — math and implementation

Deep-dive companion to
[`Construct2D_to_SOD2D/high_order_mesh/README.md`](../../Construct2D_to_SOD2D/high_order_mesh/README.md).
This step's Fortran mechanism (`MeshElasticitySolver.f90`) is shared
verbatim with `etaGrid_to_SOD2D/elasticity_run/` and
`pyHyp_to_SOD2D/`'s own elasticity step — this page derives the
mechanism itself; see
[`etaGrid_to_SOD2D/README.md` §7](../etaGrid_to_SOD2D/README.md#7-wall-spline-smoothing-elasticity_run)
for the same algorithm walked through with etaGrid's own worked example
and verification numbers.

## The problem

`p3d_to_gmsh.py` + `gmsh2sod2d.py` place a high-order mesh's *extra*
boundary nodes (beyond the linear mesh's coarse corner points) by simple
straight-line interpolation. A smoothly curved wall — the true airfoil
surface — becomes a series of short flat facets, most visible at high
curvature (leading edge) or on a coarse mesh. The wall's own *corner*
points are already exactly correct (they came straight from Construct2D);
only the newly-inserted high-order nodes between them are wrong.

## Fix, in two parts

1. **Where should each wall node go?** Fit a smooth curve through the
   wall's own corner points — an open, natural cubic spline
   (`wall_spline.py`, math in
   [`linear_mesh.md` §6](linear_mesh.md#6-wall-spline-fit-wall_splinepy)),
   parametrized by arc length.
2. **How do you move it there without breaking the mesh?** Snap every
   high-order wall node exactly onto that curve, then treat the rest of
   the mesh as an elastic solid and let it relax to accommodate the
   correction — SOD2D's own `MeshElasticitySolver`.

## 1. Fortran-side wall placement (parametric, not nearest-point)

`imposedDisplacement_elasticitySolverBufferSplineWall` places each
high-order wall node's target position **directly from its own reference
(GLL) coordinate**, not by searching for the nearest point on the spline
— a Euclidean nearest-point search is a many-to-one map near high
curvature (two nodes straddling the LE can share a nearest sample),
which caused repeatable mesh collapses in an earlier version of this
mechanism (reverse-engineered from a working Nek5000 reference,
`naca_set_e448.usr:smooth_geom0` — the same algorithm
`Construct2D_to_NEKRS/linear_mesh/smooth_2D/naca_gen_spline_info.usr`
implements for the Nek5000/NekRS pipeline).

Per wall boundary face, per node `k` with reference coordinate
`ξ_k ∈ [−1,1]` (from the Gauss–Lobatto–Legendre quadrature the spectral
element uses) mapped along the face's own "along-wall" axis:

```
s_f = s_lo + (s_hi − s_lo)·0.5·(ξ_k + 1)              # arc-length fraction, from the two matched corner arc-lengths
dt  = s_f − s_iseg
x_t = a_x + dt(b_x + dt(c_x + dt·d_x))                 # cubic segment eval, coefficients from wall_spline.py's table
y_t = a_y + dt(b_y + dt(c_y + dt·d_y))
u_buffer = scale · (target − current)                  # scale = wall_spline_correction_scale, default 1.0 (no damping)
```

The face's 4 corners are matched to the spline's own corner table by
nearest neighbor first (`corner_match_tol = 1e-4` — the run aborts
loudly, printing the max corner-match error, rather than silently
producing a bad mesh on a table/mesh mismatch), and which axis (face-
local I or J) is "along the wall" is determined from which pair of
same-row corners has a *different* matched arc-length.

## 2. Linear elasticity relaxation (`solveLinearElasticity`)

The imposed wall displacement `u_buffer` becomes a Dirichlet boundary
condition; the rest of the mesh solves the standard isotropic linear
elasticity equilibrium equation

```
∇·σ = 0,      σ = λ(∇·u)I + 2με(u),      ε(u) = ½(∇u + ∇uᵀ)
```

with the Lamé parameters computed directly from the JSON's `E`
(Young's modulus) and `nu` (Poisson's ratio):

```
μ = E / [2(1+ν)]
λ = E·ν / [(1+ν)(1−2ν)]
```

(`conjGrad_meshElasticity`, `mod_solver_meshElasticity.f90`). This is
discretized with the **same spectral-element basis** the flow solver
itself uses (GLL nodes, geometric factors `He`/`dNgp`/`gpvol`) and solved
by **preconditioned conjugate gradient** — either a block-Jacobi
preconditioner (`blck_diag_meshElasticity`, one 3×3 block per node) or a
simpler scalar-Laplacian-diagonal preconditioner, selectable via
`flag_useBlockJacobi`. Boundary conditions per node type: WALL
(`bc_type_non_slip_adiabatic_moving`) gets the full imposed `u_buffer`;
farfield/outlet-type codes are held fixed at zero displacement (see
`mod_bc_routines_meshElasticity.f90`'s dirichlet routine — any code not
explicitly one of the "moving" types falls to the same zero-displacement
`else` branch, so a boundary's exact `bc_type` string only matters for
whether it's the one *moving* boundary). `coordPar = coordPar + u` applies
the converged displacement.

### 2.1 Automatic quality fallback

After solving, `computeQuality()` reports `minQ`/`maxQ`. If `minQ < 1e-6`
(a genuinely inverted/degenerate element — the correction was too large
for the default `E, nu` to absorb without folding something over), the
solver **reverts** the displacement and calls
`assessBestElasticityParameters()`: a small grid search over `(E, nu)`
combinations (first a coarse 2-point sweep, then a finer 9-point sweep if
still failing), re-solving and re-checking quality each time, keeping
whichever parameters give the best `minQ`. Only if no combination reaches
`minQ ≥ 1e-6` does the run abort. This is why `E`/`nu` in
`MeshElasticitySolver.json` are reasonable *defaults*, not hard
requirements — the solver can self-correct for a mesh whose boundary-
layer aspect ratio makes the default parameters too stiff/soft.

## How to tell it worked

- **Max corner-match error**, printed every run — should be
  effectively 0 (floating-point noise); a table/mesh mismatch raises a
  hard error instead of a silently bad mesh.
- **`minQ`/`maxQ` before vs. after** — should stay essentially unchanged;
  the correction is sub-percent-of-chord, so a large quality swing
  signals something else is wrong.
- Load the result in ParaView and confirm the wall boundary sits on a
  smooth curve, most visible on a coarse/low-`porder` mesh where any
  remaining faceting is easy to spot.
