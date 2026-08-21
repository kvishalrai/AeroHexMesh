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

## 3. Opposite-face-distance displacement (`opposite_face_disp_factor`, in progress — not yet usable)

A second, independent `MeshElasticitySolver` displacement mode,
dispatched instead of the spline mechanism above when
`opposite_face_disp_factor` (an integer) is nonzero. Instead of snapping
onto a fitted curve, it displaces each WALL node by an integer multiple
of the distance to its own **opposite-face node** — the node directly
across the same hex element, along the element's own local axis. That
search is *not* reimplemented here: it directly reuses
`mod_wall_model.f90`'s `evalEXAtFace`, the exact same "find the node
straight across this element" routine SOD2D's own LES wall models use
for exchange-distance evaluation.

```
u_buffer = scale · factor · L · n̂
```

where `L` and `n̂` (the node's own inward, into-the-fluid normal — same
`normalsAtNodes` convention `evalBoundaryNormals`/`normalFacesToNodes`
establish elsewhere in the solver) both come straight out of
`evalEXAtFace`'s per-face-node loop; `factor` is
`opposite_face_disp_factor`; `scale` reuses `wall_spline_correction_scale`
(default `1.0`) so the existing `wall_spline_substeps` re-imposes-and-
re-solves loop (§2 above, "Sub-stepping" comment in
`solve_elasticity_ALE`) works for this mode too, without any new
substepping machinery.

**Implementation note (`connecParOrig`, not `connecParWork`):** every
other call site of `evalEXAtFace` in this codebase passes
`connecParWork` — the periodic-slave-to-master-substituted connectivity
solver assembly needs for DOF consistency. This call is a one-off
*geometric* query, not part of that assembly, and `connecParWork`
actively breaks it: for a wall node on a periodic face whose opposite-
face node is itself a periodic slave, the substitution points the lookup
at its master's coordinates — a full periodic span away. Verified on a
real run: exactly the wall nodes sitting on the periodic boundary got
`L ≈ z_span` instead of `L ≈` the true near-wall spacing, 0 elsewhere.
`connecParOrig` preserves genuine per-node geometric identity and is the
correct array for this use (confirmed against `mod_fem_precond.f90`'s
own `lMaster` construction, which builds exactly that
`connecParOrig`-genuine-id → `connecParWork`-master-id map).

### 3.1 Known limitation: fails on the wall's own geometric singularities

Unlike §1–2's mechanism, this mode is **not yet usable at `factor=1`**
(the smallest nonzero value) on any real mesh tested — the elasticity
solve's own automatic quality-fallback search (§2.1) exhaustively fails
(`Not possible to find valid configuration`, ~130 `(E,ν)` combinations
tried) on both a NACA0012 C-grid-like periodic mesh
(`etaGrid_to_SOD2D`), a Construct2D O-grid mesh, and a genuinely 3D M6
wing mesh (`pyHyp_to_SOD2D`, no periodicity, blunt trailing edge) — so
this is not specific to one mesh topology, one solver build, or a sharp
trailing edge.

#### 3.1.1 Diagnostic sequence — what was ruled out first

Before finding the actual cause, four candidate mechanisms were tested
directly against dumped data and ruled out one at a time, rather than
guessed at:

1. **Simple normal-gap closure** (the wall node closing onto a
   stationary interior). Refuted: dumping each wall node's own position
   alongside its opposite-face node's position, both before and *after*
   `conjGrad_meshElasticity` (i.e. before the `minQ<1e-6` branch reverts
   it), shows the interior node follows the wall almost exactly —
   median ratio of interior-to-wall displacement magnitude `0.998`,
   median direction cosine `1.0000` (perfectly aligned). The gap barely
   changes for most nodes (median `1.7128e-3` before → `1.7092e-3`
   after); the near-wall "sandwich" translates together rather than one
   side staying put while the other closes in.
2. **In-plane (streamwise) shear** between neighboring wall stations.
   Refuted: for every pair of nearest-neighbor wall nodes at the same
   spanwise station, the signed area of the quadrilateral formed by
   (wall₁, wall₂, opposite₂, opposite₁) was computed before and after
   the solve (a 2D projected proxy for that local patch's Jacobian).
   0 of 6604 sampled quads flipped sign or went near-zero; the worst
   shrink ratio was `0.57` (a moderate compression, not a collapse).
3. **Spanwise (Z) variation** of the imposed field. Refuted, and
   confirmed *correct* instead: at matching `(x,y)` across every
   z-station sharing that position, the imposed displacement's
   magnitude spread was exactly `0.0000e+00` and direction spread
   `3.3e-16` (floating-point noise) — the field is perfectly
   z-invariant, exactly as required for a pure spanwise extrusion with
   no periodic-boundary bug remaining.
4. Only then: **directly locating the actual failing element(s)**, by
   reusing `computeQuality`'s optional `quality_pv(nelem,ngaus)` output
   (already computed every call, just never captured — `mu_e(:,:)` is
   reused as its storage) to dump, per rank, every element below a
   quality threshold: global id, `quality`, the worst Gauss point and
   its own quality, the element's center, and whether it touches the
   wall. This is what actually found the two failure signatures below —
   guessing at candidate mechanisms from the global `minQ` scalar alone
   was not converging.

#### 3.1.2 The two failure signatures found

- **Trailing-edge tangent degeneracy** (NACA0012 meshes): the
  per-element dump found exactly **12 elements**, all with `quality`
  and `mu_e` exactly `0.0`, all touching the wall, all clustered right
  at `x ≈ 1.0, y ≈ 0` — the trailing edge, on both the upper and lower
  surface, repeated identically at each spanwise station. Dumping the
  actual imposed displacement vector at the TE tip node directly
  confirms why: `(dx, dy, dz) ≈ (+1.71e-3, -1.4e-7, 0)` — displaced
  almost **purely tangentially** (downstream, +chord),
  not into the fluid at all. A sharp/near-sharp TE has no well-defined
  surface normal at the tip — the surface tangent itself is
  discontinuous there — so `evalEXAtFace`'s per-element search still
  returns *some* direction, but it degenerates toward the local
  element's own chordwise axis instead of a true perpendicular. That
  shears an already acutely-angled wedge element along its own edge,
  inverting its Jacobian regardless of `(E,ν)`. This is a **direction**
  problem, not a magnitude one — which is exactly why substepping
  (below) doesn't fix it.
- **Wing-root band collapse** (M6 wing mesh, `pyHyp_to_SOD2D`, no
  periodicity, blunt TE — so this is a genuinely different case from
  the NACA0012 meshes above): the same per-element dump found **88**
  exactly-zero-quality elements this time, all in the single row of
  wall elements closest to the root symmetry plane (`z ≈
  0.0091–0.0101`, the innermost spanwise station), spanning almost the
  *entire chord* (`x = 0.005` near the LE out to `x ≈ 0.81`) rather
  than one corner — a qualitatively different footprint from the
  TE-only NACA0012 failure. Checking the individual displacement
  vectors at these root nodes shows nothing obviously wrong in
  isolation: magnitude ~`7e-5` (right order), `dz` exactly `0.0`
  (correct), direction transitioning smoothly from `-X`-dominant near
  the LE tip to `Y`-dominant mid-chord — a locally sane vector field
  that still collapses the whole row it acts on. **Not yet
  root-caused** — plausibly either the root band's elements already
  being thin in the spanwise direction (squeezed between the fixed
  symmetry plane and the first interior station), or a symmetry-
  boundary connectivity wrinkle analogous to the periodic-face bug
  fixed above, applied here to a different boundary type. Investigation
  paused here, not resolved.

Substepping (`wall_spline_substeps > 1`, `wall_spline_correction_scale <
1`) does **not** fix either failure mode — tried at `scale=0.3`
(imposing only 30% of the current gap per step), and the exhaustive
`(E,ν)` search still fails identically on substep 1 alone. This is
expected for the TE case (a magnitude reduction doesn't fix a
displacement pointed the wrong way) and was the direct evidence that
ruled out "just needs smaller steps" as a fix.

Two debug-only diagnostics were added alongside this mode, guarded by
`opposite_face_disp_factor /= 0` (harmless dead weight otherwise):
`dump_wall_and_opposite_coords(tag)` (per-rank wall-node ↔ opposite-node
coordinate pairs, pre- and post-solve) and
`dump_bad_quality_elements(tag, threshold)` (per-rank list of low-quality
elements: global id, quality, worst Gauss point, element center, whether
it touches the wall) — both reusable for whoever picks the root-band
investigation back up.

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
