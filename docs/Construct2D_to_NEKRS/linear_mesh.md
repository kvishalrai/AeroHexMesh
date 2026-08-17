# `Construct2D_to_NEKRS/linear_mesh/` — math and implementation

Deep-dive companion to
[`Construct2D_to_NEKRS/linear_mesh/README.md`](../../Construct2D_to_NEKRS/linear_mesh/README.md).
Unlike the SOD2D pipeline, this one runs as **4 separate CLI steps**
rather than one orchestrator — Nek5000's own tools are old,
single-purpose command-line programs, not a Python-callable library, so
there's nothing to wrap into one function call the way
`gmsh2sod2d.py`/`tool_meshConversorPar` can be.

## Step 1: `p3d_to_gmsh_nek.py` — native 2D → linear Gmsh

Structurally the same idea as `Construct2D_to_SOD2D/linear_mesh/p3d_to_gmsh.py`
(node-id aliasing for O-grid closure / C-grid wake-cut fold, AoA-based
inlet/outlet classification — see
[`Construct2D_to_SOD2D/linear_mesh.md` §3.1, §5](../Construct2D_to_SOD2D/linear_mesh.md)
for the shared math, reused here), but operates on the **native 2D**
grid directly — no spanwise pre-extrusion in Python. The 3D extrusion
instead happens later, inside Nek5000's own format-conversion tools
(step 3), because `gmsh2nek` and Nek5000's `n2to3` already know how to
build a spanwise-periodic 3D case from a 2D one; duplicating that in
Python would just be redundant.

## Step 2: `smooth_2D/naca_gen_spline_info.usr` — Fischer's wall-smoothing algorithm

The **original** algorithm this whole repo's wall-smoothing family
derives from (SOD2D's `MeshElasticitySolver.f90` mechanism was
reverse-engineered from it — see
[`Construct2D_to_SOD2D/high_order_mesh.md` §1](../Construct2D_to_SOD2D/high_order_mesh.md#1-fortran-side-wall-placement-parametric-not-nearest-point)).
Written by Paul Fischer (UIUC/ANL), adapted here to run against this
repo's own meshes and support both O-grid and C-grid.

### 2.1 `spline`/`splint` — natural cubic spline (Numerical Recipes)

A from-scratch tridiagonal (Thomas-algorithm) natural cubic spline fit —
independently implemented in Fortran77 but the **same underlying math**
as `wall_spline.py`'s own (see
[`Construct2D_to_SOD2D/linear_mesh.md` §6](../Construct2D_to_SOD2D/linear_mesh.md#6-wall-spline-fit-wall_splinepy)):
natural boundary conditions (`y2(1) = 0`, and `qn = un = 0` forcing
`y2(n) = 0`), forward elimination then back substitution for the
second-derivative array `y2`, storing coefficients implicitly rather
than as explicit per-segment `(a,b,c,d)` — `splint()` evaluates directly
from `x, y, y2` via cubic Hermite-style blending:

```
h = xa[khi] − xa[klo]
a = (xa[khi] − x)/h,   b = (x − xa[klo])/h
y = a·ya[klo] + b·ya[khi] + [(a³−a)·y2a[klo] + (b³−b)·y2a[khi]]·h²/6
```

(binary search locates the bracketing segment `klo, khi` first — the
Numerical Recipes p.88-89 reference cited directly in the code comment).

### 2.2 Gathering wall points and evaluating at every GLL node

`smooth_geom0(jlev, elev)` collects each wall element's own face-1,
first-GLL-node `(x, y)` at a fixed radial level `jlev` and spanwise plane
`elev`, builds the chord-length arc-length array `s` (same cumulative-
distance construction used everywhere else in this repo), fits `x(s)`
and `y(s)`, then — critically — evaluates the spline **at every GLL
node's own arc-length position along its segment**, not just at the
element corners:

```
sf(j, e) = s(i-1) + ds·0.5·(zgm1(j,1) + 1)      # zgm1 = the GLL reference coordinate ∈ [-1,1]
```

directly overwriting `xm1(i,jlev,1,eidx)`/`ym1(...)` — Nek5000's own
mesh-coordinate arrays. This is the **same conceptual move** as SOD2D's
`sf = s_lo + (s_hi−s_lo)·0.5·(ξ_k+1)` (§7.1 there) — map each high-order
node's own reference coordinate linearly into arc-length between its
bounding corners, then evaluate the fitted curve there — independently
arrived at in a different codebase, which is some validation that it's
the *right* way to place high-order nodes on a fitted curve.

O-grid vs. C-grid handling: for an O-grid (`nwall == n_ring`, the wall is
the *entire* ring), point `n1` is a manual duplicate of point 1, closing
the spline's fitting sequence into a loop (same "not-adjacent-so-not-
degenerate" reasoning as `wall_spline.py`'s own O-grid handling). For a
C-grid, `n1` is read from the real element just past the wall range
(the wake-cut/farfield side) — for a sharp TE, `p3d_to_gmsh_nek.py`'s own
wake-cut aliasing already makes that point coincide with point 1
(closing the curve with no special-casing needed here); for a blunt TE
it's a genuinely distinct end point, and `spline`/`splint`'s natural
(zero second-derivative) end conditions handle an open curve correctly
as-is.

### 2.3 `blend_surface`/`laplaceh` — propagating the correction into the volume

Given the **displacement** `du` the spline correction just computed at
the wall (`new position − old position`), `blend_surface(du)` extends it
into the rest of the mesh by solving a **weighted Laplace equation** —
mechanically different from SOD2D's full linear-elasticity PDE (§ below),
though conceptually the same job (relax the interior around a corrected
boundary without folding anything over).

**Boundary-layer protection**: rather than a uniform diffusivity, `h1`
(the Helmholtz operator's diffusion coefficient) is boosted near the
wall:

```
δ = volume(wall-adjacent elements) / surface_area(wall-adjacent elements)   # avg boundary-layer element thickness
δ_p = 5·δ                                                                    # "protected" thickness
h1(x) = 1 + 9·exp( −(dist_to_wall(x) / δ_p)² )                              # → 10 at the wall, → 1 far away
```

so near-wall elements resist deformation ~10× more than far-field ones —
the same *goal* as SOD2D's optional `flag_useVariableYoung_helem`
(stiffen elements with small `helem`), achieved by a spatially-varying
Helmholtz coefficient instead of a spatially-varying Young's modulus.

**Solving `∇·(h1∇u) = 0` with `u = du` on the wall**: `laplaceh` uses the
standard Dirichlet-lift decomposition `u = u0 + ub` (`ub` = boundary data
extended trivially via `dsavg`, i.e. duplicate-node-averaged onto every
element sharing that point; `u0 = 0` on Dirichlet boundaries):

```
r = -M · Ā·ub          (Ā = the Helmholtz operator applied to ub; M = the interior/Dirichlet mask)
  Ā u0 = r             (solved for u0 via hmhzpf -- a preconditioned CG-type Helmholtz solve)
u = u0 + ub
```

`Ā` here is `h1·∇² + h2·I` with `h2 = 0` — a pure (weighted) Laplace
operator, not the full elasticity tensor SOD2D's `conjGrad_meshElasticity`
solves; there's no shear/Poisson coupling between coordinate components
(`x` and `y` displacements are smoothed independently, each its own scalar
Laplace solve).

## Step 3: `arglist.sh` — 2D → 3D via Nek5000's own format converters

Chains three of Nek5000's single-purpose CLI tools:

```
re2torea   # binary .re2 -> ASCII .rea (Nek5000's older text mesh format)
n2to3      # 2D .rea -> 3D .rea: performs the actual spanwise extrusion + periodicity
reatore2   # ASCII .rea -> binary .re2 (back to the format NekRS/genmap read)
```

i.e., the spanwise extrusion genuinely happens *inside* `n2to3`, not in
any Python code in this repo — `arglist.sh` just drives the three-tool
pipeline and fills in matching mesh-size settings in the target NekRS
case.

## `compute_case_params.py` — deriving case parameters from the NMF directly

Rather than running the mesh all the way through `gmsh`/`gmsh2nek` and
inspecting the output just to find basic element counts,
`compute_wall_params()` reads them directly out of the **native 2D
NMF's own `VISCOUS` entry** (`f1 == 3`, excluding the C-grid's
`ONE_TO_ONE` wake-cut pairing of that same face) — works for O-grid or
C-grid alike, no `--mesh-type` flag needed:

```
ring_width          = idim − 1                       # NUMBER_ELEMENTS_X, both .usr styles
radial_rings         = jmax − 1                        # NUMBER_ELEMENTS_Y, step-2 .usr meaning ONLY
wall_start_element   = VISCOUS.s1                       # 1-based, matches Fortran element numbering
wall_element_count   = VISCOUS.e1 − VISCOUS.s1
```

**The one genuinely confusing subtlety** (flagged explicitly in both the
script's own docstring and the top-level README): `NUMBER_ELEMENTS_Y`
means two *different* things depending which `.usr` it's patched into —
in step 1's `.usr` (Fischer's spline-fit case) it's how many near-wall
rings to actually run the spline correction over (**not** derived from
the mesh; must be passed explicitly via `--smoothing-levels`, since
silently defaulting it to the full radial count would smooth every ring
instead of just the near-wall one(s)); in step 2's `.usr`/`.par` it's the
mesh's real, full radial ring count (auto-derived, patched
automatically). Same variable name, same-looking value in some cases,
genuinely different quantities — this is exactly the kind of thing worth
double-checking by hand before trusting a patched case file.
