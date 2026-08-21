# `etaGrid_to_SOD2D/` — math and implementation

This page is the deep-dive companion to
[`etaGrid_to_SOD2D/README.md`](../../etaGrid_to_SOD2D/README.md), which
covers *what* the pipeline does and *how to run it*. Here we derive the
actual math behind every stage and walk through where it lives in code.

The underlying η-grid concept (wall-normal/spanwise grid sizes proportional
to the local Kolmogorov scale η) is from Rouhi, Kumar, Wu, Kozul &
Lehmkuhl, "Leveraging unstructured grids for direct numerical simulations
of wall turbulence," under review, *Journal of Fluid Mechanics*
(<https://arxiv.org/abs/2605.01015>) -- this page covers this repo's own
mesh-generation and SOD2D-import math built around that concept, not the
paper's own turbulence-statistics results.

## 1. Cross-section representation

A cross-section is a set of 2D points `yz` (shape `(N, 2)`, columns
`Y, Z`) plus quad connectivity `quad_conn` (shape `(M, 4)`, 0-based
indices into `yz`) and, optionally, leftover triangle connectivity
`tri_conn`. `Y` is wall-normal distance from the airfoil surface (0 at
the wall), `Z` is the spanwise/periodic coordinate.

`load_cross_section()` (`sweep_mesh.py`) reads this from a Gmsh `.msh`
file. If it finds triangles and `recombine_to_quad=True` (the default),
it does **not** use Gmsh's `recombine()` (pairing adjacent triangles into
quads on an existing mesh) — that was found to leave stray triangles
behind ("Perfect Match failed in quadrangulation"). Instead it forces
`Mesh.SubdivisionAlgorithm = 1` and calls `gmsh.model.mesh.refine()`:
**full all-quad subdivision**, splitting every 2D element via its own
edge midpoints. This guarantees zero triangles, at the cost of ~4×
resolution everywhere. (This is also the mechanism that introduced the
~40 degenerate quads in `examples/airfoil.msh` — see §7.)

## 2. Airfoil curve: fit and resample

### 2.1 Arc-length parametrization

Given an ordered sequence of points `pts` (from a Selig-format `.dat`
file, TE → upper surface → LE → lower surface → TE), the arc-length
coordinate is the cumulative chord length:

```
s[0] = 0
s[i] = s[i-1] + |pts[i] - pts[i-1]|         (arc_length_param, sweep_mesh.py)
```

`AirfoilCurve` (`sweep_mesh.py`) then fits a cubic spline `spline(s)`
through `(s, pts)` using `scipy.interpolate.CubicSpline` with
`bc_type="not-a-knot"` (the standard choice absent extra end-condition
information — matches slope and curvature of the two end segments rather
than forcing a particular end derivative or a periodic wrap, appropriate
here since the curve is open, not closed, at the TE).

- **Position**: `curve.position(s) = spline(s)`.
- **Tangent**: central-difference derivative,
  `t(s) = (spline(s+ε) - spline(s-ε)) / (2ε)`, normalized to unit length
  (`eps=1e-6`, clamped to `[0, s_total]` at the ends so it degrades to a
  one-sided difference there).
- **Outward normal**: the tangent rotated **−90°**:
  `n(s) = (t_y(s), −t_x(s))`. The sign (−90° vs +90°) is fixed
  empirically — `sanity_check_normal_sign()` checks that the normal at
  the airfoil's own point of maximum thickness on the upper surface
  points toward `+y`, i.e. away from the body.

### 2.2 Resampling: Vinokur two-sided stretching

`airfoil_resample.py` redistributes points along a segment (the upper or
lower surface) with **independently specified first and last spacing** —
something a single geometric-growth sequence can't do (it can only pin
one end). This is Vinokur's (1983) closed-form two-sided stretching
function.

Given `n` points, and desired first/last spacing `ds0`, `ds1` **as
fractions of the segment's total arc length** (`s0 = ds0·(n−1)`,
`s1 = ds1·(n−1)`, normalizing by the uniform spacing `1/(n−1)`):

```
A = sqrt(s1 / s0)
B = 1 / sqrt(s0 · s1)
```

- **If `B > 1`** (the common case — spacings smaller than uniform):
  solve `sinh(Δy)/Δy = B` for `Δy` by Newton's method, seeded either by
  Vinokur's own small-`(B−1)` series expansion (accurate for
  `B − 1 < 0.36`) or by `log(2B) + log(log(2B))` for larger `B`. Then

  ```
  u(ξ) = 0.5 · (1 + tanh(Δy·(ξ − 0.5)) / tanh(0.5·Δy)),   ξ = i/(n−1)
  ```

- **If `B < 1`** (spacings *larger* than uniform): the same construction
  with `sin`/`tan` in place of `sinh`/`tanh`, solving `sin(Δy)/Δy = B`.

- **If `B ≈ 1`**: degenerates to uniform spacing, `u(ξ) = ξ`.

Finally the normalized arc-length fraction is

```
t = u / (A + (1 − A)·u)
```

with `t[0], t[-1]` snapped to exactly `0, 1` to avoid float drift.
`resample_segment()` then evaluates `t·total_length` through a
`CubicSpline` fit of the *original* segment, so resampled points still
lie exactly on the true input curve — only their spacing changes.

`resample_airfoil()` splits the closed input loop at the point of
maximum distance from the TE (robust to camber — doesn't assume the LE
sits at `y=0`), resamples the upper (TE→LE) and lower (LE→TE) segments
**independently** with their own point counts `n_upper`/`n_lower` but the
*same* `le_spacing`/`te_spacing` targets, then concatenates them
(dropping the duplicated shared LE point).

## 3. Body-wrap sweep

`sweep_cross_section(yz, quad_conn, tri_conn, curve, s_values)`
(`sweep_mesh.py`) places station `k` (arc-length `s_values[k]`) at:

```
X(k, j) = P_x(s_k) + Y_j · N_x(s_k)
Y(k, j) = P_y(s_k) + Y_j · N_y(s_k)
Z(k, j) = Z_j
```

where `P(s) = curve.position(s)`, `N(s) = curve.outward_normal(s)`, and
`(Y_j, Z_j) = yz[j]` is the cross-section's own local coordinate for
node `j`. This is a **ruled-surface-style offset sweep**: at `Y_j = 0`
(the wall row) every station lands exactly on the true airfoil curve;
away from the wall, each node is carried along at a fixed offset
distance `Y_j` in the *local* normal direction at that station — so the
far boundary is itself a curved, offset copy of the airfoil, not a flat
plane. Z is untouched (pure translation along the sweep), which is what
makes the wall's periodic-boundary geometry trivially consistent (see
§5).

Hex connectivity between consecutive stations `k, k+1` is built directly
from `quad_conn` by concatenating the low/high station's own node
indices (`[lo[quad_conn], hi[quad_conn]]`), and likewise 6-node prisms
from any leftover `tri_conn`.

### 3.1 Orientation correction

The 2D cross-section's own quad winding, combined with the local sweep
direction, can produce a hex with **negative signed volume** depending on
how the two interact — not a topology bug, just a handedness ambiguity
that flips sign consistently across an entire swept region (confirmed by
sampling: whole regions came out 100% negative or 100% positive, never
mixed). `signed_hex_volume()` computes, per hex, corner vectors
`e1 = p1−p0`, `e2 = p3−p0`, `e3 = p4−p0` and returns
`dot(cross(e1, e2), e3)`. `ensure_positive_orientation()` samples up to
200 elements; if the majority are negative, it swaps every hex's node
order `[4,5,6,7,0,1,2,3]` (and the analogous swap for prisms) — moving
the "top" face to the "bottom" — which flips the sign without changing
the geometry.

## 4. Wake extension (C-grid)

`build_wake.py`'s `sweep_wake_side()` extrudes the cross-section
**straight** in `+X` from a single **frozen** `(p, n)` — the position and
outward normal evaluated once at the trailing edge (`s=0` for the upper
side, `s=s_total` for the lower side) — *not* continuing to follow the
airfoil curve:

```
X(w, j) = p_x + x_offset_w + Y_j · n_x
Y(w, j) = p_y + Y_j · n_y
Z(w, j) = Z_j
```

### 4.1 Geometric offsets

`x_offsets` (the downstream station positions) come from
`geometric_offsets(s0, total, n_layers)`: a geometric series
`offsets = [0, s0, s0(1+r), s0(1+r+r²), …]` reaching `total` in exactly
`n_layers` steps. The ratio `r` has no closed form for general
`n_layers`, so it's found by bisection on

```
f(r) = s0 · (r^n − 1) / (r − 1)         (= s0 · n for r → 1)
```

searching `f(r) = total` over `r ∈ [1, ∞)` (expanding the upper bracket
by ×1.5 until `f(r_hi) ≥ total`, then bisecting to `tol=1e-10`, capped at
200 iterations). The same function sizes both the wake's downstream
growth (§ here) and, in `make_riblet_cross_sections.py`'s exploratory
code, a cross-section's own near-wall Y-clustering.

### 4.2 Node merging

Because `Y=0` always maps to `p` regardless of the local normal, the
upper-wake and lower-wake sheets' **wall** nodes are geometrically
identical at every downstream station and get merged into one shared
centerline; `Y > 0` nodes stay distinct (the genuine, physically-real gap
between the two wake sheets, widening downstream due to the airfoil's
own small trailing-edge angle). `build_full_mesh()` assembles the global
node array as: body nodes, then upper-wake nodes for stations `1..end`
(station 0 reuses the body's own station-0 nodes), then lower-wake
*non-wall* nodes for stations `1..end` (station 0 reuses the body's last
station; wall nodes at every station reuse the upper-wake side's own wall
node index instead of duplicating) — then rebuilds hex/prism
connectivity against this global numbering and re-applies
`ensure_positive_orientation()` to the upper and lower wake blocks
independently (their handedness can differ from the body and from each
other).

## 5. Boundary classification

`boundary_conditions.py` classifies purely from the cross-section's own
local `(Y, Z)` values — no dependence on 3D position, which is what lets
one classification apply consistently to every station:

```
wall_mask     = |Y − Y_min| < tol      (Y_min = 0, the wall)
outer_mask    = |Y − Y_max| < tol      (far/outer boundary → INLET)
periodic0_mask = |Z − Z_min| < tol
periodic1_mask = |Z − Z_max| < tol
```

`extract_boundary_edges()` finds the cross-section's own outer loop(s) —
edges used by exactly one 2D element — then `classify_boundary_edges()`
buckets them by which mask **both endpoints** satisfy. `edges_to_quads()`
extrudes a classified edge between two stations into a boundary quad
face: `(a@lo, b@lo, b@hi, a@hi)`. `build_boundary_faces_for_stations()`
does this across every consecutive station pair for a given sequence
(body or one wake side), with `include_wall=False` for the wake sides
(their `Y=0` line is the merged interior centerline, not a physical
wall). `run_pipeline.py`'s `collect_boundary_faces()` combines these:
WALL from the body only; INLET/PERIODIC0/PERIODIC1 from body + both wake
sides; OUTLET from the two wake sides' own final-station 2D cap elements
directly (not a classified edge — the true downstream boundary).

Physical ids: `WALL=1, INLET=2, OUTLET=3, PERIODIC=4, VolumeCode=109` —
same convention as `Construct2D_to_SOD2D/linear_mesh/`'s
`p3d_to_gmsh.py`/`mesh_extrusion.py` (see
[`reference/boundary-conditions.md`](../reference/boundary-conditions.md)).

## 6. Tagging, order elevation, conversion, partitioning

`run_pipeline.py` writes the raw mesh as Gmsh **discrete entities**
(`write_raw_msh()`: one entity per boundary group plus the volume, all
referencing the same global node array, `MSH format 2.2`), then generates
a `.geo` script that:

```
Physical Surface("WALL",1)     = {wall_tag};
Physical Surface("INLET",2)    = {inlet_tag};
Physical Surface("OUTLET",3)   = {outlet_tag};
Physical Surface("Periodic",4) = {per0_tag, per1_tag};
Physical Volume("VolumeCode",109) = {volume_tag};
Mesh.ElementOrder = {porder};
Mesh 3;
Periodic Surface {per1_tag} = {per0_tag} Translate {0, 0, span_z};
```

and runs `gmsh {file}.geo -0` — this both order-elevates (placing new
mid-edge/mid-face nodes by straight-line interpolation between existing
corners — see §7 for why that's wrong at the wall) **and** writes the
native Gmsh `$Periodic` section gmsh2sod2d.py's own parser expects
(confirmed empirically: gmsh's own `-2.2` `$Periodic` output format
happens to align with that parser's line-skipping logic). Then
`gmsh2sod2d.py -p 4 -r {porder}` converts to SOD2D's HDF5 format, and
`tool_meshConversorPar` partitions it (`mpirun -np {n} tool_meshConversorPar
input.json`).

**Standing process rule**: always check the raw mesh's `minSICN` (Scaled
Inverse Condition Number — catches local element distortion/twisting
that a simple corner-volume sign check misses) *before* submitting the
order-elevation/conversion job — `check_raw_quality.py` does this via
`gmsh.model.mesh.getElementQualities(tags, "minSICN")`.

### 6.1 `uns_per_links` must always be `True` for this pipeline

`tool_meshConversorPar`'s periodic-face node matching has two modes,
selected by the `uns_per_links` flag in its own `input.json`: a fast
structured-mesh path (assumes periodic node pairs line up by simple
index correspondence) and a slower unstructured "Pseudo-Periodic
Elements" search (matches periodic faces geometrically, node by node).
`run_pipeline.py` originally hardcoded `uns_per_links: False`, which
crashed outright on a genuinely unstructured cross-section (this
pipeline's whole premise — quad, but not index-regular):

```
Error in generate_masSlaRankPar(..)! ...
Is your mesh unstuctured and periodic? Activate flag uns_per_links!
```

Went unnoticed on smaller/simpler test cross-sections (structured
enough by coincidence to pass either path) and only surfaced on a
larger, genuinely-unstructured cross-section
(`airfoil_vish.msh` — 883,935 nodes, 3,392 pseudo-periodic elements
once matched correctly). Fixed by making `uns_per_links: True`
unconditional in `run_pipeline.py`, since this pipeline's cross-sections
are never index-regular by design — not a per-mesh setting to tune.

## 7. Wall-spline smoothing (`elasticity_run/`)

Order elevation places new high-order wall nodes by **straight-line**
interpolation between existing corners — correct for a genuinely flat
segment, wrong for a curved one, turning the smooth airfoil into a series
of flat facets (most visible at high curvature, e.g. the leading edge, or
on a coarse mesh).

The fix reuses `Construct2D_to_SOD2D/linear_mesh/wall_spline.py`
unmodified: an **open, natural cubic spline** through the wall's own
corner points (`x_wall[i], y_wall[i]`, arc-length `s[i]` via the same
chord-length parametrization as §2.1), solved via a from-scratch
tridiagonal (Thomas-algorithm) solve — see
[`Construct2D_to_SOD2D/linear_mesh.md`](../Construct2D_to_SOD2D/linear_mesh.md#wall-spline-fit)
for the full derivation. The table this writes (`n_corners`, then
`(x, y, s)` per corner, then per-segment cubic coefficients
`(a,b,c,d)_x`, `(a,b,c,d)_y`) is read directly by SOD2D's Fortran
`MeshElasticitySolver`.

### 7.1 Fortran-side placement (parametric, not nearest-point)

`imposedDisplacement_elasticitySolverBufferSplineWall`
(`MeshElasticitySolver.f90`) does **not** search for the nearest point on
the spline — a Euclidean nearest-point search is a many-to-one map near
high curvature (e.g. two nodes straddling the LE can share a nearest
sample), which caused repeatable mesh collapses in an earlier version.
Instead, for each wall boundary face:

1. Match its 4 corners to the spline's own corner table by nearest
   neighbor (`corner_match_tol = 1e-4`; the run aborts loudly on a
   mismatch rather than silently producing a bad mesh).
2. Determine which axis (face-local I or J) is the "along-wall" arc-length
   direction — whichever pair of same-J (or same-I) corners has a
   *different* matched arc-length.
3. For every node on the face (corners **and** high-order nodes, via each
   node's own Gauss–Lobatto–Legendre reference coordinate
   `ξ_k ∈ [−1, 1]`, mapped through the face's own I/J index table):

   ```
   s_f = s_lo + (s_hi − s_lo) · 0.5 · (ξ_k + 1)
   dt  = s_f − s_iseg
   x_t = a_x[iseg] + dt·(b_x[iseg] + dt·(c_x[iseg] + dt·d_x[iseg]))    (cubic segment eval)
   y_t = a_y[iseg] + dt·(b_y[iseg] + dt·(c_y[iseg] + dt·d_y[iseg]))
   u_buffer = scale · (target − current)
   ```

4. `solveLinearElasticity` then treats `u_buffer` as an imposed Dirichlet
   displacement on WALL nodes and solves the linear elasticity PDE
   (parameters `E`, `nu` from the JSON) for the rest of the mesh, via
   conjugate gradient (`conjGrad_meshElasticity`) — think of the mesh as
   an elastic solid being gently reshaped around the corrected wall
   rather than the wall nodes being moved in isolation.

This mechanism is why `examples/airfoil.msh`'s ~40 degenerate quads
(from the all-quad subdivision in §1) matter: they land at the wall's own
`Y=0` row, so the same straight-line-interpolation defect this whole
mechanism exists to fix also corrupts the spline table's own input
corner points there — not yet repaired.

## 8. Wavy-wall perturbation (`elasticity_run2/`)

A **second**, independent `MeshElasticitySolver` pass, run on the
already spline-smoothed mesh, adding a spanwise perturbation:

```
Δ(s, z) = amplitude_fraction · chord · tent(wavenumber · z / z_span) · taper(s)
```

along the wall's own local outward normal at `s`.

### 8.1 Why not extend the single-curve spline mechanism

§7's mechanism assumes **one** wall curve, shared identically at every
Z — exactly what breaks for a Z-varying target. (A bicubic-Hermite
*surface* mechanism does exist, `wall_surface_spline.py`/
`imposedDisplacement_elasticitySolverBufferSplineSurface` — built for
pyHyp's genuinely-3D wing meshes — but it forces C¹ smoothness in
*both* parametric directions, which would round off a deliberately sharp
tent peak/trough. Rather than extend Fortran to dispatch on a per-Z-plane
table (real but heavier work, never done), this is computed directly in
the same Fortran routine using the wall corner/segment table purely as a
source of local geometry, not a relocation target.)

### 8.2 Chordwise taper

```
taper(s) = sin²(π·s/s_total)
```

Zero **with zero slope** (a smooth ease, not a kink) at both `s=0` and
`s=s_total` — both ends of the arc-length parametrization are the
trailing edge (`x=1` for a unit-chord airfoil; see §2.1) — peaking at the
leading edge (`s = s_total/2`) in between.

### 8.3 Spanwise tent wave

A **unipolar** (never negative — trough sits at the *original*,
undisplaced surface, not below it) triangular tent wave, period 1,
range `[0, 1]`:

```
φ = wavenumber · (z − z_min) / z_span
φ = φ − floor(φ)                       (wrap to [0, 1))
tent(φ) = 1 − |2φ − 1|
```

`tent(0) = tent(1) = 0`; `tent(0.5) = 1`. Because `wavenumber` is
declared `integer`, `φ` at `z = z_max` is exactly an integer, so its
fractional part is exactly 0 — the wave is *exactly* zero at both
spanwise-periodic boundaries, matching the mesh's own Z-periodicity.
`wavenumber = n` gives `n` teeth across the span. `chord` is
auto-computed as the spline table's own X-extent (`max − min`
`spline_corner_x`) — no separate parameter needed.

### 8.4 Local tangent/normal (Fortran side)

Reuses §7.1's corner-matching and arc-length interpolation *purely* to
get `s` and the local tangent/normal at each wall node — **not** to
relocate it (it's already correctly placed by §7's pass):

```
t_x = b_x[iseg] + dt·(2·c_x[iseg] + 3·dt·d_x[iseg])      (d/dt of the cubic segment)
t_y = b_y[iseg] + dt·(2·c_y[iseg] + 3·dt·d_y[iseg])
(t_x, t_y) ← normalize
(n_x, n_y) = (t_y, −t_x)                                  (matches AirfoilCurve.outward_normal's sign convention exactly)

u_buffer = amplitude_fraction · chord · tent(φ) · taper(s) · (n_x, n_y)
u_buffer_z = 0                                             (always -- displacement stays in the local X-Y plane at every Z)
```

Implemented as `imposedDisplacement_elasticitySolverBufferWavyWall`
(`MeshElasticitySolver.f90`), dispatched from the existing
`imposedDisplacement_elasticitySolverBufferSplineWall` via a check on the
new `wavy_wall_amplitude_fraction` JSON field (nonzero → dispatch;
default `0.0` → disabled, original behavior unchanged). New JSON fields:
`wavy_wall_amplitude_fraction` (real, fraction of chord),
`wavy_wall_wavenumber` (integer, periods across the span). Requires
rebuilding the `build_p2` SOD2D binary to pick up.

### 8.5 Python-side verification (`wavy_wall.py`)

An independent check, run against the *output* HDF5 mesh rather than
trusting the Fortran log: for each wall node (found geometrically — see
below), recover `s` via a **Newton/Gauss–Newton search** minimizing
`|curve.position(s) − point|²`:

```
f(s)  = dot(curve.position(s) − point, curve.tangent(s))
f'(s) ≈ dot(curve.tangent(s), curve.tangent(s))            (Gauss-Newton: drop the 2nd-derivative term)
s ← s − f(s)/f'(s)
```

seeded from the nearest of the airfoil's own control points
(`initial_s_guess()`), iterated up to 25 times (converges in a handful of
steps given a point already close to the curve). `wall_target_taper()`
then recomputes the same `sin²`/tent formulas in Python and checks the
node's actual displacement against them.

**Identifying wall nodes** (needed since the partitioned HDF5 duplicates
boundary nodes across ranks, and its own boundary-face node numbering is
non-trivial to map back to global row indices): purely geometric — a
node's `(X, Y)` distance to the nearest point on the airfoil curve is
either ~`1e-8` (a wall node, floating-point noise from the spline
correction) or jumps straight to `≥ 3.4e-3` (the first interior layer) —
a clean, five-order-of-magnitude gap, so a `1e-6` threshold cleanly and
robustly separates the two, with a bounding-box prefilter first for
efficiency (most of the mesh's ~250k nodes are trivially far from the
airfoil's own small `X, Y` footprint).

## 9. Production run (`sod2d_run/`)

`BluffBodySolverIncomp` (SOD2D's incompressible bluff-body solver) on the
final mesh, wall now `bc_type_non_slip_adiabatic` (fixed — no more
"moving" displacement; the perturbed shape from §8 is already baked into
the mesh's own coordinates). No math specific to this pipeline beyond
what SOD2D's own solver implements — see the top-level README's
[CFD Solvers](../../README.md#cfd-solvers) table for SOD2D's own docs/paper.
