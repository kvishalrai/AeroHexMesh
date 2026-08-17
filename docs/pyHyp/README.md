# `pyHyp/` — math and implementation

Deep-dive companion to [`pyHyp/README.md`](../../pyHyp/README.md).
`pyHyp` itself (the hyperbolic volume-mesh extrusion) is an external
MDO Lab tool, not part of this repo — this page covers the two pieces
that *are* this repo's own: `coarsen_surface.py` and `smooth_surface.py`,
both surface-preprocessing steps that exist to make pyHyp's own marching
better-conditioned.

## Hyperbolic marching (pyHyp itself, for context)

pyHyp extrudes a volume mesh outward from a 2D surface by integrating a
hyperbolic PDE system in a pseudo-radial marching direction, starting at
first-layer height `s0` and growing according to the JSON's own
pseudo-grid/smoothing parameters (`N`, `marchDist`, `epsE`/`epsI`,
`theta`, `volCoef`, etc. — see
[`pyHyp/README.md`](../../pyHyp/README.md)'s own option table and
pyHyp's own docs for the full parameter meaning). The two preprocessing
scripts below exist because that marching is sensitive to the *input*
surface's own point distribution and minimum edge length — not because
this repo re-implements any part of the marching itself.

## `coarsen_surface.py`

**Why**: `s0` (the first marching layer's height) has to stay well under
the input surface's own smallest edge, or the very first step
self-intersects near whatever feature is that fine. Coarsening the
surface raises that smallest-edge floor, letting `s0` — and so the whole
mesh — be coarser too.

### Index selection

```
coarsen_indices(n) = sorted( {0, 2, 4, ...} ∪ {n-1} )
```

Every other point, plus the last index unconditionally — so an
even-length axis still keeps *both* endpoints (needed for block-to-block
edges to stay matched after coarsening). Applied per structured axis
(`coarsen_block`), independently for I and J.

### Protected blocks and edge propagation

`--protect-blocks` exempts specific blocks (e.g. a small, high-curvature
wingtip cap) from coarsening entirely — halving resolution there can
leave pyHyp's marched mesh (and, worse, Gmsh's later order-elevation)
producing near-degenerate elements, even when the *linear* mesh's own
quality metric doesn't flag it. `find_protected_axes()` propagates this
protection to just the **shared edge** of any neighboring block —
detected by exact point coincidence (no CGNS/B2B connectivity needed) —
not that neighbor's whole extent: a structured block's two edges along
the same axis share that axis's own point count, so only the matching
axis needs to stay uncoarsened on the neighbor's side; its other axis
coarsens normally.

## `smooth_surface.py`

**Why**: diagnosed directly on pyHyp's own `uneven_sphere.fmt` example —
15× edge-length-ratio clustering concentrated at a single block corner,
which fed a proportionally oversized post-extrusion curving correction
(downstream, in `wall_surface_spline.py`) relative to the boundary
layer's own first-layer thickness. Fixing the *input* surface's own point
distribution is more robust than correcting for it downstream. No
existing MDO Lab tool does this (checked: pyHyp itself only
coarsens/reduces, `cgnsutilities`' `rebunch` is wall-normal only, pyGeo is
FFD shape deformation, pySurf is triangulated-surface intersection).

### What didn't work first

Documented directly in the script's own header, worth repeating since it
explains *why* the current method looks the way it does:

1. **Plain bilinear TFI** on raw 3D positions — a flat/ruled interpolant
   between edges sags off a curved surface (measured 16.6% of radius off
   a unit sphere).
2. **One-shot algebraic reparametrization** — stayed on the true surface,
   visibly improved the crude edge-length-ratio metric (15.6×→3.6×), but
   more than *doubled* the metric that actually matters (the downstream
   curving-correction magnitude) — a decoupled, one-shot remap is
   "uniform enough" in the crude sense without being smoothly,
   monotonically varying cell-to-cell, which is exactly what the
   correction (a measure of local curvature-per-cell) is sensitive to.
3. **Plain (unweighted) elliptic smoothing** — fixed the wobble, but the
   downstream correction got *worse* than the original file. The
   original's extreme clustering was *accidentally* curvature-adaptive
   right where genuine curvature is elevated; uniform redistribution
   throws that adaptivity away along with the wobble.

### Curvature-adaptive equidistribution (the method actually used)

Local curvature is estimated **intrinsically** from each block's own
fitted geometry (`wall_surface_spline.py`'s bicubic Hermite per-cell fit
— see
[`etaGrid_to_SOD2D/README.md`](../etaGrid_to_SOD2D/README.md) for a
related but distinct spline mechanism, and
[`pyHyp_to_SOD2D`'s own `wall_surface_spline.py`](../../pyHyp_to_SOD2D/wall_surface_spline.py)
for this one) — so it applies identically to a sphere, a wingtip cap, or
any other shape with no per-case tuning.

**Curvature** (`surface_curvature_at`), via finite-difference first/
second partial derivatives of the bicubic patch and the standard
first/second fundamental forms:

```
P_s, P_t              = ∂P/∂s, ∂P/∂t                    (central difference, step eps)
P_ss, P_tt, P_st       = ∂²P/∂s², ∂²P/∂t², ∂²P/∂s∂t      (central difference)
n = (P_s × P_t) / |P_s × P_t|                             (unit normal)

E = P_s·P_s,  F = P_s·P_t,  G = P_t·P_t                   (first fundamental form)
L = P_ss·n,   M = P_st·n,   N = P_tt·n                    (second fundamental form)

H = (EN + GL − 2FM) / [2(EG − F²)]                         (mean curvature)
K = (LN − M²) / (EG − F²)                                  (Gaussian curvature)
κ_total = sqrt(max(4H² − 2K, 0))                            (= sqrt(k1² + k2²), total curvature magnitude)
```

**Monitor field** (`node_monitor_field`) — per-node redistribution
weight, i.e. local spring stiffness:

```
κ_cell ← clip(κ_cell, 0, percentile_99.5(κ_cell))          # outlier-robust (finite-difference noise on an unevenly-sampled fit)
ref = median(κ_cell)
m_cell = min(1 + β·(κ_cell / ref),  max_ratio)              # β = --curvature-beta, default 2.0; max_ratio default 4.0
m_node = mean of the up-to-4 surrounding cells' m_cell
```

`β = 0` gives `m ≡ 1` everywhere (recovers plain, uniform elliptic
smoothing exactly — used as a regression check). `max_ratio` is a
**separate** safety valve independent of `β`/curvature magnitude, so a
curvature outlier can't runaway-cluster one block's points into
near-coincidence with an unrelated neighboring block's point — which
would break pyHyp's own connectivity-detection tolerance.

**Edge redistribution** (`redistribute_curvature_adaptive`) — standard 1D
curvature-weighted equidistribution: integrate the monitor along a fine
cubic-spline resampling of the edge, then invert to place points at equal
increments of *weighted* arc length —

```
W(s) = ∫₀ˢ m(σ) dσ                (cumulative weighted arc length, trapezoidal on a 4000-point fine sampling)
s_i  = W⁻¹( i · W(total)/(n−1) )   (equal steps in W, inverted back to s via interpolation)
```

monotonic by construction (higher monitor ⇒ smaller local spacing, but
never reorders points); `monitor ≡ 1` reproduces plain uniform arc-length
spacing exactly.

**Interior relaxation** (`_elliptic_relax`) — per iteration, each
interior point is damped-relaxed toward its curvature-**weighted** 4
structured neighbors (this is exactly the discrete equilibrium of a
spring network with local stiffness set by curvature — the same
spring-analogy principle used by mesh-motion methods like Batina 1990,
applied here to surface-point redistribution instead of volume-mesh
motion), then **reprojected** onto the true surface.

**Reprojection** (`project_to_surface`) — local Gauss–Newton refinement
of the point's own `(cell, s, t)` location in the fine per-cell Hermite
fit, starting from the *previous* iteration's own location (cheap, since
each step's displacement is small):

```
r = target − P(s,t)
J = [∂P/∂s, ∂P/∂t]                              (3×2 Jacobian, central difference)
δ = (JᵀJ + 10⁻¹²I)⁻¹ Jᵀr                          (regularized normal equations, 2×2 solve)
(s,t) ← (s,t) + δ
```

`s, t` are allowed to run outside `[0,1]` mid-solve — near a heavily
clustered region, a point can need to migrate across several cells
before settling; any overflow converts into an exact cell-index shift
(cells are C⁰-continuous at `s=1`/`t=1`) and the Newton solve resumes in
the new cell, repeated until the point stops changing cells. (A first
version clamped `s,t` to `[0,1]` and allowed only one cell-shift per
outer iteration — too weak: points needing several cells' worth of
migration got stuck at a cell boundary, showing up as a *new* wobble
mid-block instead of fixing the one at the corner.)

### Options that change the ground-truth surface

- **`--c1-seams`**: makes the surface every point reprojects onto
  genuinely C¹ (tangent-plane, not just position) continuous across
  block seams — without it, each block's own tangent estimate is
  one-sided (from its own interior only), which two independently-fit
  neighboring blocks generally disagree on (measured ~7.7–9.4° of
  surface-normal mismatch on this tool's own sphere test, present in the
  *original* file, untouched by any of the smoothing above since none of
  it touches how tangents are estimated). Fixed the same way
  `wall_surface_spline.py`'s own `build_ghost_points` fixes it for the
  separate downstream CFD-facing fit: find each edge's neighbor's own
  next-interior point (a "ghost" point, via the same connectivity-free
  coincidence detection used for corners/edges) and use a central
  difference across the seam instead of a one-sided one.
- **`--relax-corners`**: additionally lets shared corners move (a real
  geometry change, not just redistribution) — iteratively averages each
  shared corner with its edge-adjacent neighbors from *every* block that
  shares it (symmetric, order-independent, so every sharing block agrees
  on the result — required for pyHyp's own connectivity detection to
  still find them coincident afterward), reprojected each step onto
  whichever participating block's fine fit has the lowest index (a fixed,
  deterministic reference). Runs as a pre-pass before per-block
  edge+interior smoothing, which then treats the relaxed corners as its
  new fixed endpoints.

### Known topological limit

The corner where 3+ blocks meet (e.g. a cubed-sphere corner) is a
genuinely singular point in the block topology — pinned exactly at its
original position throughout ordinary curvature-weighted smoothing (which
concentrates points *toward* it, since curvature there is genuinely
elevated, but doesn't move it). Only `--relax-corners` changes this, and
only by averaging with neighbors, not by any curvature-driven placement
of its own.
