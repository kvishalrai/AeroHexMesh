# Envelope pipeline for the 30P-30N three-element airfoil: plan

Status: **steps 0–3 implemented and working on a coarse envelope (raw surface
points, no refinement yet); steps 4–6 not started.** See "Progress" below.
First written 2026-10-04 as a
handover so the work can continue with another agent; tightened the same day
after checking it against the Construct2D source, the existing
`xcut_pipeline/`, the `Construct2D_to_SOD2D` converter and the raw geometry
(section "Review findings" at the end gives the evidence). The target sketch is
next to this file: `30p30n_envelope_sketch.pdf`.

## Progress (2026-10-04)

Done, not yet committed:

- **Step 0.** The XCUT edits were copied from `nek-wing` into
  `Construct2D/src/` (`edge_grid`, `menu`, `surface_grid`,
  `hyperbolic_surface_grid`, `vardef`; the `nek-wing` copy was checked to be
  a strict superset of this repo's). `util.f90` now writes the p3d and stats
  at `es25.16`. Rebuilt and checked: the OAT15 XCUT curve comes back
  bit-identical and at true scale.
- **Step 1** `step1_build_envelope.py`, **step 2** `step2_run_construct2d.py`,
  **step 3** `step3_check_mesh.py`, shared helpers in `env_common.py`.

Commands (build Construct2D first with `make` in `Construct2D/`):

```
python3 step1_build_envelope.py --out <dir>/30p30n_env_v1.dat
python3 step2_run_construct2d.py --curve <dir>/30p30n_env_v1.dat \
    --construct2d ../construct2d --out-prefix <dir>/30p30n_env_v1_c2d
python3 step3_check_mesh.py --meta <dir>/30p30n_env_v1_c2d.meta.json \
    --out-prefix <dir>/30p30n_env_v1_check
```

Result with the defaults (JMAX=100, RADI=15, YPLS=1, RECD=1.7e6): a 696 × 100
mesh with 0 folded cells, the wall row bit-identical to the curve, and a worst
cell angle of 41°. The same holds with `--far-end free` and the fixed
Construct2D (lines ending near y=-0.52, 0 folded cells, worst angle 41.6°). Step 3 still fails one gate: xi-growth reaches 1.25–1.45.
The worst values are in the outer field (j ≈ 80–90); near the walls (j ≤ 30)
the worst is 1.24, on the lower downstream line just off the flap TE. The
wall row itself has ratios up to 1.16, inherited from the raw spacing.

What the first runs found (now built into the scripts):

- **C10 Open cut: fixed in Construct2D.** The C-grid hyperbolic marcher
  solved only i=1..imax-1 and set column imax to the mirror of column 1 about
  y=0 at every level (`x(imax)=x(1)`, `y(imax)=-y(1)`). The XCUT patch did not
  change this, and OAT15 hid it because its cut sat on y=0. With the
  downstream lines ending near y=-0.52, the first run had 703 folded cells.
  `hyperbolic_surface_grid.f90` now marches i=1 and i=imax independently
  when the curve is open (CGRD, NWKE=0, two distinct end points). Each end
  column is a constant-offset copy of its own neighbour column, with the
  same optional XDAMP pull. Checked: the closed OAT15 XCUT curve gives a
  byte-identical mesh to before; the envelope with lines ending near y=-0.52
  now has 0 folded cells. Step 1's `--far-end free` (default) relies on the
  fix. `--far-end mirror` keeps the old y=0 shape for builds without it.
- **Open curves for OAT15 too.** xcut step 1 has `--open` (no closure point),
  and xcut steps 3–5 then run with `--nelm-clean 0`. Compared with
  closed-then-trimmed, the near-wall mesh is the same (within 3e-5 for
  j ≤ 20). The outflow cut columns stay nearly vertical (x drift 0.01–0.10
  against 1.5), and the worst cell angle in the five columns at each end is
  87° against 81°.
- **Construct2D's .nmf is wrong for an open curve.** It still welds i=1 to
  i=imax with a ONE_TO_ONE and marks all of j=1 VISCOUS. Step 2 keeps it
  only as `*_c2d_raw.nmf`. The real one comes from the merge step, as in
  xcut step 3/5.
- **Construct2D asks y/n** before accepting a C-grid for an open curve (an
  open curve looks like a blunt TE). Step 2 answers `y`. All other settings
  go through `grid_options.in`, with `topo='CGRD'` set explicitly.
- **Tangent landings leave long slivers, not just thin corners.** In the
  middle half of each upper bridge, the gap to the wall it lands on is only
  0.0017 (slat → main) and 0.0004 (main → flap). Step 1 reports this as a
  warning. It is the main input to the landing trade-off before step 4.

## Goal

Mesh the 30P-30N (slat + main + flap) with Construct2D plus block-structured
fill meshes, merge everything into one `p3d`/`nmf`, and take it through
`Construct2D_to_SOD2D` to a high-order SOD2D mesh. This is the same route
`xcut_pipeline/` takes for the blunt-TE OAT15.

Construct2D extrudes only one smooth curve. Here that curve is the **outer
envelope** of the whole configuration. Everything inside the slots is left out
of it and filled with separate block-structured meshes.

## Step 0 (prerequisite): the Construct2D build

The whole plan depends on the **XCUT variant** of Construct2D (NWKE=0). That
variant exists only as **uncommitted** edits in
`nek-wing/foil_mesh/Construct2D/src/` (`edge_grid.f90`, `menu.f90`,
`surface_grid.f90`, `hyperbolic_surface_grid.f90`, `util.f90`, `vardef.f90`).
The copy in this repo (`AeroHexMesh/Construct2D/src/`) is the stock code.

The stock code does two things that break this pipeline:

- On `GRID`, it calls `transform_airfoil`: it subtracts `min(x)` and divides
  x and y by `0.5*(x(1)+x(n))`. Verified: the OAT15 XCUT curve
  (`oat15_slc_xcut_spline_v2.dat`, ends at x=20) built with the stock source
  comes out at exactly 1/20 scale, and its first and last points are moved to
  `(RADI+0.5, 0)`.
- With BUFF it refuses an open curve ("Buffer airfoil must be closed").

The XCUT variant skips the rescale when NWKE=0, accepts an open curve whose
i=1 and i=imax stay apart (a real TE gap), and marches i=1 and i=imax as cut
columns.

**Tasks before step 1:**

1. Port the XCUT edits into `AeroHexMesh/Construct2D/src/` and commit them. Do
   not port the `.bak`/`.o`/`.mod` files.
2. While doing so, change the p3d coordinate format in `write_srf_grid` from
   `es17.8` to `es25.16`. Today the mesh keeps only 8 significant digits, so a
   surface point read back from the p3d differs from the source by about
   1e-9 relative. Without this change, "bit-identical surface points" cannot
   hold past step 2 (see constraint C5).
3. Run step 2 non-interactively by piping the menu commands into
   `construct2d` (stdin works: `printf 'SOPT\nNWKE\n0\nQUIT\n...GRID\nBUFF\nQUIT\n' | construct2d curve.dat`),
   and save the exact command string in the meta file.

## The envelope (the single Construct2D curve)

Index order follows the XCUT convention (as in `oat15_*_xcut*`): **i=1 is the
far end of the lower downstream line**, the curve runs along the lower side to
the front, around, and back along the upper side, and **i=imax is the far end
of the upper downstream line**.

- **Lower side (i increasing):** lower downstream line → flap lower TE corner →
  flap lower surface → lower bridge → main lower surface from the cove's
  bottom corner forward → lower bridge → slat hook tip → slat bottom → slat
  leading edge.
- **Upper side:** slat leading edge → slat outer surface →
  slat trailing tip → upper bridge → main upper surface → main trailing edge →
  upper bridge → flap upper surface → flap upper TE corner → upper downstream
  line.
- **Downstream end: left open.** The two downstream lines end apart. They
  start at the flap TE corners (gap 0.00588) and are not joined. The XCUT variant accepts this, and it is
  the same layout as the OAT15 TE box. Earlier drafts closed the far ends
  with a short segment; the open end is simpler and needs no trim.

Real surfaces that are **not** on the envelope (they are walls of the blocks):

| Wall | Belongs to block |
|---|---|
| slat cove, hook tip → trailing tip | slat–main |
| main nose, between the two landing points | slat–main |
| main cove, shelf + vertical wall | main–flap |
| flap nose, between the two landing points | main–flap |

### Decision D1 (decided 2026-10-04): no upstream channel

The sketch routed two lines from the slat bottom far upstream, which kept the
slat bottom inside a block. That channel is **dropped**. The slat bottom is on
the envelope, and the lower bridge goes straight from the main element to the
hook tip.

Why: the C-grid has only one cut (at the downstream end), so the channel's far
end would have been a cap in the middle of the curve. That cap cannot be
trimmed, and the channel's dog-leg shape put two concave pockets on the
envelope, where hyperbolic marching (HYPR) folds. The hook tip does not need
the channel. Measured with a straight bridge from main index 103:

- the envelope turns **23.7° convex** at the hook tip (0° if the bridge
  arrives tangent to the slat bottom);
- the slat–main block's corner at the hook tip is **131.5°** (155° with a
  tangent arrival). The slat's own 24.8° hook wedge is solid, not block.

The sketch's green upstream lines no longer apply. Its downstream lines and its
blue/red slot outlines still do.

## The block-structured regions

1. **Slat–main:** slat cove, upper bridge, main nose, lower bridge.
2. **Main–flap:** main cove, lower bridge, flap nose, upper bridge.
3. **Downstream channel:** between the two lines off the flap TE (TE-box style,
   one block, a VISCOUS face at the flap TE, FARFIELD at the far end).

Each bridge and each line is an interface shared between the Construct2D mesh
and a block mesh, so the points along it must match one-for-one.

**The block layout has to be sketched before step 1**, not after step 3. The
point counts on every bridge and line are fixed in step 1, and they follow
from the block layout (constraint C3). A coarse first pass of steps 1–3 can
use provisional counts, but the counts must be recomputed before any block is
built.

- **Slat–main:** the cove (75 intervals, hook tip → trailing tip) and the main
  nose (55 intervals, index 103 → 158) face each other. A single H-block needs
  equal counts, so the nose needs ≥ 20 inserted points (rule C2). The two
  bridges are 0.094 (upper) and 0.115 (lower) long. Block corners: 153° at the
  slat trailing tip and 131.5° at the hook tip (both fine), **1.3° at the
  upper landing** (main 158) and 53° at the lower landing (main 103), measured
  with straight bridges (see "Landing points"). The cove
  shear layer runs from the hook tip toward the slat TE, and it is the main
  target for 30P-30N slat-noise LES. Lay out the block lines along it.
- **Main–flap:** six sides. Main TE → shelf (40) → cove top corner (90°) →
  vertical wall (17) → cove bottom corner → lower bridge (0.200 long) → flap
  nose (index 66 → 139, 73 intervals) → upper bridge (0.046 long) → main TE.
  Block corners: 165° at the main TE, **0° at the upper landing** (flap 139),
  91.5° at the cove bottom corner, 19.5° at the lower landing (flap 66).
  The two bridges differ in length by 4.4×, so one H-block with matched
  bridge counts would have a 4.4× spacing jump across it. Plan for three
  sub-blocks, for example: a cove-box under the shelf, a slot block between
  the cove bottom corner and the flap nose, and a gap block under the main TE.
  Fix the internal interface counts on paper first.

## Pipeline steps (new folder, `xcut_pipeline/` left untouched)

0. **Port and commit the XCUT Construct2D build** (see step 0 above).
1. **Build the envelope curve** and a `.meta.json` recording every index range,
   the landing indices, the source-file index ranges of the block walls, and
   the spline parameters of every surface (rule C2).
2. **Construct2D:** HYPR, CGRD, NWKE=0, BUFF, scripted (step 0, task 3). The
   run parameters (YPLS, RECD, CFRC, JMAX, the HYPR smoothing settings,
   XDAMP) go into the meta file.
3. **Check, don't trim.** The downstream end is open and there is no
   upstream cap, so there is nothing degenerate to trim. Check the two cut
   columns (i=1, imax) and the whole block against the quality gates (C8). Trim only if those
   columns fail, as in OAT15's `--nelm-clean`.
4. **Build the blocks.** Take every interface point from the step-2 p3d, never
   from the source files (C5).
5. **Merge** all blocks into one `p3d`/`nmf`. The xcut step5 handles exactly
   2 blocks; this needs an N-block merge that checks every ONE_TO_ONE by
   coordinates.
6. **SOD2D conversion** (needs work, see C9).

Start light: get steps 0–3 working on a coarse envelope. Lay out the blocks on
paper before designing any block mesh.

## Geometry facts (measured from the raw files)

Files: `sample_airfoils/30p30n_source/30P-30N-{Main,Slat,Flap}-raw.dat`, plain
`x y` per line. Indices below are 0-based, after dropping the duplicate closing
point of main and slat.

| Element | Points | Feature | Index | (x, y) |
|---|---|---|---|---|
| Slat (closed) | 200 | trailing tip (sharp, 171.8°) | 0 | (0.0188, 0.0051) |
| | | hook tip (sharp, wedge 24.8°) | 75 | (-0.0267, -0.1087) |
| | | leading edge (rounded) | 143 | (-0.0854, -0.0996) |
| Main (closed) | 220 | trailing edge (sharp, 172.2°) | 0 | (0.8740, 0.0309) |
| | | cove top corner (90°) | 40 | (0.6999, 0.0320) |
| | | cove bottom corner (79°) | 57 | (0.6999, -0.0172) |
| | | leading edge | 122 | (0.0438, -0.0167) |
| Flap (open TE) | 242 | lower TE corner | 0 | (1.1283, -0.1458) |
| | | leading edge | 97 | (0.8715, 0.0058) |
| | | upper TE corner | 241 | (1.1309, -0.1405) |

Index directions:

- **Slat:** 0→75 is the cove, 75→143 the bottom, 143→199→0 the outer surface.
- **Main:** 0→40 is the cove shelf, 40→57 the cove wall, 57→122 the lower
  surface (aft to nose), 122→219→0 the upper surface (nose to trailing edge).
- **Flap:** 0→97 is the lower surface (TE to nose), 97→241 the upper surface.
  The TE gap is 0.00588.

Raw surface resolution (spacing between neighbouring points):

| Element | min | max | max neighbour ratio |
|---|---|---|---|
| Slat | 5.0e-4 | 7.0e-3 | 1.13 |
| Main | 1.1e-3 | **4.9e-2** | 1.16 |
| Flap | 4.3e-4 | 1.1e-2 | 1.08 |

The raw points are smooth, but far too coarse for a wall-resolved LES surface
(the main element has 5% chord spacing mid-chord). BUFF meshes exactly the
points it is given, so all surface refinement has to happen in step 1 (C2).

## Landing points

Angles below are for straight bridges. The "block corner" is the angle inside
the slot block, between the bridge and the wall that continues into the
block.

| Bridge | From → to | Envelope turn at the landing | Block corner at the landing | Envelope turn at the departure | Block corner at the departure |
|---|---|---|---|---|---|
| slat upper | slat tip → main 158, (0.1077, 0.0353) | 0.3° (tangent) | **1.3°** | 18.5° convex | 153° |
| flap upper | main TE → flap 139, (0.9177, 0.0172) | 0.9° (tangent) | **0°** | 6.9° convex | 165° |
| slat lower | main 103, (0.0591, -0.0315) → hook tip | 53° concave | 53° | 23.7° convex | 131.5° |
| flap lower | main 57 (cove bottom) → flap 66, (0.9002, -0.0224) | — | 19.5° | — | 91.5° |

- **Departures are not a problem.** Their block corners are wide, and the
  envelope kinks there are convex, which hyperbolic marching handles. Making
  a departure tangent (a Hermite bridge leaving along the upstream element's
  surface) would remove the kink and widen the block corner to about 172°.
  That is optional, and worth trying if the extruded cells at the slat tip
  look skewed.
- **Landings are the real trade-off.** A tangent landing gives a smooth
  envelope but a near-zero block corner (1.3°, 0°). A steeper landing helps
  the block but puts a concave kink in the envelope. The two lower bridges
  show this the other way round: a straight bridge gives the block 53° and
  19.5° but a concave envelope turn (53° at main 103). An S-blend that lands
  tangent removes the concave turn but takes those block corners to near zero
  too.
- **Starting choice:** favor the envelope, since it drives the whole outer
  mesh. Use tangent upper landings and an S-blend on the lower bridges, then
  handle the four thin corners inside the blocks. If step 4 cannot get those
  corners through the quality gates (C8), move each landing a few degrees off
  tangent. Log the trade-off in the meta file (C6).
- Lower landing points sit just aft of each nose, where surface curvature has
  dropped off.

## Requirements for the downstream lines

- Smooth, near-parallel, ending in a straight far section (user requirement).
- **Downstream pair:** attaches at the two flap TE corners. Build it as one
  centerline from the flap TE midpoint along the TE bisector, a Hermite turn
  (`--turn-length`) onto a straight far section, and straight on to
  `--xfar`. Both lines share the centerline's point distribution exactly.
  `--far-end free` (default): the far section runs at `--far-angle`, and the
  lines are exact ±(TE gap)/2 translates of the centerline. `--far-end
  mirror`: the far section runs along y=0, and the offset rotates to (0, ±h),
  so the far ends are mirror images about y=0 (for builds without C10's
  fix). Ends open.

## Constraints step 1 must respect so the later steps work

- **C1 Matching counts:** the two bridges of each slot need the same number of
  points if a single H-block spans them, and so do the two lines of a channel.
  The real counts come from the block layout, so sketch it first.
- **C2 Surface refinement by insertion only.** Fit one spline per element
  through the raw points, and **insert** new points on it. Never move or
  replace a raw point. Raw points then stay bit-identical (they are a subset),
  and the counts can still be raised to any target. The same spline must be
  used for the envelope surfaces, the block walls, and the SOD2D high-order
  wall projection (C9). If the converter projects high-order nodes onto a
  different curve than the one the linear mesh was built on, the wall moves.
  Store the spline knots and control points in the meta file.
- **C3 Spacing continuity along the curve:** at every junction (surface ↔
  bridge, bridge ↔ surface, flap corner ↔ downstream line) the neighbouring
  spacings differ by ≤ 1.2×. The bridge and line distributions are chosen
  from their end spacings (geometric/tanh), not from a fixed count.
- **C4 Normal spacing across the interfaces.** Construct2D uses one first-cell
  height (from YPLS/RECD/CFRC) for the whole j=1 row, so bridges and lines get
  wall-type clustering too. Every block cell touching a bridge or line must
  start at that same height and grow at the same rate, so there is no size
  jump across the interface. Record the height in the meta file.
- **C5 Interface points come from the mesh.** Every node shared between the
  Construct2D block and a block mesh is read from the step-2 p3d. That
  includes bridge and line points, and the endpoints where block walls meet
  the envelope: slat tip, hook tip, slat LE, main TE, cove bottom corner,
  landing points and flap corners. Do not recompute these from the source.
  With the es25.16 output (step 0) they agree with the source to round-off.
- **C6 Corners are measured and logged.** Step 1 prints the envelope's turn
  angle at every junction and each block's interior angle at every corner, and
  refuses to write if a concave envelope turn exceeds a set limit (first trial:
  30° within one local spacing).
- **C7 Self-check:** refuse to write a curve that crosses itself or cuts
  through any wall left for the blocks. Also refuse if the curve passes closer
  to any block wall than about 0.2 × the local slot gap.
- **C8 Quality gates** for every block (Construct2D and the blocks): no
  negative Jacobians, minimum cell angle ≥ 20° away from the wedge corners,
  growth ratio ≤ 1.2 in both directions, and cell-size ratio across every
  interface ≤ 1.3. Put the numbers in the stats output and fail the step when
  a gate fails.
- **New files only:** never overwrite existing outputs.

## C9 Downstream converter gap (`Construct2D_to_SOD2D`)

`linear_mesh/wall_spline.py` and `mesh_extrusion.py` accept only a 1-block
C/O-grid or the 2-block OAT15 shape (`unsupported multiblock configuration`).
This mesh has three separate walls, on both closed elements and the open flap,
spread across about 5–7 blocks. The converter needs:

- general N-block NMF reading, and ONE_TO_ONE node aliasing for every pair of
  blocks;
- VISCOUS ranges grouped per element, then assembled in order into one wall
  curve per element (slat and main closed, flap open with its TE face);
- high-order wall projection onto **that element's** spline (C2), not one
  global wall spline.

Plan this alongside step 4. Without it the merged mesh cannot be converted.

## Where things are

This plan lives in `AeroHexMesh/Construct2D/envelope_pipeline/`. In the same repo:

- `sample_airfoils/30p30n_source/`: the three raw element files.
- `visualize_curve.py`: views any curve `.dat` with its `.meta.json`; works
  for the envelope as-is.
- `xcut_pipeline/xcut_common.py`: reusable helpers (centerline spline + offset,
  Hermite points, geometric growth, p3d/nmf and meta I/O, quality stats).
- `../Construct2D_to_SOD2D/linear_mesh/`: the converter that C9 extends.

Outside this repo: `nek-wing/foil_mesh/Construct2D/src/` holds the XCUT
Construct2D edits (step 0), and also the two rejected earlier attempts:
`xcut_pipeline/step1_build_multielement_curve.py` (branch-cut topology) and
`xcut_pipeline/step1_build_slotted_curve.py` with
`sample_airfoils/30p30n_slotted.*` (curve threaded through the slots). Both put
the cove and nose surfaces on the Construct2D curve. Their small helpers
(Hermite bridge, arc indexing) can be copied if useful; the topology cannot.

## Review findings (2026-10-04, evidence for the changes above)

- Stock `AeroHexMesh/Construct2D/src` built and run with BUFF/NWKE=0 on
  `oat15_slc_xcut_spline_v2.dat`: wall row = input / 20 (max error 4e-10);
  first and last points moved to (15.5, 0).
- The OAT15 meshes in `nek-wing/foil_mesh/Construct2D/`
  (`oat15_xcut_full.p3d`, `oat15_slc_xcut_spline.p3d`) are at true scale (wall
  x 0…20 and 0…40). The outer boundary is not a circle, so they were made with
  HYPR. Conclusion: they came from the patched build, which is not in this repo.
- `create_farfield_xcut` (the patched elliptic far-field) centres its arc on
  `min(x)` of the curve, i.e. the slat leading-edge region once the channel is
  dropped. HYPR does not use it either way.
- Angles and lengths quoted above were computed directly from the raw files.
