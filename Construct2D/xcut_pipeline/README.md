# Blunt-TE XCUT meshing pipeline (zone (a))

Turns a blunt-trailing-edge airfoil (e.g. `sample_airfoils/oat15.dat` or
`oat15_slc.dat`) into a 2-block Construct2D-style mesh: a CGRD/NWKE=0/BUFF
main mesh with prescribed straight wake lines, plus a thin gap-filler
block that closes the physical TE face. Every step writes new files —
nothing here overwrites an existing file in place.

Run from anywhere with `python3 <path-to-script> ...` (each script finds
`xcut_common.py` next to itself automatically).

## Workflow

**Step 1 — build the curve** (`step1_build_curve.py`)

```
python3 step1_build_curve.py \
    --airfoil sample_airfoils/oat15_slc.dat \
    --lfar 20.0 --nwake 50 \
    --out sample_airfoils/oat15_xcut.dat
```

Reads your airfoil geometry, auto-strips any pre-existing TE-closure
artifact the file may already carry (both `oat15.dat` and
`oat15_slc.dat` in `sample_airfoils/` close themselves with a shared
point at the TE — this is detected and removed automatically; you'll
see it reported), then appends wake lines (geometric growth, matching
the airfoil's own local TE spacing as the first step) from x=1 out to
`--lfar`, closing the whole curve with a shared point at `(--lfar,
0.0)` so Construct2D sees it as closed (required for BUFF).

**Hard requirement: the airfoil surface coordinates are never
modified — not even the two TE corner points.** Every point between
(and including) the two TE corners is written back bit-for-bit
identical to the source file (verified with a round-trip test; the
write format uses full float64 precision, `%.17g`, specifically so
this holds even for source files with more decimal digits than a
naive fixed-width format would preserve). Equidistant/parallel wake
lines are achieved without touching the corners: both wake lines are
built relative to the TRUE corner midpoint `y_mid = 0.5*(y_bot+y_top)`
(generally not exactly 0) and the TRUE half-gap `h =
0.5*(y_top-y_bot)`, so they land exactly on the corners' own
unmodified y-values while still staying exactly parallel/equidistant
to each other — no nudging, no approximation.

- `--lfar` is an **absolute x-coordinate** (chord normalized to 1.0) —
  e.g. `--lfar 20.0` means the wake extends to x=20, i.e. 19 chords
  beyond the TE. If you meant a wake *length* instead, pass `1.0 +
  that length`.
- `--nwake` is points **per side** (bottom and top each get this many,
  including the TE corner as the first one).
- `--wake-type straight` (default) — constant-y lines, exactly as
  validated in every run so far.
- `--wake-type spline` — a single smooth cubic-Hermite **centerline**
  from the TE-bisector direction (averaged from the actual local
  upper/lower surface tangents at the TE) to horizontal far
  downstream, then offset by ±h for the two sides. Because both sides
  are the same curve shifted by a constant, they come out **exactly**
  parallel and equidistant (gap = 2h everywhere) and with identical
  point distributions, by construction — not just approximately.
  `--pull-length` (default 1.0, in chords) controls how far the TE
  angle's influence reaches before the curve flattens out; it's an
  absolute length, not a fraction of the wake length, since both
  endpoints sit at y=0 and a distance-scaled tangent would badly
  overshoot on a long wake. Console output reports the TE-bisector
  angle and the peak centerline offset so you can sanity-check the
  shape before running Construct2D.

Writes `<out>` and `<out-without-ext>.meta.json` (records exactly
which curve indices will be the TE corners once meshed — later steps
read this instead of re-detecting anything).

**Step 2 — generate the mesh (you, interactively)**

Load the step-1 `.dat` in Construct2D and run **CGRD + NWKE=0 + BUFF**
as usual. This writes the raw `.p3d`/`.nmf`.

**Step 3 — trim the wake-tip closure** (`step3_trim_mesh.py`)

```
python3 step3_trim_mesh.py \
    --p3d oat15_xcut.p3d --meta sample_airfoils/oat15_xcut.meta.json \
    --stats-p3d oat15_xcut_stats.p3d \
    --nelm-clean 1 \
    --out-prefix oat15_xcut_trimmed
```

The artificial closure point from step 1 leaves 1 (or more) degenerate
sliver elements next to i=1 and i=imax. `--nelm-clean N` removes N
columns from each end. N=1 is what's been needed so far for
Lfar=20/nwake=50 on OAT15 — inspect the raw mesh near the wake tips if
you change those parameters substantially and aren't sure.

`--stats-p3d` is optional: pass Construct2D's `*_stats.p3d` mesh-quality
file (skew angle, xi-growth, eta-growth) from the same step-2 run, and
it's trimmed the same N columns off each end so you can still inspect
quality on the mesh you'll actually use, written to
`<out-prefix>_stats.p3d`.

Writes `<out-prefix>.p3d`, `.nmf` (single-block, same BC style as any
other Construct2D CGRD output), `<out-prefix>_stats.p3d` (if
`--stats-p3d` was given), and `.meta.json` (carries forward the
now-shifted bottomTE/topTE column indices).

**Step 4 — build the TE-face gap-filler box** (`step4_build_te_box.py`)

```
python3 step4_build_te_box.py \
    --trimmed-p3d oat15_xcut_trimmed.p3d \
    --meta oat15_xcut_trimmed.meta.json \
    --out oat15_te_box.p3d
```

Reuses the trimmed mesh's own bottom/top wake-arm wall rows as the
box's two long edges (already fold-free, already shares nodes with
the main mesh). The number of points spanning the gap is chosen
automatically so spacing is **uniform** and matches the main mesh's
own first wall-normal (j) grid spacing near the TE — not a fixed
count. Refuses to write anything if the fold check finds any
degenerate cells.

Also **always** computes and writes `<out-without-ext>_stats.p3d` —
skew angle, xi-growth, eta-growth for the box block, in the same
format as Construct2D's own `*_stats.p3d`. The box never goes through
Construct2D (it's built directly here), so it never gets a quality
file from Construct2D the way the main mesh does; this fills that gap
with a Python port of Construct2D's own `compute_quality_stats`
(verified against real Construct2D output — see
`xcut_common.compute_quality_stats`'s docstring).

**Step 5 — merge into the final 2-block mesh** (`step5_merge_blocks.py`)

```
python3 step5_merge_blocks.py \
    --trimmed-p3d oat15_xcut_trimmed.p3d \
    --box-p3d oat15_te_box.p3d \
    --meta oat15_te_box.meta.json \
    --out-prefix oat15_final
```

Coordinate-verifies both `ONE_TO_ONE` interfaces before writing
anything (refuses on mismatch), then writes the final `.p3d`
(2-block) and `.nmf`:

| Face | Block | Type | Notes |
|---|---|---|---|
| i=1, i=imax | 1 | FARFIELD | wake tips, x=Lfar |
| j=1, i=1..bottomTE | 1 | ONE_TO_ONE → block 2 | reversed direction |
| j=1, bottomTE..topTE | 1 | VISCOUS | the real airfoil surface |
| j=1, topTE..imax | 1 | ONE_TO_ONE → block 2 | same direction |
| j=jmax | 1 | FARFIELD | outer boundary |
| i=1 | 2 | VISCOUS | the TE face itself |
| i=imax | 2 | FARFIELD | far end of the gap strip, x=Lfar |

If a `_stats.p3d` sits next to **both** `--trimmed-p3d` and
`--box-p3d` (block 1's from Construct2D via step 3's `--stats-p3d`,
block 2's from step 4's own quality computation), this step also
merges them into `<out-prefix>_stats.p3d` — one file covering both
blocks' quality, using the same multi-block convention as the mesh
itself (`xcut_common.write_p3d_function_multiblock`). Silently skipped
(with a note) if either is missing — always optional, never required.

## Viewing the result (`postpycess.py`)

`postpycess.py` (in the Construct2D root) reads both single- and
multi-block `.p3d`/`_stats.p3d` files, so a final 2-block mesh from
this pipeline can be loaded directly by prefix (e.g. `oat15_final`)
like any other Construct2D case. Grid and contour plots draw every
block on the same axes; any block with its own quality data (from a
merged `_stats.p3d` as above) is colored/contoured using that data,
with one shared color scale across all colored blocks.

## Notes / things worth knowing

- `--lfar`/`--nwake` reproduce exactly the OAT15 case worked out
  interactively (Lfar=20, nwake=50, h≈0.002498, growth ratio≈1.167) —
  step 1 was cross-checked bit-for-bit against that earlier, hand-built
  run.
- The reversed-range `ONE_TO_ONE` syntax (`S2 > E2`) mirrors the
  convention already used in this project's own `naca0012_sharp.nmf`.
- Every `ONE_TO_ONE` connection is checked against actual coordinates
  (not just index arithmetic) before any file is written.
