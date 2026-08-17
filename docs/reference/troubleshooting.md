# Troubleshooting / known caveats (unified reference)

Gotchas and known issues scattered across individual READMEs and code
comments, collected in one place.

## Unverified: inlet/outlet AoA sign convention

`Construct2D_to_SOD2D/linear_mesh/` and `Construct2D_to_NEKRS/linear_mesh/`
both classify farfield boundary points as inlet/outlet via
`dot(outward_direction, (cos(AoA), sin(AoA))) < 0`. This convention
(flow tilts toward `+y` for positive AoA) behaves consistently across
O-grid and C-grid, but has **not** been checked against an
independently-known-good result. See
[`boundary-conditions.md`](boundary-conditions.md#inlet-outlet-classification-shared-math).

## MPI crashes at `num_partitions ≥ 3` on a large mesh

Two *unrelated* crashes, both needing their own workaround, both applied
together in every SOD2D GPU job script in this repo — see
[`cluster-setup.md`](cluster-setup.md#slurm-job-scripts-gpu-steps) for
the exact flags. If you hit a segfault running SOD2D on a large
partitioned mesh, check whether one of `UCX_TLS=...` or `--mca io
^ompio` got dropped before assuming it's a new bug.

## NekRS: two settings `compute_case_params.py` does *not* automate

- **Viscosity/AoA** in `naca.par` — hand-set, must match whatever
  `Re`/`angle_of_attack` the mesh was actually generated with.
- **Inflow angle is hardcoded** as a literal number inside `naca.oudf`'s
  GPU (OKL) kernel — not read from `naca.par`, not derived from
  anything. Changing AoA means editing **two** places by hand.

See [`Construct2D_to_NEKRS/high_order_mesh.md`](../Construct2D_to_NEKRS/high_order_mesh.md#two-settings-this-step-does-not-automate).

## NekRS: `NUMBER_ELEMENTS_Y` means two different things

Same variable name, genuinely different meaning depending which `.usr`
it's patched into — see
[`Construct2D_to_NEKRS/linear_mesh.md`](../Construct2D_to_NEKRS/linear_mesh.md#compute_case_paramspy--deriving-case-parameters-from-the-nmf-directly).
Worth double-checking by hand before trusting a patched case file.

## pyHyp: first-layer height `s0` vs. smallest surface edge

`s0` has to stay well under the input surface's smallest edge, or the
very first marching step self-intersects near whatever feature is that
fine — this is *why* `pyHyp/coarsen_surface.py` exists (raises the
smallest-edge floor by roughly halving point count). If pyHyp's own
extrusion fails near a specific feature (e.g. a wingtip), check whether
that feature's own local edge length is disproportionately small
relative to `s0` first, before assuming it's a marching-parameter issue.
Symmetric caveat: coarsening a small, high-curvature block (e.g. a
wingtip cap) too aggressively can itself introduce near-degenerate
elements in the marched/order-elevated mesh even when the *linear*
mesh's own quality metric doesn't flag it — use `--protect-blocks` to
exempt specific blocks.

## `MeshElasticitySolver.json`: wrong `bc_type` on the wall silently does nothing

The wall's `bouCodes` entry must be `bc_type_non_slip_adiabatic_moving`
(not plain `bc_type_non_slip_adiabatic`) for the elasticity solver's own
Dirichlet routine to impose the correction at all — see
[`boundary-conditions.md`](boundary-conditions.md#bc_type-strings-sod2d-json-configs).
Using the wrong string doesn't error; it just silently leaves the wall
un-corrected (held at zero displacement, same as every other
unrecognized code).

## `etaGrid_to_SOD2D`: `examples/airfoil.msh` has known degenerate quads

~40 near-zero-area quads from the quad-subdivision step that removed the
original file's stray triangles (traced, not yet repaired) — see
[`etaGrid_to_SOD2D/README.md`'s own Example input section](../../etaGrid_to_SOD2D/README.md#example-input).
Prefer `make_demo_cross_section.py`'s clean synthetic cross-section for
anything where mesh quality matters (e.g. isolating a sweep-algorithm
bug from a real-mesh quality issue).

## Login node kills background jobs silently

See [`cluster-setup.md`](cluster-setup.md#login-node-resource-limits) —
always submit real work via SLURM, use `python3 -u` when backgrounding.

## Raw-mesh quality: always check `minSICN` before order-elevating

Standing process rule in `etaGrid_to_SOD2D/` (and good practice
generally): order elevation and downstream conversion/partitioning are
comparatively expensive and hide the *raw* mesh's own quality problems —
check `gmsh.model.mesh.getElementQualities(tags, "minSICN")` on the raw,
pre-elevation mesh first (`etaGrid_to_SOD2D/check_raw_quality.py`).
`minSICN` catches local element distortion/twisting that a simple
corner-based signed-volume check misses.
