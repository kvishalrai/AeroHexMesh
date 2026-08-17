# `Construct2D/` — math and implementation

Deep-dive companion to [`Construct2D/README.md`](../../Construct2D/README.md).
Construct2D's own 2D grid-generation algorithm (elliptic or hyperbolic
PDE-based structured-grid generation from an airfoil curve) is vendored
upstream Fortran, out of scope to re-derive here — see its own
`doc/user_manual.pdf`. This page covers what *this repo's* own pieces do:
driving Construct2D non-interactively, and the two auxiliary Python
scripts.

## Non-interactive (batch) driving

Construct2D is normally an interactive menu program, but it auto-loads
settings from a fixed-name `grid_options.in` file in the working
directory if one is present (`src/menu.f90:set_defaults`) — covering
every surface/volume/output option non-interactively. The only menu
interaction left is driving the main loop itself:
`GRID → SMTH (generate the smoothed grid) → QUIT`, piped via stdin. Each
`generate_*.sh` script in this repo writes its own `grid_options.in`
(a Fortran namelist — `&SOPT`/`&VOPT`/`&OOPT` blocks) then runs:

```bash
printf 'GRID\nSMTH\nQUIT\n' | ./construct2d <airfoil>.dat
```

### Namelist fields that matter for the rest of the pipeline

| Field | Meaning |
|---|---|
| `topo` | `'OGRD'` or `'CGRD'` — the grid topology; downstream tooling (`Construct2D_to_SOD2D`, `Construct2D_to_NEKRS`) branches on this via `mesh_type`, must match. |
| `slvr` | `'HYPR'` (hyperbolic) or `'ELIP'` (elliptic) — which PDE-based grid generator Construct2D itself runs. |
| `nsrf` | Number of points around the airfoil surface — directly sets the O-grid's `imax` (`topo='OGRD'` ⇒ `imax = nsrf`). |
| `jmax` | Number of points in the wall-normal direction, wall to far boundary. |
| `radi` | Far-boundary distance, in chords. |
| `ypls` | Target `y+` for the first wall-normal cell (drives the near-wall spacing). |
| `nwke` | Number of points in the wake region (C-grid only). |

The output `.p3d` (grid coordinates) + `.nmf` (Neutral Map File —
boundary meaning) pair is exactly what every downstream pipeline in this
repo consumes as its own starting point (see e.g.
[`Construct2D_to_SOD2D/linear_mesh.md` §1](../Construct2D_to_SOD2D/linear_mesh.md#1-inputs-plot3d-grid--neutral-map-file)).

## `generate_grid_previews.py`

Regenerates the 4 preview PNGs shown in `Construct2D/README.md`, via
`postpycess.py`'s own `read_grid()`/`plot_grid()` (matplotlib, `Agg`
backend — no display needed, safe to run on a login node/CI).

## `postpycess.py`

A vendored, Python-3-ported copy of Construct2D's own bundled CFD
postprocessor (GPLv3, © 2013–2018 Daniel Prosser) — reads Plot3D
grid/function files and can plot contours or airfoil-surface data
interactively. Used here purely as a plotting library
(`generate_grid_previews.py`'s own dependency); its own CFD-postprocessing
features (contour plots of a flow solution) aren't exercised by anything
else in this repo.
