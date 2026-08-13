# Applying the wall correction, and running NekRS

`linear_mesh/`'s `p3d_to_gmsh_nek.py` + `gmsh2nek` place the mesh's
high-order boundary nodes on the wall by straight-line interpolation
between the linear Construct2D corner points — the curved wall boundary
is faceted, not truly smooth. `linear_mesh/smooth_2D/` already fit a cubic
spline through the wall's corner points and resampled it at the mesh's
actual node positions **on the 2D mesh only**, storing the result
compactly (one entry per element, for every element in the 2D plane) in
`geom_3x3_naca_gen_spline_info.dat`.

This directory (`example_nekrs/`) is a NekRS case whose `usrdat2` calls
`get_smooth_data`, which reads that same file and stamps it onto every
element of the real 3D mesh — repeating across every spanwise level,
since the mesh is just a periodic Z-extrusion of one 2D pattern — before
running the actual flow solve. No new spline fitting happens here; this
step is pure propagation of the correction `smooth_2D/` already computed.

## Where the input files come from

Only the hand-written `.usr`/`.par`/`.udf`/`.oudf`/`airfoil0.sh` files are
tracked here — the mesh and geometry-correction data are produced
upstream and copied in (see `../linear_mesh/README.md`'s
[Running it end to end](../linear_mesh/README.md#running-it-end-to-end)):

| File (gitignored, not tracked) | Produced by |
|---|---|
| `naca.re2` | `linear_mesh/arglist.sh`'s 3D extrusion — write directly here via its `output_re2` argument |
| `geom_3x3_naca_gen_spline_info.dat` | `linear_mesh/smooth_2D/`'s `nek5000` run — copy in by hand |

Everything else here (`naca.usr`, `naca.par`, `naca.udf`, `naca.oudf`,
`ci.inc`, `ci.oudf`, `intercept.oudf`, `my_drag.oudf`, `naca.nek5000`,
`nekrs.upd`, `airfoil0.sh`) is hand-written case configuration and is
tracked normally.

Files this folder should **not** keep tracked: anything NekRS itself
writes when you run it (`.cache/`, `base.fld`, `naca0.f*` restart/output
fields, `newrea.out`, `cmake.log`, SLURM `out.o`/`error.e`) — those are
gitignored; rerun the case to regenerate them.

## Before running

`linear_mesh/compute_case_params.py` (via `arglist.sh`) already keeps
`naca.usr`'s `NUMBER_ELEMENTS_X/Y/Z` `#define`s and `naca.par`'s matching
`[CASEDATA] number_elements_x/y/z` in sync with the actual mesh
automatically. Two things it does **not** touch, and that are easy to
leave inconsistent by hand:

- `naca.par`'s `[VELOCITY] viscosity` (set to `-Re`) and
  `[CASEDATA] p_aoa`/`p_u0` — the actual physics you want to run, not
  mesh-derived.
- `naca.oudf`'s `codedFixedValueVelocity` — it's a hardcoded
  `aoa`/`cos`/`sin`, **not** wired to `naca.par`'s `p_aoa`
  (`kernelInfo["defines/p_*"]` set in `UDF_LoadKernels` aren't visible in
  a coded-BC kernel's own compile context, a NekRS quirk, not a bug to
  fix here) — keep both in sync by hand if you change the angle of
  attack.

## Running it

`sbatch airfoil0.sh` from this folder — see the script for module/queue
setup (site-specific; it currently assumes an existing `NEKRS_HOME`
install and BSC MareNostrum 5's module system).

## Verifying the result

- Build/mesh-load log lines to check for: `building nekInterface for
  lx1=..., lelt=N and lelg=N` should show `N` matching the mesh's actual
  element count (`arglist.sh`'s "elements written" line from the 3D
  extrusion step), and `loading nek ... done` means `usrdat2`/
  `get_smooth_data` ran without error.
- Once timestepping starts, sanity-check `t-cd-cl` (drag/lift) and the
  velocity/pressure min-max lines NekRS prints every step: `Cd`/`Cl`
  should be O(1) and physically plausible for the case, not blown up;
  `dt` should stay roughly stable (a collapsing `dt` alongside huge
  velocity/pressure values is solution divergence, not a mesh/BC problem
  — checked separately from the items above).
