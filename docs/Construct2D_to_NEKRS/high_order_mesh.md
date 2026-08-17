# `Construct2D_to_NEKRS/high_order_mesh/` — math and implementation

Deep-dive companion to
[`Construct2D_to_NEKRS/high_order_mesh/README.md`](../../Construct2D_to_NEKRS/high_order_mesh/README.md).
This step does **no new smoothing math of its own** — it's the
mechanism that takes `linear_mesh/`'s already-corrected 2D geometry
(§2.2–2.3 of
[`linear_mesh.md`](linear_mesh.md#step-2-smooth_2dnaca_gen_spline_infousr--fischers-wall-smoothing-algorithm))
and copies it onto the real 3D mesh, at whatever polynomial order the
production case actually uses.

## Why a separate interpolation step at all

`linear_mesh/smooth_2D/`'s case is deliberately run at a **fixed, small
order** (`lx1 = 4`, enforced by an explicit abort check in
`smooth_geom`: `if (lx1.ne.4) call exitti(...)`) — the spline correction
itself doesn't depend on solver order, so there's no need to redo the
(relatively expensive, iterative) `laplaceh` Helmholtz solve at
production order/mesh size. Instead, the corrected 2D geometry is sampled
once, at a fixed 4×4-per-element GLL grid, and saved to
`geom_3x3_naca_gen_spline_info.dat` — a per-element table of `(x, y)`
pairs. The real production case then just **interpolates** this fixed
reference sampling onto its own (possibly different-order) elements.

## `get_smooth_data` (`example_nekrs/naca.usr`)

Called from `usrdat2` (the standard NekRS hook for modifying mesh
coordinates before the run starts) and again from `userchk` on the first
timestep (`istep.eq.0`, guarding against a restart re-reading a mesh that
was already corrected once):

1. Rank 0 reads the whole `geom_3x3_naca_gen_spline_info.dat` table
   (`4*4` `(x,y)` pairs per 2D element) and broadcasts it to every rank.
2. For each of the mesh's own elements `e` (real, 3D — `nelv` of them),
   find its corresponding 2D reference element via
   `e2d = mod1(lglel(e), lelxy)` — i.e. wrap the global element number
   modulo the 2D section's own element count `lelxy = NUMBER_ELEMENTS_X ×
   NUMBER_ELEMENTS_Y`, since every spanwise plane repeats the same 2D
   element layout.
3. `map_m_to_n(xm1(...), lx1, xy4(...), 4, ...)` — a Lagrange-basis
   change-of-interpolant-order operation (NekRS/Nek5000's own standard
   utility for resampling a field from an `m`-point to an `n`-point GLL
   grid) — interpolates the fixed 4-point reference `(x, y)` onto the
   element's own `lx1`-point grid, overwriting `xm1`/`ym1` directly (the
   mesh's own coordinate arrays) for the **first** spanwise (`iz=1`)
   layer.
4. **Copy, don't recompute, the other spanwise layers**:
   `xm1(:,:,iz,e) = xm1(:,:,1,e)` for every `iz = 2..lz1`. Because this
   pipeline extrudes spanwise by pure translation (§2,
   [`linear_mesh.md`](linear_mesh.md)), the corrected 2D `(x,y)` shape is
   identical at every spanwise plane — no interpolation or recomputation
   needed there, just a copy.
5. `geom_reset(1)` recomputes the mesh's geometric factors (Jacobians,
   etc.) from the now-directly-overwritten coordinate arrays — required
   any time `xm1`/`ym1`/`zm1` are modified in place like this.

## Two settings this step does *not* automate

Flagged explicitly in the README, worth repeating here since they're
easy to silently desync from the actual case physics:

- **Viscosity/AoA** live in `naca.par`, hand-set to match whatever `Re`/
  `angle_of_attack` the mesh was actually generated with —
  `compute_case_params.py` never touches these (it only derives
  element-count parameters from the NMF, §
  [`linear_mesh.md`](linear_mesh.md#compute_case_paramspy--deriving-case-parameters-from-the-nmf-directly)).
- **The inflow angle is hardcoded** as a literal number inside
  `naca.oudf`'s GPU (OKL) kernel — not read from `naca.par`, not derived
  from anything — so changing AoA means editing two places by hand, not
  one.
