# p3d-to-nekrs

Converts a 2D airfoil mesh (O-grid or C-grid) into a wall-spline-smoothed,
spanwise-periodic 3D mesh for NekRS, via Gmsh and Nek5000's own mesh tools.

This is the **middle step** of a 3-step pipeline:

1. **Construct2D** (separate, external) — generates the 2D Plot3D mesh
   (`.p3d`) and its Neutral Map File (`.nmf`) for a given airfoil and mesh
   topology (O-grid or C-grid). Not included here; run it yourself.
2. **This directory** — converts the 2D mesh to Gmsh format with
   wall/inlet/outlet boundary classification, fits and applies a wall
   spline correction to it, and extrudes it spanwise into a periodic 3D
   mesh.
3. **`../high_order_mesh/example_nekrs/`** — applies that same correction
   to the real 3D mesh and runs NekRS.

## Pipeline stages

```
2D .p3d/.nmf   →  p3d_to_gmsh_nek.py     →  linear .msh (wall/inlet/outlet tagged)
               →  Gmsh (.geo)            →  order-elevated .msh
               →  gmsh2nek               →  2D .re2
               →  smooth_2D/ (Fischer's smooth_geom0)  →  geom_3x3_*.dat (wall
                                                            spline correction)
               →  arglist.sh (re2torea/n2to3/reatore2) →  3D .re2, spanwise-periodic
```

There's no single orchestrating script for all of this (unlike the SOD2D
pipeline's `run_pipeline.py`) — see [Running it end to end](#running-it-end-to-end)
below for the actual command sequence.

## Prerequisites

- Python 3 with `numpy`.
- [Gmsh](https://gmsh.info/), invoked as a CLI.
- MPI, for Nek5000.
- `gmsh2nek`, `genmap`, `re2torea`, `reatore2`, `n2to3` built from the
  `../CFD_code/Nek5000` submodule — see the top-level
  [`Construct2D_to_NEKRS/README.md`](../README.md#cfd_code).

## Computing case parameters

`smooth_2D/naca_gen_spline_info.usr` and `../high_order_mesh/example_nekrs/naca.usr`
(+ `naca.par`) both need to agree on several mesh-derived numbers — how
many elements are in one ring, and for a C-grid, exactly which sub-range
of that ring is the real wall (the rest being the wake-cut seam). Don't
work these out by hand or by inspecting Gmsh output: **`compute_case_params.py`**
reads them directly from the airfoil's `.nmf` file (its `VISCOUS` boundary
entry already encodes the wall's exact element range, for O-grid or
C-grid alike) and can patch them straight into the case files:

```bash
python3 compute_case_params.py <airfoil>.nmf                                  # just print
python3 compute_case_params.py <airfoil>.nmf --smoothing-levels 1 \
    --patch-step1-usr smooth_2D/naca_gen_spline_info.usr                       # step 1
python3 compute_case_params.py <airfoil>.nmf --nz 3 \
    --patch-step2-usr ../high_order_mesh/example_nekrs/naca.usr \
    --patch-step2-par ../high_order_mesh/example_nekrs/naca.par               # step 2
```

The two `.usr` styles need genuinely different things patched, not just
the same numbers in different files — `NUMBER_ELEMENTS_Y` means "how many
near-wall rings to actually spline-correct" in `smooth_2D`'s (an
algorithmic choice, typically `1`, hence `--smoothing-levels` with no safe
default) versus "the mesh's total radial ring count" in `example_nekrs`'s
(mesh-derived, computed automatically). See the script's own docstring
(`python3 compute_case_params.py --help`) for the full picture.
`NUMBER_ELEMENTS_Z` (spanwise resolution) isn't derivable from the 2D mesh
at all — it's a free physics choice, passed to `arglist.sh` below.

## Running it end to end

```bash
# 1. Generate the 2D grid (external, not part of this repo)
#    -> {airfoil}.p3d, {airfoil}.nmf   (e.g. Construct2D/naca0012_sharp.p3d)

# 2. Convert to a linear 2D Gmsh mesh, order-elevate, and hand off to
#    gmsh2nek. write_geo_file_2d() isn't wired into the CLI (it needs the
#    physical-group ids p3d2gmsh_nek() only knows after conversion), so
#    this is a short inline script rather than one command:
python3 - <<'EOF'
from p3d_to_gmsh_nek import p3d2gmsh_nek, write_geo_file_2d
out, groups = p3d2gmsh_nek(
    p3d_file="../../Construct2D/naca0012_sharp.p3d",
    angle_of_attack=5.0,          # must match naca.par's p_aoa later
    mesh_type="CGRD",             # or "OGRD"
    output_file="airfoil.msh",
)
write_geo_file_2d("airfoil.msh", "airfoil.geo", groups, "airfoil_o2.msh", order=2)
EOF
gmsh airfoil.geo -0

gmsh2nek   # interactive: dimension=2, airfoil_o2, solid=0, periodic pairs=0, airfoil_o2

# 3. Step 1: fit + apply the wall spline correction on the 2D mesh
python3 compute_case_params.py ../../Construct2D/naca0012_sharp.nmf \
    --smoothing-levels 1 \
    --patch-step1-usr smooth_2D/naca_gen_spline_info.usr
cp airfoil_o2.re2 smooth_2D/naca_gen_spline_info.re2
cd smooth_2D
genmap   # interactively: naca_gen_spline_info, tolerance (e.g. 0.01)
./makenek naca_gen_spline_info
mpirun -np 1 ./nek5000   # writes geom_3x3_naca_gen_spline_info.dat;
                          # "EXIT: smooth_geom SUCCESS!" + nonzero exit is
                          # expected, that's how this case signals success
cd ..

# 4. Extrude the (unsmoothed) 2D mesh to 3D, and patch step 2's case files
#    in the same call -- NUMBER_ELEMENTS_Z (spanwise resolution) is a free
#    choice, set once here rather than separately in two places
./arglist.sh airfoil_o2.re2 ../../Construct2D/naca0012_sharp.nmf \
    ../high_order_mesh/example_nekrs/naca.re2 3 \
    ../high_order_mesh/example_nekrs/naca.usr \
    ../high_order_mesh/example_nekrs/naca.par
```

Continue with `../high_order_mesh/README.md` for step 2 (applying the
correction to the 3D mesh and running NekRS).

## Physical-id convention

Same fixed ids as the SOD2D pipeline's `p3d_to_gmsh.py`, defined once in
`p3d_to_gmsh_nek.py`:

| id | name |
|---|---|
| 1 | WALL |
| 2 | INLET |
| 3 | OUTLET |
| 109 | domain (surface) |

`gmsh2nek` doesn't use physical *names* for anything, only the tag number
— it becomes `boundaryID` in the `.re2` file, and real Nek5000 BC types
(`W  `, `v  `, `O  `) get assigned from that in the case's `.usr` file
(`example_nekrs/naca.usr`'s `usrdat2`). Wall/inlet/outlet are classified
the same way as the SOD2D pipeline: at each farfield-type point, the local
outward direction is compared against the free-stream direction implied
by `angle_of_attack`. For a C-grid, the two wake-end faces are always
OUTLET regardless of AoA.

## Mesh topology notes

- **O-grid**: the periodic i-closure (i=0 / i=idim-1 are the same physical
  point) is realized purely by node-id aliasing in `p3d_to_gmsh_nek.py` —
  no boundary face or Gmsh-side periodic surface is needed for it.
  `smooth_2D`'s wall-spline fit treats the wall as a closed loop, spanning
  the whole ring.
- **C-grid**: the wake-cut fold (the two i-sub-ranges the `.nmf` marks
  `ONE_TO_ONE` on the wall face) is merged the same way, by node-id
  aliasing. The real wall is only a *sub-range* of the ring — the rest is
  the wake-cut seam. `smooth_2D`'s wall-spline fit handles this as a
  generally open curve (natural/free spline end conditions); for a sharp
  trailing edge specifically, the node aliasing above makes the wall's
  start and end points coincide anyway, so it closes up without any
  special-casing. See `smooth_2D/naca_gen_spline_info.usr`'s
  `WALL_START_ELEMENT`/`WALL_ELEMENT_COUNT` and the comments directly
  above `smooth_geom0`.

## Directory layout

```
p3d_to_gmsh_nek.py     # Plot3D+NMF -> linear Gmsh .msh, wall/inlet/outlet classification
compute_case_params.py # derives ring/wall-range parameters from a .nmf, patches case files
arglist.sh              # 2D .re2 -> 3D .re2 (re2torea/n2to3/reatore2), reports/patches case params
base.rea, tail.rea      # Nek5000 .rea header/footer templates arglist.sh sandwiches around the mesh
smooth_2D/              # Fischer's smooth_geom0 wall-spline case, adapted for this repo's meshes
  naca_gen_spline_info.usr  # the algorithm itself -- see Mesh topology notes above
  naca_gen_spline_info.par  # small Nek5000 run-parameter file
  SIZE, makenek             # Nek5000 case build config
  README                    # this case's own step-by-step notes (incl. Fischer's original CASE I/II history)
```
