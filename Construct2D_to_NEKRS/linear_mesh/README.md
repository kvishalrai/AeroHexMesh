# p3d-to-nekrs

This stage takes the flat 2D grid Construct2D produced and turns it into a
3D mesh with a smooth wall surface, ready for NekRS's own case in
`../high_order_mesh/example_nekrs/`. If you haven't read it yet, the
[top-level README's "big picture" section](../README.md#the-big-picture)
explains *why* each of these steps exists before you dive into the *how*
below.

## What's in a `.p3d`/`.nmf` pair?

Construct2D writes two files per mesh, and this whole stage starts from
them:

- **`.p3d`** (Plot3D format) — just a big grid of `(x, y)` coordinates,
  arranged in rows and columns like a spreadsheet. One "direction" of the
  grid runs around the airfoil (circumferentially); the other runs away
  from it, from the wall out to the far boundary.
- **`.nmf`** (Neutral Map File) — a small text file that says what each
  edge of that grid of points *means*: which edge is the solid wall,
  which is the inflow, which is the outflow, and (for a C-grid) which two
  strips of points are actually the same physical wake line, folded onto
  itself.

Nothing in this pipeline modifies these two files — they're read once, at
the very start, and everything downstream is derived from them.

## The algorithm

This whole pipeline (this directory *and* `../high_order_mesh/`) is four
steps, each one taking the previous step's output as input. This
directory does the first three; `../high_order_mesh/example_nekrs/` does
the fourth and last:

```
INPUT: airfoil.p3d, airfoil.nmf           (from Construct2D)

STEP 1  Convert to a mesh format Nek5000's tools understand,
        tagging which boundary is the wall/inlet/outlet.
        →  a 2D mesh file (.re2), still with a jagged (faceted) wall

STEP 2  Fit a smooth curve through the wall's corner points, and
        record where that curve actually is at every mesh node.
        →  a small data file with the correction (not applied to
           the mesh yet -- that's step 4)

STEP 3  Stretch the 2D mesh into a 3D mesh by repeating it along a
        straight line (the wing's span), telling the solver the flow
        repeats at every cross-section, and fill in the matching
        settings (how many mesh elements, etc.) in the NekRS case
        that will use this mesh.
        →  a 3D mesh file (.re2), wall still jagged, and a
           correctly-configured NekRS case ready for it

STEP 4  (in ../high_order_mesh/ -- see its own README) Stamp the
        smooth-wall correction from step 2 onto the real 3D mesh
        from step 3, then run the actual airflow simulation.
```

Unlike the SOD2D pipeline (which has one script, `run_pipeline.py`, that
does everything), these steps are run one at a time — see
[Running it end to end](#running-it-end-to-end) below for the exact
commands. Why: Nek5000's mesh tools are old, single-purpose command-line
programs (each one only reads/writes one specific file format), not a
single library you can call from Python — so building the mesh really is
a short pipeline of separate programs, not a design choice made here.

## Prerequisites

- Python 3 with `numpy`.
- [Gmsh](https://gmsh.info/), invoked as a CLI.
- MPI, for Nek5000.
- `gmsh2nek`, `genmap`, `re2torea`, `reatore2`, `n2to3` built from the
  `../../CFD_solvers/Nek5000` submodule — see the top-level
  [`Construct2D_to_NEKRS/README.md`](../README.md#cfd_solversnek5000-and-cfd_solversnekrs).

## Computing case parameters

Two files need to agree with each other and with the actual mesh, or the
correction from step 2 will get stamped onto the wrong elements later:
`smooth_2D/naca_gen_spline_info.usr` (step 2) and
`../high_order_mesh/example_nekrs/naca.usr` + `naca.par` (step 4). Both
need to know, in the language of "mesh elements" rather than physical
coordinates: how many elements make up one full ring of the mesh, and —
for a C-grid, where the wall is only *part* of that ring — exactly which
elements are the real wall versus the wake-cut region on either side of
it.

Don't work these numbers out by hand (e.g. by opening the mesh in a
viewer and counting). **`compute_case_params.py`** reads them straight
out of the `.nmf` file — its `VISCOUS` entry already records the wall's
exact element range, for O-grid or C-grid alike — and can write them
directly into the case files for you:

```bash
python3 compute_case_params.py <airfoil>.nmf                                  # just show the numbers
python3 compute_case_params.py <airfoil>.nmf --smoothing-levels 1 \
    --patch-step1-usr smooth_2D/naca_gen_spline_info.usr                       # fill in step 2's case
python3 compute_case_params.py <airfoil>.nmf --nz 3 \
    --patch-step2-usr ../high_order_mesh/example_nekrs/naca.usr \
    --patch-step2-par ../high_order_mesh/example_nekrs/naca.par               # fill in step 4's case
```

One subtlety: a setting called `NUMBER_ELEMENTS_Y` appears in *both*
`.usr` files, but it means two different things, so the script needs to
be told explicitly rather than guessing:

- In step 2's case, it means "how many rings next to the wall should
  actually get the smoothing correction" — a modeling choice, not
  something the mesh tells you, typically `1` (just the ring touching the
  wall). That's why you must pass `--smoothing-levels` yourself.
- In step 4's case, it means "how many rings does the mesh have in
  total" — that *is* derivable from the mesh, so the script fills it in
  automatically.

`NUMBER_ELEMENTS_Z` (how many mesh elements along the wing's span) isn't
derivable from the 2D mesh at all — it's your own choice of spanwise
resolution, given directly to `arglist.sh` in step 3 below. Run
`python3 compute_case_params.py --help` for the complete picture.

## Running it end to end

```bash
# STEP 0 -- Generate the 2D grid with Construct2D itself (a separate tool,
# not part of this repo). You should end up with a pair of files like
# Construct2D/naca0012_sharp.p3d and Construct2D/naca0012_sharp.nmf.

# STEP 1 -- Convert to a linear (straight-sided) 2D mesh, add the extra
# curved-edge nodes Nek5000's tools need (Gmsh's "order elevation"), and
# hand off to gmsh2nek. This is a short inline Python script rather than
# one command because write_geo_file_2d() needs the boundary-tag ids that
# p3d2gmsh_nek() only produces after conversion -- there's no way to know
# them beforehand:
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

gmsh2nek   # interactive prompts: dimension=2, airfoil_o2, solid=0, periodic pairs=0, airfoil_o2
           # -> writes airfoil_o2.re2, the 2D mesh in Nek5000's own format

# STEP 2 -- Fit the smooth wall curve on the 2D mesh. First fill in this
# case's settings, then feed it the mesh from step 1, then build and run it:
python3 compute_case_params.py ../../Construct2D/naca0012_sharp.nmf \
    --smoothing-levels 1 \
    --patch-step1-usr smooth_2D/naca_gen_spline_info.usr
cp airfoil_o2.re2 smooth_2D/naca_gen_spline_info.re2
cd smooth_2D
genmap   # interactive prompts: naca_gen_spline_info, tolerance (e.g. 0.01)
./makenek naca_gen_spline_info
mpirun -np 1 ./nek5000
    # -> writes geom_3x3_naca_gen_spline_info.dat, the wall correction data.
    #    "EXIT: smooth_geom SUCCESS!" followed by a nonzero shell exit is
    #    expected -- that's just how this program signals success, not a
    #    crash. Look for the success message and the output file, not the
    #    exit code.
cd ..

# STEP 3 -- Stretch the (still jagged-walled) 2D mesh into a 3D mesh, and
# fill in step 4's case settings (in ../high_order_mesh/example_nekrs/),
# in one call. NUMBER_ELEMENTS_Z (how finely to resolve the wing's span)
# is a free choice you make here -- 3 in this example -- rather than
# something derived from the 2D mesh:
./arglist.sh airfoil_o2.re2 ../../Construct2D/naca0012_sharp.nmf \
    ../high_order_mesh/example_nekrs/naca.re2 3 \
    ../high_order_mesh/example_nekrs/naca.usr \
    ../high_order_mesh/example_nekrs/naca.par
```

Continue with `../high_order_mesh/README.md` for step 4, the final one:
applying the wall correction to this 3D mesh and actually running NekRS.

## How the wall gets tagged (physical-id convention)

Every boundary edge/element in the mesh gets labeled with one of these
fixed numeric ids, defined once in `p3d_to_gmsh_nek.py` (same convention
the SOD2D pipeline's `p3d_to_gmsh.py` uses):

| id | name | meaning |
|---|---|---|
| 1 | WALL | the solid airfoil surface |
| 2 | INLET | boundary where flow enters the domain |
| 3 | OUTLET | boundary where flow leaves the domain |
| 109 | domain | not a boundary -- the 2D surface (interior) itself |

Gmsh and `gmsh2nek` don't care about human-readable names, only these
numbers — the id becomes `boundaryID` in the `.re2` file, and it's what
the NekRS case later uses (in `naca.usr`'s `usrdat2`) to decide the real
physics boundary condition for each face (no-slip wall, velocity inlet,
pressure outlet).

Inlet vs. outlet isn't hardcoded — it's worked out geometrically: at each
point on the outer/far boundary, the code compares the local outward
direction against the free-stream direction implied by
`angle_of_attack`, and calls it an inlet if flow would be entering there,
outlet if leaving. For a C-grid, the two short edges at the open ends of
the "C" (where the wake gets cut) are always OUTLET, regardless of angle
of attack.

## O-grid vs. C-grid, in this code

Both grid shapes are handled by the *same* conversion script; the
difference is entirely in the shape of the input `.nmf` file. See the
[top-level README](../README.md#the-big-picture) for the plain-language
picture; here's the implementation detail for each:

- **O-grid**: the mesh closes on itself — the first and last "column" of
  points around the airfoil are actually the same physical points. Rather
  than storing them twice, the code just reuses the same point ids for
  both ("node-id aliasing"), so the mesh comes out seamless with no extra
  bookkeeping. The wall in this case is the *entire* ring going around
  the airfoil, and `smooth_2D`'s wall-spline fit treats it as one closed
  loop.
- **C-grid**: the two strips of points where the wake gets "cut" (marked
  `ONE_TO_ONE` in the `.nmf`) are merged the same way, by node-id
  aliasing — not by any special Gmsh feature. The real wall here is only
  a *sub-range* of the ring (the rest of that ring is the wake-cut seam,
  not solid surface), so `smooth_2D`'s wall-spline fit treats it as an
  open curve in general. For an airfoil with a sharp trailing edge
  specifically, the node aliasing above happens to make the wall's start
  and end points the same physical point anyway, so the curve closes up
  on its own with no special-casing needed. See
  `smooth_2D/naca_gen_spline_info.usr`'s `WALL_START_ELEMENT`/
  `WALL_ELEMENT_COUNT` and the comments directly above `smooth_geom0` for
  exactly how that works in code.

## Directory layout

```
p3d_to_gmsh_nek.py     # STEP 1: Plot3D+NMF -> linear Gmsh mesh, with wall/inlet/outlet tags
compute_case_params.py # reads ring/wall-range settings from a .nmf, fills them into case files
arglist.sh              # STEP 3+4: 2D mesh -> 3D mesh, then fills in step 4's case settings
base.rea, tail.rea      # small Nek5000 template files arglist.sh needs internally (see its own comments)
smooth_2D/              # STEP 2: Paul Fischer's wall-smoothing case, adapted for this repo's meshes
  naca_gen_spline_info.usr  # the smoothing algorithm itself -- see "O-grid vs. C-grid" above
  naca_gen_spline_info.par  # small Nek5000 run-parameter file
  SIZE, makenek             # Nek5000 case build settings
  README                    # this case's own step-by-step notes
```
