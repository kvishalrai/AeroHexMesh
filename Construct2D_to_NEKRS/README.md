# Construct2D → NekRS

This is the **middle and final stage** of the pipeline described in the
[top-level README](../README.md): it takes the 2D `.p3d`/`.nmf` grid that
Construct2D produces and turns it into a spanwise-periodic, wall-spline-
smoothed 3D mesh that NekRS can run on, via Nek5000's own mesh tools.

```
2D .p3d/.nmf   ──►  linear_mesh/         ──►  high_order_mesh/example_nekrs/
(Construct2D)       2D convert, smooth        apply the 2D correction to the
                     the wall, extrude 3D      real 3D mesh, run NekRS
```

Currently set up and verified for a **C-grid** mesh; the tooling also
supports O-grid (see [`compute_case_params.py`](#computing-case-parameters)
below), it just isn't the case checked in here.

No large generated mesh/geometry/binary files are committed anywhere in
this pipeline (mirrors the top-level README's policy) — everything below
regenerates them from the small `.usr`/`.par`/script files that are.

## Stages

1. **[`linear_mesh/`](linear_mesh)** — converts the native 2D Construct2D
   grid into a linear Gmsh mesh with wall/inlet/outlet boundary groups
   (`p3d_to_gmsh_nek.py`), which Gmsh order-elevates and Nek5000's
   `gmsh2nek` turns into a `.re2`. `smooth_2D/` (Fischer's `smooth_geom0`
   algorithm, adapted for this repo's own meshes — see
   [Credits](../README.md#credits)) fits a cubic spline through the wall's
   coarse control points and resamples it at the mesh's actual high-order
   nodes, replacing the faceted wall approximation with the true smooth
   curve; that correction is stored compactly (`geom_3x3_*.dat`), one
   entry per element, for every element in the 2D plane. `arglist.sh` then
   extrudes the *original* (unsmoothed) 2D mesh into a periodic-in-Z 3D
   mesh via Nek5000's `re2torea`/`n2to3` tools.
2. **[`high_order_mesh/example_nekrs/`](high_order_mesh/example_nekrs)** —
   a NekRS case. Its `usrdat2`/`get_smooth_data` reads the compact
   correction from step 1 and stamps it onto every element of the real 3D
   mesh (repeating across all spanwise levels, since the mesh is just a
   periodic Z-extrusion of one 2D pattern), then runs the actual flow
   solve.

## Computing case parameters

Both `smooth_2D/naca_gen_spline_info.usr` and
`high_order_mesh/example_nekrs/naca.usr` (+ `naca.par`) need to agree on
several mesh-derived numbers — how many elements are in one ring, and for
a C-grid, exactly which sub-range of that ring is the real wall (the rest
being the wake-cut seam). Don't work these out by hand or by inspecting
Gmsh output: **`linear_mesh/compute_case_params.py`** reads them directly
from the airfoil's `.nmf` file (its `VISCOUS` boundary entry already
encodes the wall's exact element range, for O-grid or C-grid alike) and
can patch them straight into the case files. See the script's own
docstring (`python3 compute_case_params.py --help`) for the full
step1/step2 distinction — they need different things patched, not just the
same numbers in different files.

## Running it end to end

```bash
cd linear_mesh

# 1. Generate the 2D grid (external, not part of this repo)
#    -> {airfoil}.p3d, {airfoil}.nmf   (e.g. Construct2D/naca0012_sharp.p3d)

# 2. Convert to a linear 2D Gmsh mesh, order-elevate, and hand off to
#    gmsh2nek. write_geo_file_2d() isn't wired into the CLI (it needs the
#    physical-group ids p3d2gmsh_nek() only knows after conversion), so
#    this is a short inline script rather than one command:
python3 - <<'EOF'
from p3d_to_gmsh_nek import p3d2gmsh_nek, write_geo_file_2d
out, groups = p3d2gmsh_nek(
    p3d_file="../Construct2D/naca0012_sharp.p3d",
    angle_of_attack=5.0,          # must match naca.par's p_aoa later
    mesh_type="CGRD",             # or "OGRD"
    output_file="airfoil.msh",
)
write_geo_file_2d("airfoil.msh", "airfoil.geo", groups, "airfoil_o2.msh", order=2)
EOF
gmsh airfoil.geo -0

gmsh2nek   # interactive: dimension=2, airfoil_o2, solid=0, periodic pairs=0, airfoil_o2

# 3. Step 1: fit + apply the wall spline correction on the 2D mesh
python3 compute_case_params.py ../Construct2D/naca0012_sharp.nmf \
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

# 4. Extrude the (unsmoothed) 2D mesh to 3D, and patch step-2's case files
#    in the same call -- NUMBER_ELEMENTS_Z (spanwise resolution) is a free
#    choice, set once here rather than separately in two places
./arglist.sh airfoil_o2.re2 ../Construct2D/naca0012_sharp.nmf \
    ../high_order_mesh/example_nekrs/naca.re2 3 \
    ../high_order_mesh/example_nekrs/naca.usr \
    ../high_order_mesh/example_nekrs/naca.par

# 5. Step 2: copy the spline correction in, then run NekRS
cp smooth_2D/geom_3x3_naca_gen_spline_info.dat \
    ../high_order_mesh/example_nekrs/geom_3x3_naca_gen_spline_info.dat
cd ../high_order_mesh/example_nekrs
sbatch airfoil0.sh
```

Double-check before running step 5: `naca.par`'s `[VELOCITY] viscosity`
(set to `-Re`) and `[CASEDATA] p_aoa`/`p_u0` match your intended physics —
these are free choices `compute_case_params.py` never touches, only the
mesh-derived `NUMBER_ELEMENTS_*` are automated. Also check `naca.oudf`'s
`codedFixedValueVelocity`: it's a hardcoded `aoa`/`cos`/`sin`, not wired to
`naca.par`'s `p_aoa` (`kernelInfo` defines set in `UDF_LoadKernels` aren't
visible in a coded-BC kernel's own compile context) — keep both in sync by
hand.

## `CFD_code/`

`CFD_code/Nek5000` and `CFD_code/nekRS` are included as git submodules
(see the top-level README's [Credits](../README.md#credits)). Both stages
above depend on tools built from Nek5000:

| Tool | Used by | Build notes |
|---|---|---|
| `gmsh2nek` | step 2 above (2D Gmsh → `.re2`) | `tools/maketools gmsh2nek` |
| `genmap` | `smooth_2D/` (step 3) | `tools/maketools genmap` |
| `re2torea`, `reatore2`, `n2to3` | `arglist.sh` (step 4) | `tools/maketools re2torea reatore2 n2to3` (`re2torea`/`reatore2` share one makefile — one `maketools` call builds both) |

NekRS itself (`CFD_code/nekRS`) is a much larger build (OCCA JIT, GPU
backend) — see its own [build docs](CFD_code/nekRS/README.md) or use a
site-provided install (e.g. an existing `NEKRS_HOME`) as `high_order_mesh/example_nekrs/airfoil0.sh`
does.

## Prerequisites

- Python 3 with `numpy`.
- [Gmsh](https://gmsh.info/), invoked as a CLI.
- MPI, for Nek5000/NekRS.
- Nek5000 tools built as above, and a NekRS install — see
  [`CFD_code/`](#cfd_code).

Generated mesh, geometry-correction, and run output files (`*.re2`,
`geom_3x3*.dat`, `*.ma2`, compiled binaries, `.cache/`, logs, restart
fields) are gitignored — regenerate them by running the pipeline above
rather than expecting them to be present after a clone.
