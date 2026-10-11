# Construct2D → STREAmS-2

This pipeline turns a Construct2D airfoil C-grid into the grid file that
[STREAmS-2](https://github.com/STREAmS-CFD/STREAmS-2) reads for its
curvilinear airfoil setup.

It is much shorter than the SOD2D and NekRS pipelines, because STREAmS-2 is
a structured finite-difference code:

- **No 3D mesh.** STREAmS-2 reads a single 2D block and extrudes it in z
  itself (`nzmax` in `singleideal.ini`, or `nzmax = 1` for a 2D run).
- **No partitioning.** STREAmS-2 splits the block across MPI ranks at run
  time (`x_split`, `y_split` in `[mpi]`).
- **No high-order wall smoothing.** The finite-difference grid points
  *are* the solution points, so there are no extra nodes to put back onto
  the curved wall.

What's left is to write the Construct2D grid in STREAmS-2's binary format,
oriented the way its airfoil boundary conditions expect.

```
airfoil .dat  ──►  Construct2D (C-grid)  ──►  p3d_to_streams.py  ──►  grid2d.xyz
                   {name}.p3d / .nmf           orient, close the       STREAmS-2
                                               wake cut exactly        grid_type = 6
```

## Running it

```bash
./generate_grid.sh <case_dir> <airfoil.dat>
```

`<case_dir>` must contain a Construct2D `grid_options.in`; the airfoil
path is relative to `<case_dir>` (or absolute). The script:

1. (Re)builds `../Construct2D/construct2d` from scratch if the binary is
   missing or older than any file in `Construct2D/src/`. Construct2D's
   Makefile doesn't track module dependencies, so incremental builds after
   a source change can fail.
2. Runs Construct2D in `<case_dir>` with the menu sequence `GRID` →
   `SMTH` → `QUIT`. This writes `{name}.p3d`, `{name}.nmf` and
   `{name}_stats.p3d`, with `{name}` taken from `grid_options.in`. The
   console output goes to `construct2d.log`.
3. Converts `{name}.p3d` to `<case_dir>/grid2d.xyz` with
   `p3d_to_streams.py` and prints the checks listed below.

The conversion can also be run on its own:

```bash
python3 p3d_to_streams.py naca0012_re23k.p3d grid2d.xyz [--ref old_grid2d.xyz]
```

### What `p3d_to_streams.py` does

- **Reads** the 2D ASCII Plot3D file (`imax jmax` header, then x and y).
- **Orients** the grid as in STREAmS-2's `examples/airfoil_*`:
  - i runs lower wake → pressure side → leading edge → suction side →
    upper wake. The grid is flipped if Construct2D wrote it the other way.
  - j runs from the wall (j = 1) to the far field.
- **Closes the wake cut exactly.** The two copies of the wake cut on
  j = 1 (points i and nx+1−i) are set to their common midpoint, so they
  are bit-identical.
- **Writes** `grid2d.xyz` as a stream binary: `int32 nx, ny, 1`, then
  `float64 x(i,j)`, then `float64 y(i,j)`, with i fastest.
- **Prints** the trailing-edge index, the number of wall points, chord and
  thickness, the first wall-normal spacing, the surface spacing range, the
  far-field radius, the maximum stretching ratios, and the sign of the
  Jacobian, which must be the same everywhere.

## Example: `examples/naca0012_re23k/`

This is the grid of the NACA 0012 transonic-buffet case (M = 0.8,
Re_c = 23 000, α = 3°, laminar 2D DNS). The case reproduces Jia, Li & Bae
(2024), *J. Aerosp. Eng.* 37(6), 04024082, doi:10.1061/JAEEEZ.ASENG-4978.

```bash
./generate_grid.sh examples/naca0012_re23k naca0012.dat     # ~1 min
```

| Setting | Value | Result |
|---|---|---|
| `nsrf`, `nwke` | 1200, 300 | nx = nsrf + 2·nwke = 1800, trailing edge at i = 301 |
| `jmax` | 512 | ny = 512 |
| `lesp`, `tesp` | 5e-4, 6e-4 | surface spacing 5.0e-4 to 2.6e-3 c |
| `radi` | 100 | far field at 113–154 c from mid-chord |
| `ypls`, `recd`, `cfrc` | 0.220895, 23000, 1.0 | first wall cell 1.5e-4 c |
| `slvr`, `topo` | `HYPR`, `CGRD` | hyperbolic C-grid |

`naca0012.dat` is the NACA 0012 with its trailing edge closed by
subtracting x·t(1) from the half-thickness t(x). It has 1501 cosine-spaced
points per side, from this generator:

```python
import numpy as np
t = lambda x: 0.6*(0.2969*np.sqrt(x) - 0.1260*x - 0.3516*x**2 + 0.2843*x**3 - 0.1015*x**4)
yt = lambda x: t(x) - x*t(1.0)
x = 0.5*(1 - np.cos(np.linspace(0, np.pi, 1501)))
xs = np.r_[x[::-1], x[1:]]; ys = np.r_[yt(x[::-1]), -yt(x[1:])]
with open('naca0012.dat', 'w') as f:
    f.write('NACA0012 closed TE (sod2d)\n')
    for a, b in zip(xs, ys): f.write(f'{a:.12f} {b:.12f}\n')
```

### Setting the wall spacing

Construct2D sets the first wall-normal cell from `ypls` with a turbulent
flat-plate skin-friction estimate:

```
Cf = (2 log10(Re) - 0.65)^(-2.3),    y0 = ypls · cfrc / (Re · sqrt(Cf/2))
```

`cfrc` defaults to 0.5. To get a target `y0` (in chords):

```bash
python3 -c "import math; Re=23000; y0=1.5e-4; cfrc=1.0; Cf=(2*math.log10(Re)-0.65)**-2.3; print(y0*Re*math.sqrt(Cf/2)/cfrc)"
```

For laminar flows this "y+" is only a convenient knob, not the actual y+.
Check the spacing `p3d_to_streams.py` reports.

In STREAmS-2 the explicit time step is usually limited by wall-normal
diffusion in the first cell, so dt ∝ y0². On this grid 1.5e-4 c keeps the
wall cell below one viscous unit and allows about a 10× larger time step
than 5e-5 c.

## Using the grid in STREAmS-2

Copy `grid2d.xyz` into the run directory. The settings in
`singleideal.ini` that have to match the grid:

```ini
[grid]
grid_dim  = 2
grid_type = 6          ; read grid2d.xyz
nxmax     = 1800       ; nx printed by p3d_to_streams.py
nymax     = 512        ; ny
nzmax     = 1          ; 2D run (or the number of spanwise points)

[bc]
xmin = 4               ; wake outflow (lower)
xmax = 4               ; wake outflow (upper)
ymin = 11              ; airfoil wall + wake cut
ymax = 3               ; far field

[mpi]
x_split = 120          ; must divide nxmax
y_split = 1            ; the airfoil sponge needs a single block in j
```

The airfoil wall (`ymin = 11`) is isothermal in upstream STREAmS-2. The
fork used here adds `adiabatic_wall = 1` under `[flow]` (see below).

## `CFD_solvers/STREAmS-2/`

STREAmS-2 is included as a git submodule at the repo root,
[`CFD_solvers/STREAmS-2`](../CFD_solvers). It points to the fork
[kvishalrai/STREAmS-2](https://github.com/kvishalrai/STREAmS-2), branch
`aerohexmesh`, which is upstream STREAmS-2 plus three commits:

1. **Adiabatic airfoil wall.** `[flow] adiabatic_wall = 1` extrapolates
   the wall temperature with zero normal gradient to third order,
   T_w = (18 T₂ − 9 T₃ + 2 T₄)/11, and mirrors the ghost temperatures.
   The default 0 keeps the isothermal wall. CPU and GPU backends.
2. **Runtime airfoil forces every N steps.**
   `[curvi] enable_forces_runtime = N` writes `airfoil_forces_runtime.dat`
   every N steps instead of every step.
3. **Build configuration.**
   - A `GPU_CC` make variable sets the compute capability for `nv-cuda`
     builds without `CUDA_VER`, e.g. `GPU_CC = cc90` for H100.
   - New MareNostrum 5 templates in `code/make-templates/`:
     `marenostrum5_cuda_nvhpc243_cc90` (ACC),
     `marenostrum5_gpp_cpu_intel_impi` (GPP) and
     `marenostrum5_grace_cpu_gnu` (Grace).

To build it:

```bash
cd CFD_solvers/STREAmS-2/code
cp make-templates/makefile.inc.marenostrum5_gpp_cpu_intel_impi makefile.inc
make
```

See STREAmS-2's own [README](../CFD_solvers/STREAmS-2/README.md) and
documentation for everything else.

## Prerequisites

- gfortran and make, to build Construct2D.
- Python 3 with `numpy`.
- STREAmS-2 itself, built from the submodule as above, to run the case.

Generated files are not tracked: `*.p3d`, `*.nmf`, `grid2d.xyz` and
`construct2d.log` in `examples/*/`. Regenerate them with
`generate_grid.sh`.
