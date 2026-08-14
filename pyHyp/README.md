# pyHyp: hyperbolic surface extrusion

This is a second, independent front end for building 3D meshes, alongside
the `Construct2D_to_<SOLVER>/` pipelines described in the
[top-level README](../README.md). It solves a different problem:

- `Construct2D_to_<SOLVER>/` starts from a single 2D airfoil cross-section
  and **spanwise-extrudes** it — the same shape repeated along a straight
  line. Good for studying a 2D airfoil section in a 3D solver, but it can't
  represent a real 3D shape like a tapered or swept wing.
- `pyHyp/` starts from an actual 3D surface mesh (e.g. a wing, wall-to-wall,
  with taper/sweep/twist) and grows a volume mesh outward from that surface,
  layer by layer, along the local surface normal — a technique called
  **hyperbolic marching**. Each new layer's spacing and shape are computed
  from the previous layer, not fixed in advance, so the mesh naturally
  follows the surface's real curvature out to a far-field boundary.

This directory only covers **surface → volume mesh generation**, using the
[pyHyp](https://github.com/mdolab/pyhyp) library (MDO Lab). Importing the
resulting mesh into SOD2D or NekRS is not implemented yet.

## Prerequisites

pyHyp is a compiled Fortran/PETSc tool (via f2py), not a plain pip package —
building it means building against a working CGNS + PETSc + MPI stack
first, similar in spirit to how this repo already builds Nek5000/SOD2D from
source. The recipe below uses a dedicated conda environment, since CGNS
isn't available as a system module on MN5:

```bash
module load anaconda/2023.07
conda create -n pyhyp-env -c conda-forge --solver=libmamba -y \
    python=3.11 cgns petsc=3.20 petsc4py mpi4py openmpi compilers cmake numpy pip
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate pyhyp-env
pip install mdolab-baseclasses scipy tabulate
```

(`--solver=libmamba` matters here — conda's classic solver can hang or get
killed on a dependency tree this size; libmamba solves the same environment
in well under a minute.)

`cgnsutilities` and `pyhyp` are **not** on PyPI (`pip install cgnsutilities`
fails outright) — both ship a Fortran/f2py extension that has to be built
against this environment's CGNS/PETSc first:

```bash
git clone --branch v2.9.0 --depth 1 https://github.com/mdolab/cgnsutilities.git
cd cgnsutilities
cp config/defaults/config.LINUX_GFORTRAN.mk config/config.mk
CGNS_HOME=$CONDA_PREFIX make
pip install --no-deps .
cd ..

git clone --branch v2.6.3 --depth 1 https://github.com/mdolab/pyhyp.git
cd pyhyp
cp config/defaults/config.LINUX_GFORTRAN.mk config/config.mk
CGNS_HOME=$CONDA_PREFIX PETSC_DIR=$CONDA_PREFIX PETSC_ARCH= make
pip install --no-deps .
```

Verify with `python3 -c "from pyhyp import pyHyp"` — a successful build
prints "Module hyp was successfully imported" partway through `make`.

Neither `pyhyp` nor `cgnsutilities` is vendored into this repo (both are
built from source, like the rest of this section) — clone them wherever you
build your environment, not into this directory.

### Cluster-specific module/environment setup

`example_m6_wing/generate_mesh.sh` runs whatever shell command is in the
`AEROHEXMESH_PYHYP_MODULE_SETUP` environment variable before calling pyHyp,
the same externalized-setup pattern used by the `Construct2D_to_<SOLVER>/`
pipelines. **Set it yourself** if you're not on MN5:

```bash
export AEROHEXMESH_PYHYP_MODULE_SETUP="module load anaconda/2023.07 && source \"\$(conda info --base)/etc/profile.d/conda.sh\" && conda activate pyhyp-env"
```

**On BSC MareNostrum 5**, this is already the default — no need to set it,
as long as your environment is named `pyhyp-env` as above.

### Run it on a compute node, not the login node

The hyperbolic marching itself is genuinely CPU-heavy — the login node
kills any single process past 300s of CPU time, which silently cuts off a
run partway through (no error message, no traceback, just the process
disappearing). Submit it as a SLURM job instead:
`example_m6_wing/generate_mesh.job` does this (`sbatch generate_mesh.job`)
using the CPU debug queue; adapt the `#SBATCH` lines if you're not on MN5.

## Usage

```bash
python3 generate_volume_mesh.py \
    --config config.json \
    --input-file surface.cgns \
    --file-type CGNS \
    --output-file volumeMesh.xyz
```

`--config` is a JSON file of pyHyp's own options (grid, pseudo-grid, and
smoothing parameters — `N`, `s0`, `marchDist`, `epsE`, etc.; see
[pyHyp's options docs](https://mdolab-pyhyp.readthedocs-hosted.com/en/latest/options.html)
for the full list) — everything except the input/output file, which are
CLI flags instead so the same config can be reused across different surface
meshes. `--file-type` is `CGNS` (default) or `PLOT3D`, matching the input
surface mesh's own format. Output is always written as Plot3D.

## Example

`example_m6_wing/` extrudes a volume mesh around the ONERA M6 wing (a
standard CFD validation geometry), from a bundled Plot3D surface mesh
(`m6_small.fmt`, taken from pyHyp's own test suite) — no download needed.
Good for checking your build works before pointing it at your own surface:

```bash
cd example_m6_wing
sbatch generate_mesh.job
```

Confirms it worked the same way `Construct2D_to_SOD2D/high_order_mesh/`
does: the run log shows all `N` grid levels completing (81 here) without
error, `Min Quality` staying well above 0, and the output mesh
(`volumeMesh.xyz`) opens in ParaView/Tecplot as a smooth hex volume
wrapping the wing surface out to a farfield boundary.

## Directory layout

```
generate_volume_mesh.py   # CLI driver: config JSON + surface mesh -> Plot3D volume mesh
example_m6_wing/
  m6_small.fmt             # bundled Plot3D surface mesh (ONERA M6 wing)
  config.json               # pyHyp grid/pseudo-grid/smoothing options for this example
  generate_mesh.sh           # runs generate_volume_mesh.py on the example
  generate_mesh.job          # SLURM wrapper around generate_mesh.sh
```
