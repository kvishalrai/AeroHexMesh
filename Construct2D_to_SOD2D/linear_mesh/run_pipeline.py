"""Orchestrates the full p3d -> Gmsh -> SOD2D mesh pipeline for one airfoil.

Given a user-supplied 2D Plot3D + Neutral Map File pair (produced
externally with Construct2D) and a flow_config*.json describing the mesh
(mesh_type, spanwise extrusion, ...; the 2D mesh dimensions are read
straight from the NMF, not from config), this script:

1. extrudes the 2D mesh spanwise into 3D and remaps the NMF (mesh_extrusion.py)
2. converts the extruded Plot3D mesh to Gmsh format (p3d_to_gmsh.py)
3. writes the periodic Gmsh .geo file and runs Gmsh on it
4. exports the Gmsh mesh to SOD2D's HDF5 format (sod2d_tools/gmsh2sod2d.py)
5. partitions the mesh for parallel SOD2D runs (sod2d_tools/tool_meshConversorPar)
"""

from pathlib import Path
import json
import os
import subprocess

from mesh_extrusion import (
    sod2d_mesh,
    write_ext_nmf,
    write_geo_file,
    write_partition_input_json,
    read_plot3d_2d,
    PERIODIC_ID,
)
from p3d_to_gmsh import p3d2gmsh
from wall_spline import build_wall_spline_table

# Environment setup for the external tools this pipeline shells out to
# (Gmsh -- via python3's own "gmsh" package, MPI/HDF5 for partitioning).
# Defaults to the modules that give a working Gmsh + MPI + HDF5 + python3
# environment on BSC MareNostrum 5. On any OTHER system, this default is
# meaningless -- set AEROHEXMESH_MODULE_SETUP yourself instead of editing
# this file:
#
#   export AEROHEXMESH_MODULE_SETUP="module load gmsh openmpi hdf5 python3"
#
# (whatever module names/commands actually load a working Gmsh + MPI +
# HDF5 + python3 on your system). Set it to "" for a no-op if you've
# already activated everything yourself (e.g. a virtualenv) before
# running this script.
MODULE_SETUP = os.environ.get(
    "AEROHEXMESH_MODULE_SETUP",
    "module purge && module load bsc/1.0 nvidia-hpc-sdk/24.3 "
    "hdf5/1.14.1-2-nvidia-nvhpcx mkl/2025.2 python/3.12.1-gcc cmake/3.30.5",
)
# Extra setup needed only for steps that run a python3 script requiring
# numpy/h5py (gmsh2sod2d.py, the partitioner's input.json driver): on MN5,
# the base module set's own python3 lacks those, so swap in anaconda.
PYTHON_MODULE_SETUP = os.environ.get(
    "AEROHEXMESH_PYTHON_MODULE_SETUP", "module unload python\nmodule load anaconda"
)


def valid_file(path, min_size=1024):
    """Return True when path is a regular file larger than min_size bytes."""
    path = Path(path)
    return path.is_file() and path.stat().st_size > min_size


def read_sod2d_config(config_file):
    """Load the SOD2D configuration located next to this script."""
    script_dir = Path(__file__).resolve().parent
    config_path = script_dir / config_file

    if not config_path.is_file():
        raise FileNotFoundError(
            f"SOD2D configuration file was not found: {config_path}"
        )

    with config_path.open("r", encoding="utf-8") as file:
        return json.load(file)


def run_bash(command, work_dir):
    """Run a Bash command in work_dir and raise an error if it fails."""
    subprocess.run(
        command,
        cwd=work_dir,
        shell=True,
        executable="/bin/bash",
        check=True,
    )


def sod2d_airfoil(
    env_id,
    airfoil_file,
    work_dir,
    config_file,
):
    """
    Generate and partition a SOD2D mesh for an airfoil.

    Expects a user-supplied 2D Plot3D + Neutral Map File pair
    (``{airfoil_name}.p3d`` / ``{airfoil_name}.nmf``) already generated
    with Construct2D and placed in work_dir; this function does not run
    Construct2D itself.

    Parameters
    ----------
    env_id : int
        Environment identifier used for logging.
    airfoil_file : str or Path
        Airfoil filename or basename expected by the mesh utilities.
    work_dir : str or Path
        Directory in which mesh-generation files will be created.
    config_file : str or Path
        JSON configuration file (e.g. flow_config_ogrd.json /
        flow_config_cgrd.json) located next to this script, unless an
        absolute path is supplied.
    """
    sod2d_config = read_sod2d_config(config_file)

    script_dir = Path(__file__).resolve().parent
    sod2d_tools = script_dir / "sod2d_tools"

    work_dir = Path(work_dir).expanduser().resolve()
    work_dir.mkdir(parents=True, exist_ok=True)

    gmsh2sod2d = sod2d_tools / "gmsh2sod2d.py"
    tool_mesh = sod2d_tools / "tool_meshConversorPar"

    if not gmsh2sod2d.is_file():
        raise FileNotFoundError(
            f"gmsh2sod2d.py was not found: {gmsh2sod2d}"
        )

    if not tool_mesh.is_file():
        raise FileNotFoundError(
            f"Mesh partitioning executable was not found: {tool_mesh}"
        )

    # These utilities use the airfoil name to construct output filenames.
    airfoil_name = Path(airfoil_file).stem

    original_cwd = Path.cwd()

    print(f"[environment {env_id}] Working directory: {work_dir}")
    print(f"[environment {env_id}] Generating mesh for: {airfoil_name}")

    try:
        # Some imported mesh utilities may create output relative to the
        # current directory, so temporarily enter the environment directory.
        os.chdir(work_dir)

        mesh_type = sod2d_config["mesh_type"]

        # Extrude the user-supplied 2D Plot3D mesh spanwise into 3D.
        sod2d_mesh(
            airfoil_file,
            work_dir,
            sod2d_config["z_spanwise_len"],
            sod2d_config["z_spanwise_planes"],
            mesh_type,
        )

        # Generate the connectivity map, remapped from the user-supplied
        # 2D NMF onto the extruded 3D block. idim/jmax/boundaries_2d are
        # read from the 2D NMF's own header/body here, not from config --
        # reused below rather than re-parsing the file.
        idim, jmax, boundaries_2d = write_ext_nmf(
            airfoil_file,
            work_dir,
            sod2d_config["z_spanwise_planes"],
            mesh_type,
        )

        # Fit an open cubic spline through the wall's original linear corner
        # points (pre-extrusion 2D .p3d/.nmf, untouched by sod2d_mesh()) and
        # write a compact corner+coefficient table, consumed by SOD2D's
        # MeshElasticitySolver to place each high-order boundary node
        # directly and parametrically (arc-length interpolation, no nearest-
        # point search) onto the wall's own shape without moving any corner
        # node.
        x2d, y2d = read_plot3d_2d(work_dir / f"{airfoil_name}.p3d")
        wall_spline_file = work_dir / f"{airfoil_name}_wall_spline.dat"
        build_wall_spline_table(
            x2d, y2d, boundaries_2d, mesh_type, wall_spline_file
        )
        print(f"[environment {env_id}] Wrote wall spline table: {wall_spline_file}")

        # Convert the Plot3D mesh to Gmsh format. groups is the realized
        # (dim,id,name) list -- inlet/outlet ids depend on per-element
        # geometry (see p3d_to_gmsh.GmshFile._classify_inlet_outlet) and
        # can't be known ahead of time.
        _, groups = p3d2gmsh(
            p3d_file=f"{airfoil_name}_ext.p3d",
            idim=idim,
            angle_of_attack=sod2d_config["angle_of_attack"],
            mesh_type=mesh_type,
            map_file=f"{airfoil_name}_ext.nmf",
            output_file=f"{airfoil_name}.msh",
        )

        # Generate the periodic Gmsh geometry file.
        write_geo_file(
            airfoil_file,
            work_dir,
            sod2d_config["z_spanwise_len"],
            sod2d_config["porder"],
            groups,
        )

        # Process the Gmsh geometry.
        run_bash(
            f"""
set -e

{MODULE_SETUP}

gmsh airfoil_per.geo -0
""",
            work_dir,
        )

        # Convert the Gmsh mesh to the SOD2D HDF5 format.
        run_bash(
            f"""
set -e

{MODULE_SETUP}
{PYTHON_MODULE_SETUP}

python3 "{gmsh2sod2d}" \
    "{airfoil_name}_per" \
    -p "{PERIODIC_ID}" \
    -r "{sod2d_config['porder']}"
""",
            work_dir,
        )

        # Create the partitioner's input.json file.
        write_partition_input_json(
            airfoil_file,
            work_dir,
            sod2d_config["num_partitions"],
        )

        input_json = work_dir / "input.json"
        if not input_json.is_file():
            raise FileNotFoundError(
                f"Partition input file was not generated: {input_json}"
            )

        # Partition the SOD2D mesh.
        run_bash(
            f"""
set -e

{MODULE_SETUP}
{PYTHON_MODULE_SETUP}

export OMPI_MCA_io=romio321
export OMP_NUM_THREADS=1
export SLURM_CPU_BIND=none

echo "Working directory: $(pwd)"
echo "Mesh partitioner: {tool_mesh}"
ls -lh input.json

mpirun --bind-to none \
    -np "{sod2d_config['num_partitions']}" \
    "{tool_mesh}" \
    input.json
""",
            work_dir,
        )

        print(f"[environment {env_id}] Mesh generation completed successfully.")

    finally:
        os.chdir(original_cwd)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Generate and partition a SOD2D mesh for an airfoil")
    parser.add_argument("--env-id", type=int, default=0)
    parser.add_argument("--airfoil-file", default="airfoil")
    parser.add_argument("--work-dir", default=None, help="Default: runs/env_<env-id>")
    parser.add_argument("--config", required=True, help="e.g. flow_config_ogrd.json or flow_config_cgrd.json")
    args = parser.parse_args()

    work_dir = args.work_dir or f"runs/env_{args.env_id}"

    sod2d_airfoil(
        env_id=args.env_id,
        airfoil_file=args.airfoil_file,
        work_dir=work_dir,
        config_file=args.config,
    )
