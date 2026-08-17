# Cluster / environment setup (unified reference)

The `AEROHEXMESH_*_MODULE_SETUP` env-var pattern recurs in three
pipelines with the same shape but different variable names — this page
collects them in one place, plus the recurring MPI-crash workarounds.

## The `AEROHEXMESH_*_MODULE_SETUP` pattern

Each pipeline's own orchestrator shells out to tools (Gmsh, MPI, HDF5,
pyHyp, …) that need to already be on `PATH`. Rather than hardcoding a
`module load` line with no way to override it, each pipeline runs
whatever shell command is in a dedicated environment variable right
before the step that needs it — set your own, or rely on the MN5
(MareNostrum 5) default baked in as a fallback.

| Pipeline | Env var(s) | MN5 default |
|---|---|---|
| `Construct2D_to_SOD2D/linear_mesh/` | `AEROHEXMESH_MODULE_SETUP`, `AEROHEXMESH_PYTHON_MODULE_SETUP` | `module purge && module load bsc/1.0 nvidia-hpc-sdk/24.3 hdf5/1.14.1-2-nvidia-nvhpcx mkl/2025.2 python/3.12.1-gcc cmake/3.30.5`, then `module unload python && module load anaconda` for the two steps needing `numpy`/`h5py` |
| `Construct2D_to_NEKRS/linear_mesh/` | same pattern, see its own README | Nek5000-toolchain equivalent |
| `pyHyp/` | `AEROHEXMESH_PYHYP_MODULE_SETUP` | conda-based (pyHyp isn't on PyPI — built from source; see `pyHyp/README.md`'s own build recipe) |

Set either half to `""` for a no-op (e.g. if you've already activated
everything yourself — a virtualenv, `conda activate`, …) before running
the pipeline.

## SLURM job scripts (GPU steps)

Every SOD2D GPU step (`MeshElasticitySolver`, production solves, across
`Construct2D_to_SOD2D/high_order_mesh/`, `etaGrid_to_SOD2D/*_run/`,
`pyHyp_to_SOD2D/`) uses the same shape:

```bash
module purge
module load nvidia-hpc-sdk/24.3 hdf5/1.14.1-2-nvidia-nvhpcx
export UCX_TLS=rc,cuda_copy,cuda_ipc,sm,self
mpirun -np 4 --mca io ^ompio --map-by ppr:4:node:PE=20 --report-bindings \
    ./mn5_bind.sh $SOD2D_SRC_DIR/sod2d <SolverType>
```

- **`mn5_bind.sh`**: pins each of the node's 4 MPI ranks to its own GPU,
  CPU NUMA domain, and network interface card
  (`CUDA_VISIBLE_DEVICES`/`numactl --membind`/`UCX_NET_DEVICES` per
  `OMPI_COMM_WORLD_LOCAL_RANK`) — without it, ranks contend for whichever
  GPU/memory/NIC happens to be nearest, silently slowing every transfer.
- **`UCX_TLS=rc,cuda_copy,cuda_ipc,sm,self`**: works around a UCX
  rendezvous-protocol segfault (`ucp_rkey_pack_memh`) seen during
  halo-exchange at `num_partitions ≥ 3` on a large mesh — forcing this
  specific transport set avoids it. Not a mesh/partitioning bug.
- **`--mca io ^ompio`** (equivalently `export OMPI_MCA_io=romio321` for
  the CPU-only partitioner step): works around a *separate* crash —
  reading a large mesh HDF5 file's `Parallel_data` group can segfault
  inside Open MPI's native `ompio` component during the initial parallel
  collective read; forcing the older ROMIO component avoids it. A
  different code path (file I/O, not point-to-point messaging) from the
  UCX issue above — fixing one does not fix the other, both are needed
  together.
- **`--bind-to none` / `SLURM_CPU_BIND=none`**: NVHPC's own MPI
  implementation doesn't support `srun` directly — must be launched via
  `mpirun --bind-to none`, with `SLURM_CPU_BIND=none` set so SLURM's own
  binding doesn't fight `mpirun`'s.

## Login-node resource limits

Both foreground and `nohup`-backgrounded Python/Gmsh/SOD2D processes get
killed (exit 137, SIGKILL) on the login node for anything beyond
trivially small computations — **always submit real work via SLURM**,
even for something that "should" be quick. Use `python3 -u` for
unbuffered output when backgrounding anything — Python fully buffers
stdout when not attached to a TTY, so a killed background process can
leave an empty log file even after doing real work.

## `--mem` is not a valid SLURM request on this cluster

`sbatch --mem=...` is rejected outright
(`You cannot submit a job requesting memory parameters, memory is
automatically set for each asked cpu`). Request more memory indirectly,
via `--cpus-per-task` (2G/core default on the general-purpose queue).
