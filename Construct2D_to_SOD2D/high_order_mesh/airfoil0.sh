#!/bin/bash
# Submits the wall-smoothing run (step 4 of the pipeline -- see
# ../README.md) to the SLURM job scheduler.
#
# If you haven't used a scheduler like SLURM before: on a shared
# supercomputer, you don't run your program directly -- you can't just
# grab a GPU and start using it, because everyone else on the cluster is
# competing for the same hardware. Instead you write a request ("how many
# nodes/GPUs, for how long") as a script like this one, hand it to SLURM
# with `sbatch airfoil0.sh`, and SLURM runs the actual command at the
# bottom once matching resources are free. The `#SBATCH` lines below are
# that request, read by SLURM itself, not by bash.
#
# This script is specific to BSC's MareNostrum 5 cluster (module names,
# queue names, account) -- adapt it if you're running somewhere else.

### Job name, as it will appear in the queue
#SBATCH --job-name=sod2d_ale

### Run from this directory, and put SLURM's own log files here
#SBATCH -D .
#SBATCH --output=out.o
#SBATCH --error=error.e

### How much hardware to request: one full node, running 4 MPI ranks (one
### per GPU). ntasks-per-node * cpus-per-task must equal 80, the number
### of CPU cores per node on this cluster's GPU partition.
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=20
#SBATCH --gres=gpu:4

### Which queue to submit to, and for how long. Two options below: the
### debug queue (fast to get a slot, but capped at 20 minutes -- active
### by default, good for testing that a run actually works) and the
### normal queue (up to 1 day, for a real production-size mesh). To
### switch, comment out the acc_debug pair and uncomment the acc_bsccase
### pair instead.
##SBATCH --qos=acc_bsccase
##SBATCH --time=1-00:00:00
#SBATCH --qos=acc_debug
#SBATCH --time=00:20:00
#SBATCH --account=bsc21

### Load a working compiler/MPI/GPU software environment (MN5 modules)
module purge
module load nvidia-hpc-sdk/24.3 hdf5/1.14.1-2-nvidia-nvhpcx

SOD2D_SRC_DIR=../../CFD_solvers/sod2d_gitlab/build_gpu/src/app_sod2d/

# Two unrelated crashes have been seen running this on a large mesh with
# 3+ partitions, and both workarounds below are needed together -- if
# you hit a segfault here, check whether one of these got dropped rather
# than assuming it's a new bug:
#
# 1. A UCX (the network transport library MPI uses here) rendezvous-
#    protocol segfault (ucp_rkey_pack_memh) during halo-exchange (the
#    communication between neighboring mesh partitions). Not a real
#    mesh/partitioning bug -- forcing a different set of transports
#    avoids it. See the project_ucx_partitioning_fix.md memory for the
#    full diagnosis.
export UCX_TLS=rc,cuda_copy,cuda_ipc,sm,self

# 2. A *separate* crash: reading the large mesh HDF5 file's Parallel_data
#    group can segfault inside Open MPI's native "ompio" component
#    (mca_common_ompio_file_read_at_all -> H5D__mpio_select_read) during
#    the initial parallel collective read. Forcing the older ROMIO
#    component for MPI-IO instead avoids it. This is a different code
#    path (file I/O, not point-to-point messaging) from the UCX crash
#    above -- fixing one does not fix the other.
mpirun -np 4 --mca io ^ompio --map-by ppr:4:node:PE=20 --report-bindings ./mn5_bind.sh $SOD2D_SRC_DIR/sod2d MeshElasticitySolver
