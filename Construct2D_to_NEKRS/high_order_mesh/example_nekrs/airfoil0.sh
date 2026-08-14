#!/bin/bash
# Submits the actual NekRS simulation (step 4 of the pipeline -- see
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
#SBATCH --job-name=nekrs_ddes

### Run from this directory, and put SLURM's own log files here
#SBATCH -D .
#SBATCH --output=out.o
#SBATCH --error=error.e

### How much hardware to request. NekRS needs exactly one GPU here
### ("--gres=gpu:1"); ntasks-per-node * cpus-per-task must equal 80,
### the number of CPU cores per node on this cluster's GPU partition.
##SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=20
#SBATCH --gres=gpu:1

### Which queue to submit to, and for how long. Two options below: the
### debug queue (fast to get a slot, but capped at 2 hours -- active by
### default, good for testing that a case actually runs) and the normal
### queue (up to 2 days, for a real production run). To switch, comment
### out the acc_debug pair and uncomment the acc_bsccase pair instead.
##SBATCH --qos=acc_bsccase
#SBATCH --time=2-00:00:00
#SBATCH --qos=acc_debug
#SBATCH --time=02:00:00
#SBATCH --account=bsc21

### Only needed if this job must wait for an earlier one to finish first
### (e.g. a mesh-generation job) -- uncomment and fill in that job's own
### id (shown by `squeue` or when you submitted it) if you need this:
##SBATCH --dependency=afterany:<job-id-to-wait-for>

### Load a working compiler/MPI/GPU software environment. You'd expect a
### module set named after NekRS itself, but that one doesn't provide a
### working environment on this cluster -- this loads the same module set
### the SOD2D pipeline uses instead (compiler, MPI, CUDA: everything
### NekRS needs too):
#module getdefault nekrs
module purge
module load bsc/1.0 nvidia-hpc-sdk/24.3 hdf5/1.14.1-2-nvidia-nvhpcx \
    mkl/2025.2 python/3.12.1-gcc cmake/3.30.5

### Point at the NekRS install to use, and load its own helper functions
export NEKRS_HOME=$HOME/.local/nekrs_v24_wf
source $NEKRS_HOME/bin/nrsqsub_utils

### Finally, actually run the simulation. "naca" here is the case name --
### NekRS looks for naca.par/naca.udf/naca.oudf/naca.usr/naca.re2 in the
### current directory, all the files this pipeline built/copied in.
mpirun -np 1 $NEKRS_HOME/bin/nekrs --setup naca
