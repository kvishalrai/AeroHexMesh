#!/bin/bash
### Job name on queue
#SBATCH --job-name=sod2d_ale

### Output and error files directory
#SBATCH -D .

### Output and error files
#SBATCH --output=out.o
#SBATCH --error=error.e

### Run configuration
### Rule: {ntasks-per-node} \times {cpus-per-task} = 80
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=20
#SBATCH --gres=gpu:4

### Queue and account
##SBATCH --qos=acc_bsccase
##SBATCH --time=1-00:00:00
#SBATCH --qos=acc_debug
#SBATCH --time=00:20:00
#SBATCH --account=bsc21
### MN% modules
module purge
module load nvidia-hpc-sdk/24.3 hdf5/1.14.1-2-nvidia-nvhpcx

SOD2D_SRC_DIR=../CFD_code/sod2d_gitlab/build_gpu/src/app_sod2d/

# Works around a UCX rendezvous-protocol segfault (ucp_rkey_pack_memh) seen
# under certain mesh-size/partition-count combinations at num_partitions>=3;
# not a real mesh/GEMPA bug. See project_ucx_partitioning_fix.md memory.
export UCX_TLS=rc,cuda_copy,cuda_ipc,sm,self

# A large mesh HDF5 read (the Parallel_data group) can segfault inside Open
# MPI's native ompio component's collective read path
# (mca_common_ompio_file_read_at_all -> H5D__mpio_select_read); force ROMIO
# instead. This is a separate crash from the UCX_TLS one above (MPI-IO, not
# halo-exchange point-to-point) -- see the same memory note.
mpirun -np 4 --mca io ^ompio --map-by ppr:4:node:PE=20 --report-bindings ./mn5_bind.sh $SOD2D_SRC_DIR/sod2d MeshElasticitySolver
