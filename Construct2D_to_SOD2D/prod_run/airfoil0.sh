#!/bin/bash
### Job name on queue
#SBATCH --job-name=str_lps

### Output and error files directory
#SBATCH -D .

### Output and error files
#SBATCH --output=out.o
#SBATCH --error=error.e

### Run configuration
### num_partitions=4: the earlier segfault at num_partitions>=3 (UCX
### ucp_rkey_pack_memh crash during the first halo exchange) was root-caused
### to UCX's default transport selection, not a hard GEMPA/partition-count
### limit -- confirmed fixed by restricting UCX_TLS below (verified clean at
### np=3 and np=4 on both the O-grid and C-grid meshes).
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=20
#SBATCH --gres=gpu:4

### Queue and account
##SBATCH --qos=acc_bsccase
##SBATCH --time=12:00:00
#SBATCH --qos=acc_debug
#SBATCH --time=02:00:00
#SBATCH --account=bsc21
#sbatch --dependency=afterany:9794090

### MN% modules
module purge
module load nvidia-hpc-sdk/24.3 hdf5/1.14.1-2-nvidia-nvhpcx

SOD2D_SRC_DIR=../CFD_code/sod2d_gitlab/build_gpu/src/app_sod2d

export UCX_TLS=rc,cuda_copy,cuda_ipc,sm,self

mpirun -np 4 --bind-to none ./mn5_bind.sh $SOD2D_SRC_DIR/sod2d BluffBodySolverIncomp
