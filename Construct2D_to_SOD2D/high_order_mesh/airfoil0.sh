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
##sbatch --dependency=afterany:9878054

### MN% modules
module purge
module load nvidia-hpc-sdk/24.3 hdf5/1.14.1-2-nvidia-nvhpcx

SOD2D_SRC_DIR=../CFD_code/sod2d_gitlab/build_gpu/src/app_sod2d/
#SOD2D_SRC_DIR=/gpfs/scratch/bsc21/bsc021712/1.sod2d/3.naca/5.RE200K/sod2d_gitlab/build_p4/src/app_sod2d/

mpirun -np 2 --map-by ppr:4:node:PE=20 --report-bindings ./mn5_bind.sh $SOD2D_SRC_DIR/sod2d MeshElasticitySolver
