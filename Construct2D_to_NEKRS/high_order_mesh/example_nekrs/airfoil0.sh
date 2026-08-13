#!/bin/bash
### Job name on queue
#SBATCH --job-name=nekrs_ddes

### Output and error files directory
#SBATCH -D .

### Output and error files
#SBATCH --output=out.o
#SBATCH --error=error.e

### Run configuration
### Rule: {ntasks-per-node} \times {cpus-per-task} = 80
##SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=20
#SBATCH --gres=gpu:1

### Queue and account
##SBATCH --qos=acc_bsccase
#SBATCH --time=2-00:00:00
#SBATCH --qos=acc_debug
#SBATCH --time=02:00:00
#SBATCH --account=bsc21
##sbatch --dependency=afterany:9878054

### MN% modules
#module getdefault nekrs
module getdefault sod2d

export NEKRS_HOME=$HOME/.local/nekrs_v24_wf
source $NEKRS_HOME/bin/nrsqsub_utils

mpirun -np 1 $NEKRS_HOME/bin/nekrs --setup naca 
