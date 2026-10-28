#!/bin/bash
# submit from the run directory (it holds tree/ and uc/):  sbatch uc/scripts/mn5_all.sbatch
#SBATCH --account=bsc08
#SBATCH --qos=gp_debug
#SBATCH --job-name=ut_all
#SBATCH --output=ut_all.%j.log
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=112
#SBATCH --time=02:00:00
bash "$SLURM_SUBMIT_DIR/uc/scripts/mn5_all.sh"
