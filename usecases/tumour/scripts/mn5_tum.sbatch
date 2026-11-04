#!/bin/bash
# submit from the run directory (it holds tree/ and uc/):
#   sbatch --nodes=N --ntasks=2N --cpus-per-task=56 --qos=gp_bscls --time=HH:MM:SS --job-name=tum_x --export=ALL,RUNS='id1 id2' uc/scripts/mn5_tum.sbatch
#SBATCH --account=bsc08
#SBATCH --cpus-per-task=56
bash "$SLURM_SUBMIT_DIR/uc/scripts/mn5_tum.sh"
