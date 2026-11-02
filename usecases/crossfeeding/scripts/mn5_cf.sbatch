#!/bin/bash
# submit from the run directory (it holds tree/ and uc/):
#   sbatch --nodes=N --ntasks=4N --qos=gp_bscls --time=HH:MM:SS --job-name=cf_x --export=ALL,RUNS='id1 id2' uc/scripts/mn5_cf.sbatch
# Builds the binary on the first node, then runs each run as a 28-thread step (4 per node).
#SBATCH --account=bsc08
#SBATCH --cpus-per-task=28
bash "$SLURM_SUBMIT_DIR/uc/scripts/mn5_cf.sh"
