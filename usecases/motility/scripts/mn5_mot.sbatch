#!/bin/bash
# submit from the run directory (it holds tree/ and uc/):
#   sbatch --nodes=1 --ntasks=8 --cpus-per-task=14 --qos=gp_debug --time=01:00:00 --job-name=mot --export=ALL,RUNS='id1 id2 ...' uc/scripts/mn5_mot.sbatch
# Builds the binary on the node if uc/bin is empty, then runs every run of RUNS as a step with 14 cpus (the configs use 6 threads).
#SBATCH --account=bsc08
#SBATCH --cpus-per-task=14
bash "$SLURM_SUBMIT_DIR/uc/scripts/mn5_mot.sh"
