#!/bin/bash
# submit from the run directory (it holds tree/ and uc/):  sbatch uc/scripts/mn5_pg.sbatch
#SBATCH --account=bsc08
#SBATCH --qos=gp_debug
#SBATCH --job-name=col_pg_dt
#SBATCH --output=col_pg.%j.log
#SBATCH --nodes=11
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=112
#SBATCH --time=02:00:00
bash "$SLURM_SUBMIT_DIR/uc/scripts/mn5_pg.sh"
