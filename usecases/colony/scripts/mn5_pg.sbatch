#!/bin/bash
#SBATCH --account=bsc08
#SBATCH --qos=gp_debug
#SBATCH --chdir=/gpfs/scratch/bsc08/bsc008059/pcdfba_review/colony_pg_mn5
#SBATCH --job-name=col_pg_dt
#SBATCH --output=/gpfs/scratch/bsc08/bsc008059/pcdfba_review/colony_pg_mn5/col_pg.%j.log
#SBATCH --nodes=11
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=112
#SBATCH --time=02:00:00
/gpfs/scratch/bsc08/bsc008059/pcdfba_review/colony_pg_mn5/uc/scripts/mn5_pg.sh
