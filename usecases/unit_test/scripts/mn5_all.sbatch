#!/bin/bash
#SBATCH --account=bsc08
#SBATCH --qos=gp_debug
#SBATCH --chdir=/gpfs/scratch/bsc08/bsc008059/pcdfba_review/unit_test_mn5
#SBATCH --job-name=ut_all
#SBATCH --output=/gpfs/scratch/bsc08/bsc008059/pcdfba_review/unit_test_mn5/ut_all.%j.log
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=112
#SBATCH --time=02:00:00
/gpfs/scratch/bsc08/bsc008059/pcdfba_review/unit_test_mn5/uc/scripts/mn5_all.sh
