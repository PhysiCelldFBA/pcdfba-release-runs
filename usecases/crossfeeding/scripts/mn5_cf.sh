#!/bin/bash
B=${SLURM_SUBMIT_DIR:?submit from the run directory}   # holds tree/ and uc/
module purge >/dev/null 2>&1; module load gcc/13.2.0 >/dev/null 2>&1
if [ ! -x $B/uc/bin/crossfeeding ]; then
  cd $B/tree || exit 1
  make -j56 PHYSICELL_CPP=g++ PROGRAM_NAME=crossfeeding > $B/build.log 2>&1 || { echo BUILD FAILED; tail -20 $B/build.log; exit 1; }
  mkdir -p $B/uc/bin; cp crossfeeding $B/uc/bin/crossfeeding
fi
cd $B/uc; chmod +x scripts/run_task.sh bin/crossfeeding
echo "binary $(sha256sum bin/crossfeeding | cut -c1-16) commit $(cat $B/tree/COMMIT)"
for r in $RUNS; do srun -N1 -n1 -c28 --exact --job-name=$r scripts/run_task.sh $r & done
wait; echo "all done $(date -Is)"
