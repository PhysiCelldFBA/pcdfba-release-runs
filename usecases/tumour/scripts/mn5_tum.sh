#!/bin/bash
# build the tumour binary from the exported pinned code (tree/COMMIT) if uc/bin is empty, then run RUNS ("id1 id2"), each as a step with the threads of its config
B=${PCDFBA_MN5_DIR:-${SLURM_SUBMIT_DIR:?submit from the run directory (it holds tree/ and uc/)}}
module purge >/dev/null 2>&1; module load gcc/13.2.0 >/dev/null 2>&1
if [ ! -x $B/uc/bin/cancer-core-metabolism ]; then
  cd $B/tree || exit 1
  make -j56 PHYSICELL_CPP=g++ PROGRAM_NAME=cancer-core-metabolism > $B/build.log 2>&1 || { echo BUILD FAILED; tail -20 $B/build.log; exit 1; }
  mkdir -p $B/uc/bin; cp cancer-core-metabolism $B/uc/bin/cancer-core-metabolism
fi
cd $B/uc; chmod +x scripts/run_task.sh bin/cancer-core-metabolism
echo "binary $(sha256sum bin/cancer-core-metabolism | cut -c1-16) commit $(cat $B/tree/COMMIT)"
for r in $RUNS; do srun -N1 -n1 -c${SLURM_CPUS_PER_TASK:-56} --exact --job-name=$r scripts/run_task.sh $r & done
wait; echo "all done $(date -Is)"
