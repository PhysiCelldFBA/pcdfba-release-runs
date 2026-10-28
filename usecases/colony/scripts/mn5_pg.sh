#!/bin/bash
# build the colony binary from the exported pinned code (tree/COMMIT) on the first node, then run the pre-grown runs, one node each
B=${PCDFBA_MN5_DIR:-${SLURM_SUBMIT_DIR:?submit from the run directory or set PCDFBA_MN5_DIR}}   # holds tree/ and uc/
module purge >/dev/null 2>&1; module load gcc/13.2.0 >/dev/null 2>&1
cd $B/tree || exit 1
make -j112 PHYSICELL_CPP=g++ PROGRAM_NAME=bacterial-colony > $B/build.log 2>&1 || { echo BUILD FAILED; tail -20 $B/build.log; exit 1; }
cp bacterial-colony $B/uc/bin/bacterial-colony; cd $B/uc; chmod +x scripts/run_task.sh bin/bacterial-colony
echo "binary $(sha256sum bin/bacterial-colony | cut -c1-16) commit $(cat $B/tree/COMMIT)"
for r in $(tail -n +2 runs_def.csv | grep -E '^COL-(PG|DT|DX|NEC)' | cut -d, -f1); do
  srun -N1 -n1 -c112 --exclusive --job-name=$r scripts/run_task.sh $r &
done
wait; echo "all done $(date -Is)"
