#!/bin/bash
# build the unit-test binary from the exported pinned code (tree/COMMIT), then run every planned run: one-thread runs 96 at a time, then the multi-thread runs one by one
B=${PCDFBA_MN5_DIR:-${SLURM_SUBMIT_DIR:?submit from the run directory or set PCDFBA_MN5_DIR}}   # holds tree/ and uc/
module purge >/dev/null 2>&1; module load gcc/13.2.0 >/dev/null 2>&1
cd $B/tree || exit 1
make -j112 PHYSICELL_CPP=g++ PROGRAM_NAME=dfba-unit-test > $B/build.log 2>&1 || { echo BUILD FAILED; tail -20 $B/build.log; exit 1; }
cp dfba-unit-test $B/uc/bin/dfba-unit-test; cd $B/uc
chmod +x scripts/run_task.sh bin/dfba-unit-test
all=$(cat runs_def.csv runs.csv 2>/dev/null | tail -n +2 | grep -v '^run_id' | cut -d, -f1 | sort -u)
one=""; many=""
for r in $all; do [ -f config/$r.xml ] || continue; t=$(grep -o '<omp_num_threads>[0-9]*' config/$r.xml | grep -o '[0-9]*$'); if [ "${t:-1}" = "1" ]; then one="$one $r"; else many="$many $r"; fi; done
echo "one-thread: $(echo $one | wc -w), multi-thread: $(echo $many | wc -w)"
printf '%s\n' $one | xargs -P 96 -I{} scripts/run_task.sh {}
printf '%s\n' $many | xargs -P 1 -I{} scripts/run_task.sh {}
echo "all done $(date -Is)"
