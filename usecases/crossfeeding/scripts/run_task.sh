#!/bin/bash
# one cross-feeding run from the use-case directory (cwd). usage: run_task.sh RUN_ID   refuses a non-empty runs/RUN_ID; threads from the config
n=$1; [ -n "$(ls -A runs/$n 2>/dev/null)" ] && { echo "REFUSING $n"; exit 1; }
mkdir -p runs/$n
th=$(grep -o '<omp_num_threads>[0-9]*' config/$n.xml | grep -o '[0-9]*$')
( echo "start=$(date -Is) binary=$(sha256sum bin/crossfeeding | cut -c1-16) host=$(hostname)"; OMP_NUM_THREADS=${th:-28} ./bin/crossfeeding config/$n.xml; rc=$?; echo "end=$(date -Is) rc=$rc" ) > runs/$n/run.log 2>&1
