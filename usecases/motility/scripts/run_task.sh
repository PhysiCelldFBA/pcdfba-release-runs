#!/bin/bash
# one motility run from the use-case directory (cwd). usage: run_task.sh RUN_ID    refuses a non-empty runs/RUN_ID
n=$1; [ -n "$(ls -A runs/$n 2>/dev/null)" ] && { echo "REFUSING $n"; exit 1; }
mkdir -p runs/$n
( echo "start=$(date -Is) binary=$(sha256sum bin/ecoli-dfba-motility | cut -c1-16) host=$(hostname)"; ./bin/ecoli-dfba-motility config/$n.xml; rc=$?; echo "end=$(date -Is) rc=$rc" ) > runs/$n/run.log 2>&1
