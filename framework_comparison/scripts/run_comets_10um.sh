#!/bin/bash
# usage: run_comets_10um.sh <dir> ; env (MOTION, TIMESTEP_MIN, T_END_H, LOG_MIN, FLUX_H, NTHREADS, OUT_TAG, ...) passed through
export COMETS_HOME=/home/ohayoun/benchmark_tools/comets/comets_linux/comets_2.12.5/
export PATH=$PATH:$COMETS_HOME
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$COMETS_HOME/lib/or-tools
export PYTHONPATH=$PYTHONPATH:$COMETS_HOME/lib/cometspy-master/
source /home/ohayoun/benchmark_tools/venv/bin/activate
export JAVA_HOME=/usr
cd "$1" && for f in cf_ref_box_biomass.json C_beijerinckii_model_iCB925_trimmed.xml M_barkeri_model_iMG746_trimmed.xml C_beijerinckii_model_iCB925_trimmed_nocaps.xml M_barkeri_model_iMG746_trimmed_nocaps.xml; do ln -sf /home/ohayoun/benchmark_tools/cf_test/$f .; done
python comets_cfref_10um.py > run.log 2>&1
