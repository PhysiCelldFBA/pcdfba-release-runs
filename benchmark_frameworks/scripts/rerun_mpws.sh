#!/bin/bash
# Rerun one of the four benchmark runs on mpws (bsclife006) with the environment of env/.
# usage, from benchmark_frameworks/ of a checkout on mpws (Steps 3 and 4 of the top README done, so that the shipped
# models and cells_CF3-REF_s50.csv exist):
#     bash scripts/rerun_mpws.sh <run_id> runs/<run_id>
# run_id: BM-ACE-COMETS | BM-ACE-BACARENA | BM-CF-COMETS | BM-CF-BACARENA (runs.csv). The output folder must not exist.
# Stages the inputs, checks each against config/SHA256SUMS (the files of the recorded runs), runs the tool with the settings of runs.csv,
# then the extraction script; the small outputs read by make_benchmark.py are then copied by hand to data/<folder>/.
# BENCHMARK_TOOLS (default /home/ohayoun/benchmark_tools) holds COMETS 2.12.5, the Python venv and R_libs (env/).
set -euo pipefail
RUN=$1; OUT=$2
B=${BENCHMARK_TOOLS:-/home/ohayoun/benchmark_tools}
C=$B/comets/comets_linux/comets_2.12.5
HERE=$(cd "$(dirname "$0")/.." && pwd); REPO=$(dirname "$HERE"); S=$HERE/scripts
[ -e "$OUT" ] && { echo "$OUT exists; not overwriting" >&2; exit 1; }
mkdir -p "$OUT"; OUT=$(cd "$OUT" && pwd); cd "$OUT"

stage() { cp "$1" "${2:-.}"; }
comets_env() {   # as run_comets_10um.sh / the COMETS ACE launch
  export COMETS_HOME=$C/ PATH=$PATH:$C/ LD_LIBRARY_PATH=${LD_LIBRARY_PATH:-}:$C/lib/or-tools
  export PYTHONPATH=${PYTHONPATH:-}:$C/lib/cometspy-master/ JAVA_HOME=/usr
  set +u; source "$B/venv/bin/activate"; set -u
}
verify() {   # every staged file must be listed in config/SHA256SUMS with the same hash
  for f in *; do grep -q "  $f\$" "$HERE/config/SHA256SUMS" || { echo "$f not in config/SHA256SUMS" >&2; exit 1; }; done
  grep -E "  ($(ls | paste -sd'|'))\$" "$HERE/config/SHA256SUMS" | sha256sum -c --quiet - || { echo "input checksum mismatch" >&2; exit 1; }
}
NOCAPS=$REPO/usecases/unit_test/models_full

case $RUN in
  BM-ACE-COMETS)
    stage "$S/comets_ecoli_aceref_v2.py"; stage "$REPO/usecases/acetate/config/shipped/Ecoli_core.xml"
    verify
    comets_env; python comets_ecoli_aceref_v2.py > run.log 2>&1 ;;
  BM-ACE-BACARENA)
    stage "$S/run_bacarena_aceref_v6.R"
    verify
    TSTEP=0.01 GRID_N=17 DEATHRATE=0 OUT_SFX=_dt0.01_n17_d0 Rscript run_bacarena_aceref_v6.R > run_dt0.01_n17_d0.log 2>&1 ;;
  BM-CF-COMETS)
    stage "$S/comets_cfref_10um.py"; stage "$S/extract_comets_cf_qual.py"
    stage "$REPO/usecases/crossfeeding/config/cells_CF3-REF_s50.csv"
    stage "$NOCAPS/C_beijerinckii_model_iCB925_trimmed_nocaps.xml"; stage "$NOCAPS/M_barkeri_model_iMG746_trimmed_nocaps.xml"
    verify
    comets_env
    SCENARIO=fig6 MOTION=convection TIMESTEP_MIN=0.25 OUT_TAG=fig6_conv python comets_cfref_10um.py > run.log 2>&1
    python extract_comets_cf_qual.py . fig6_conv > extract.log 2>&1 ;;
  BM-CF-BACARENA)
    stage "$S/run_bacarena_fig6_80x80.R"; stage "$S/extract_bacarena_fig6_qual.R"
    for f in cb_modelorg.rds mb_modelorg.rds nocaps_bounds_cb.csv nocaps_bounds_mb.csv; do stage "$HERE/config/$f"; done
    verify
    BOUNDARY=closed INIT_POOL=0 SBML_BOUNDS=1 OUT_SFX=_closed_nopool Rscript run_bacarena_fig6_80x80.R > run.log 2>&1
    Rscript extract_bacarena_fig6_qual.R > extract.log 2>&1 ;;
  *) echo "unknown run_id $RUN" >&2; exit 1 ;;
esac
echo "$RUN done in $OUT"; ls "$OUT"
