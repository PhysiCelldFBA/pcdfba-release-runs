#!/bin/bash
# Run the dfba_unit_test checks (from the pinned code) on the runs of this use case. The checks expect the old folder names, so a
# scratch directory of symlinks maps runs/<run_id> to them. usage: scripts/run_checks.sh   (from the use-case directory; PYTHON=... for cobra)
PY=${PYTHON:-/opt/anaconda3/bin/python}; UC=$(pwd); CK=$UC/checks; CODE=$UC/../../code/sample_projects_intracellular/fba/dfba_unit_test/scripts
rm -rf $CK; mkdir -p $CK; rc=0
for r in $UC/runs/UT-DXDT-A_*; do n=$(basename $r); ln -s $r $CK/output_dxdt_inactive_${n#UT-DXDT-A_}; done
for r in $UC/runs/UT-DXDT-B_*; do n=$(basename $r); ln -s $r $CK/output_dxdt_capped_${n#UT-DXDT-B_}; done
for c in n1 n4 n8; do ln -s $UC/runs/UT-VOX-${c^^}_s50 $CK/output_voxel_$c; done; ln -s $UC/runs/UT-VOX-MIXED_s50 $CK/output_voxel_mixed
for r in $UC/runs/UT-ORDER_*; do n=$(basename $r); ln -s $r $CK/output_order_${n#UT-ORDER_}; done
echo "### dx/dt (A-D)"; $PY $CODE/check_dxdt.py --dir $CK --model $UC/config/shipped/Ecoli_core.xml || rc=1
echo "### elemental balance (E)"; $PY $CODE/elemental_check.py $UC/config/shipped/Ecoli_core.xml
echo "### voxel sharing"; for c in n1 n4 n8 mixed; do (cd $CK && $PY $CODE/check_voxel_nonnegative.py output_voxel_$c) || rc=1; done
(cd $CK && $PY $CODE/check_voxel_draw_equal.py output_voxel_n4 output_voxel_mixed) || rc=1
echo "### update order"; (cd $CK && $PY $CODE/check_update_order.py) || rc=1
exit $rc
