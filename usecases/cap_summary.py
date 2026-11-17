"""Cap number R and capped fractions of the reference runs (Supplementary Methods S7.7, Table cap_summary).
usage: /opt/anaconda3/bin/python usecases/cap_summary.py          (run from the repository root)
Reads the per-snapshot tables written by colony/scripts/cap_check.py (colony, motility; label cap_check_max_<run>.csv) and by
crossfeeding/scripts/analyse_cf.py (cf3nc50_cap.csv, main-figure consortium run). For every substrate, snapshots with t > 0 only
(the t = 0 snapshot has no recorded bound): the largest R over cells and snapshots, the largest median R, the largest fraction of cells
whose recorded bound is below 0.99 of the Michaelis-Menten bound, and the lowest median ratio.
The 'capped' fraction can be above zero with R << 1 when the concentration of the snapshot differs from the one the bound was set from (C. beijerinckii, glucose)."""
import pandas as pd
pd.set_option("display.width", 200)
def show(label, f, by):
    d = pd.read_csv(f); d = d[d.t > 0]
    g = d.groupby(by).agg(R_max=("max_R", "max"), R_med_max=("median_R", "max"), capped_max_pct=("frac_capped_pct", "max"), median_ratio_min=("median_ratio", "min"))
    print("==", label); print(g.round(3).to_string())
show("colony COL-REF_s50 (dx 2, dt 0.0005, 600 min)", "usecases/colony/results/cap_check_max_COL-REF_s50.csv", "substrate")
show("motility high supply MOT4-HIGH_s50", "usecases/motility/results/cap_check_max_MOT4-HIGH_s50.csv", "substrate")
show("motility low supply MOT3-LOW_s50", "usecases/motility/results/cap_check_max_MOT3-LOW_s50.csv", "substrate")
show("consortium CF3-NC50_s50 (main figure)", "usecases/crossfeeding/results/cf3nc50_cap.csv", ["species", "substrate"])
