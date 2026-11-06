"""Write runs_def.csv and overrides_def.csv: the definitive motility runs on the pinned code (MOT-*).
Two variants of the metabolism-driven motility example, each from its shipped settings file (they differ in Km(glucose), glucose_threshold_low and motility_cost_at_vmax only):
  MOT-REF  motility-only: config/shipped/PhysiCell_settings_motility.xml (Km 0.5 mM, threshold 0.5 mM, cost 48.9), seeds 50, 51, 52 and a repeat of seed 50 (MOT-REF-REP)
  MOT-GR   growth recovery: config/shipped/PhysiCell_settings_motility_proliferation.xml (Km 0.02 mM, threshold 0.35 mM, cost 66.4), same seeds and repeat (MOT-GR-REP)
All runs: 600 min, dt_diffusion = intracellular_dt = 0.001 min (the harmonised configuration, PR #16), dx 5 um, 6 threads, ledger on, no SVG, snapshots every 5 min. The cell starts from config/shipped/cells.csv (one cell at x = -140 um).
usage: design_def.py"""
import csv, os
UC = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); P = "/PhysiCell_settings"
CD = P + "/cell_definitions/cell_definition/phenotype/intracellular/settings"
runs, ov = [], []
run = lambda rid, eid, seed=50: runs.append((rid, eid, seed)); put = lambda s, x, v: ov.append((s, x, v))
put("*", CD + "/sbml_filename", "./config/shipped/Ecoli_core.xml")
put("*", "+" + CD, "<dfba_conservation_ledger>true</dfba_conservation_ledger>")
put("*", P + "/initial_conditions/cell_positions/folder", "./config/shipped")
put("*", P + "/save/SVG/enable", "false")
for eid in ("MOT-REF", "MOT-GR"):
    for s in (50, 51, 52): run(f"{eid}_s{s}", eid, s)
    run(f"{eid}-REP_s50", eid, 50)
with open(f"{UC}/runs_def.csv", "w", newline="") as f: w = csv.writer(f, lineterminator="\n"); w.writerow(["run_id", "experiment_id", "seed"]); w.writerows(runs)
with open(f"{UC}/overrides_def.csv", "w", newline="") as f: w = csv.writer(f, lineterminator="\n"); w.writerow(["scope", "xpath", "value"]); w.writerows(ov)
print(len(runs), "runs,", len(ov), "override rows")
