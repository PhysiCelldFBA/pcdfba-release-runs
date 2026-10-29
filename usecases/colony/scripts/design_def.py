"""Write runs_def.csv and overrides_def.csv: the colony experiments that exist as runs (the pre-grown study).
Values are absolute. usage: design_def.py"""
import csv, os
UC = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); P = "/PhysiCell_settings"
runs, ov = [], []
run = lambda rid, eid, seed=50: runs.append((rid, eid, seed)); put = lambda s, x, v: ov.append((s, x, v))
CD = P + "/cell_definitions/cell_definition/phenotype/intracellular"
SBML = CD + "/settings/sbml_filename"; LEDGER = ("+" + CD + "/settings", "<dfba_conservation_ledger>true</dfba_conservation_ledger>")
# pre-grown colony: R0 200 um, 10 min, snapshots every 2 min, necrosis x10 (shipped 1.67e-5)
PG = ["COL-PG", "COL-DT", "COL-DX", "COL-NEC", "COL-R0"]
for e in PG:
    put(e, SBML, "./config/shipped/Ecoli_core.xml"); put(e, *LEDGER)
    put(e, P + "/overall/max_time", 10); put(e, P + "/save/full_data/interval", 2); put(e, P + "/save/SVG/enable", "false")
    put(e, P + "/user_parameters/colony_radius", 200)
    put(e, CD + "/death_model/death_rate_increase", 1.67e-5 if e == "COL-NEC" else 1.67e-4)
for rid, seed in (("COL-PG_s50", 50), ("COL-PG_s51", 51), ("COL-PG_s52", 52), ("COL-PG_s50_rep", 50)): run(rid, "COL-PG", seed)
for dt in (0.0005, 0.002, 0.004, 0.01):
    rid = f"COL-DT_dt{dt}_s50"; run(rid, "COL-DT"); mech = max(dt, 0.001)
    put(rid, P + "/overall/dt_diffusion", dt); put(rid, CD + "/settings/intracellular_dt", dt)
    put(rid, P + "/overall/dt_mechanics", mech); put(rid, P + "/overall/dt_phenotype", mech)
    put(rid, P + "/cell_definitions/cell_definition/phenotype/mechanics/cell_cell_repulsion_strength", round(1000 * 0.001 / mech, 6))
# finer dt, shortened to fit the 2 h debug queue (cost ~ 1/dt): snapshots every minute
for dt, tmax in ((0.00025, 6), (0.0001, 2)):
    rid = f"COL-DT_dt{dt}_s50"; run(rid, "COL-DT")
    put(rid, P + "/overall/max_time", tmax); put(rid, P + "/save/full_data/interval", 1)
    put(rid, P + "/overall/dt_diffusion", dt); put(rid, CD + "/settings/intracellular_dt", dt)
    put(rid, P + "/overall/dt_mechanics", 0.001); put(rid, P + "/overall/dt_phenotype", 0.001)
    put(rid, P + "/cell_definitions/cell_definition/phenotype/mechanics/cell_cell_repulsion_strength", 1000.0)
# proposed (not yet run): voxel size with the cap inactive, x1 necrosis to 40 min, and the small colony
for dx in (1, 4):
    rid = f"COL-DX_dx{dx}_dt0.0005_s50"; run(rid, "COL-DX")
    put(rid, P + "/overall/max_time", 6); put(rid, P + "/domain/dx", dx); put(rid, P + "/domain/dy", dx)
    put(rid, P + "/overall/dt_diffusion", 0.0005); put(rid, CD + "/settings/intracellular_dt", 0.0005)
    put(rid, P + "/overall/dt_mechanics", 0.001); put(rid, P + "/overall/dt_phenotype", 0.001)
    put(rid, P + "/cell_definitions/cell_definition/phenotype/mechanics/cell_cell_repulsion_strength", 1000.0)
run("COL-NEC_x1_t40_s50", "COL-NEC"); put("COL-NEC_x1_t40_s50", P + "/overall/max_time", 40); put("COL-NEC_x1_t40_s50", P + "/save/full_data/interval", 5)
rid = "COL-R0_R100_dt0.0005_s50"; run(rid, "COL-R0")
put(rid, P + "/user_parameters/colony_radius", 100); put(rid, P + "/overall/max_time", 6)
put(rid, P + "/overall/dt_diffusion", 0.0005); put(rid, CD + "/settings/intracellular_dt", 0.0005)
put(rid, P + "/overall/dt_mechanics", 0.001); put(rid, P + "/overall/dt_phenotype", 0.001)
put(rid, P + "/cell_definitions/cell_definition/phenotype/mechanics/cell_cell_repulsion_strength", 1000.0)
for dx in (1, 4):
    rid = f"COL-DX_dx{dx}_s50"; run(rid, "COL-DX"); put(rid, P + "/domain/dx", dx); put(rid, P + "/domain/dy", dx)
run("COL-NEC_x1_s50", "COL-NEC")
with open(f"{UC}/runs_def.csv", "w", newline="") as f: w = csv.writer(f); w.writerow(["run_id", "experiment_id", "seed"]); w.writerows(runs)
with open(f"{UC}/overrides_def.csv", "w", newline="") as f: w = csv.writer(f); w.writerow(["scope", "xpath", "value"]); w.writerows(ov)
print(len(runs), "runs,", len(ov), "override rows")
