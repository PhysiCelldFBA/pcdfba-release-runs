"""Write runs_def.csv and overrides_def.csv: the colony experiments that exist as runs (the pre-grown study).
Values are absolute. usage: design_def.py"""
import csv, os
UC = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); P = "/PhysiCell_settings"
runs, ov = [], []
run = lambda rid, eid, seed=50: runs.append((rid, eid, seed)); put = lambda s, x, v: ov.append((s, x, v))
CD = P + "/cell_definitions/cell_definition/phenotype/intracellular"
SBML = CD + "/settings/sbml_filename"; LEDGER = ("+" + CD + "/settings", "<dfba_conservation_ledger>true</dfba_conservation_ledger>")
# pre-grown colony: R0 200 um, 10 min, snapshots every 2 min, necrosis x10 (shipped 1.67e-5)
PG = ["COL-PG", "COL-DT", "COL-DX", "COL-NEC", "COL-R0", "COL-QTR"]
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
def quarter_rows(scope):   # symmetry planes at x = 0 and y = 0: domain starts at 0, no Dirichlet on xmin/ymin, only the x >= 0, y >= 0 part of the initial disc
    put(scope, "+" + P + "/user_parameters", '<quarter_colony type="bool" units="" description="only the quadrant x >= 0, y >= 0 of the initial disc (quarter-domain runs)">true</quarter_colony>')
    put(scope, P + "/domain/x_min", 0); put(scope, P + "/domain/y_min", 0)
    for v in ("oxygen", "glucose", "acetate", "lactate", "CO2"):
        for face in ("xmin", "ymin"): put(scope, P + f"/microenvironment_setup/variable[@name='{v}']/Dirichlet_options/boundary_value[@ID='{face}']/@enabled", "False")
# COL-QTRG: growth from the shipped start (colony_radius 10 um, shipped dt 0.001 and necrosis) to 240 min in a +-150 um box: full colony against its quarter, 2 seeds each
put("COL-QTRG", SBML, "./config/shipped/Ecoli_core.xml"); put("COL-QTRG", P + "/overall/max_time", 240); put("COL-QTRG", P + "/save/full_data/interval", 20); put("COL-QTRG", P + "/save/SVG/enable", "false")
put("COL-QTRG", P + "/parallel/omp_num_threads", 32)
for a in ("x_min", "y_min"): put("COL-QTRG", P + "/domain/" + a, -150)
for a in ("x_max", "y_max"): put("COL-QTRG", P + "/domain/" + a, 150)
for seed in (50, 51):
    run(f"COL-QTRG_full_s{seed}", "COL-QTRG", seed); run(f"COL-QTRG_qtr_s{seed}", "COL-QTRG", seed); quarter_rows(f"COL-QTRG_qtr_s{seed}")
# COL-QTRL: the same, to 340 min (about 20k cells, colony radius about 100 um, still inside the box; 420 min would reach about 60k cells and 170 um)
put("COL-QTRL", SBML, "./config/shipped/Ecoli_core.xml"); put("COL-QTRL", P + "/overall/max_time", 340); put("COL-QTRL", P + "/save/full_data/interval", 20); put("COL-QTRL", P + "/save/SVG/enable", "false")
put("COL-QTRL", P + "/parallel/omp_num_threads", 32)
for a in ("x_min", "y_min"): put("COL-QTRL", P + "/domain/" + a, -150)
for a in ("x_max", "y_max"): put("COL-QTRL", P + "/domain/" + a, 150)
for seed in (50, 51):
    run(f"COL-QTRL_full_s{seed}", "COL-QTRL", seed); run(f"COL-QTRL_qtr_s{seed}", "COL-QTRL", seed); quarter_rows(f"COL-QTRL_qtr_s{seed}")
# COL-QTRN: quarter colony only, to 500 min in a +-250 um box (legacy full run: first dead cells between 240 and 360 min, 13.7k dead at 480 min, radius 158 um at 480 min)
put("COL-QTRN", SBML, "./config/shipped/Ecoli_core.xml"); put("COL-QTRN", P + "/overall/max_time", 500); put("COL-QTRN", P + "/save/full_data/interval", 20); put("COL-QTRN", P + "/save/SVG/enable", "false")
put("COL-QTRN", P + "/parallel/omp_num_threads", 112)
for a in ("x_min", "y_min"): put("COL-QTRN", P + "/domain/" + a, -250)
for a in ("x_max", "y_max"): put("COL-QTRN", P + "/domain/" + a, 250)
for seed in (50, 51):
    run(f"COL-QTRN_qtr_s{seed}", "COL-QTRN", seed); quarter_rows(f"COL-QTRN_qtr_s{seed}")
# COL-QTRB: as COL-QTRN but in the box of the legacy full run (+-500 um, quarter 0..500) and to 600 min: removes the box as a difference to colony_ref_s51
put("COL-QTRB", SBML, "./config/shipped/Ecoli_core.xml"); put("COL-QTRB", P + "/overall/max_time", 600); put("COL-QTRB", P + "/save/full_data/interval", 20); put("COL-QTRB", P + "/save/SVG/enable", "false")
put("COL-QTRB", P + "/parallel/omp_num_threads", 112)
for a in ("x_min", "y_min"): put("COL-QTRB", P + "/domain/" + a, -500)
for a in ("x_max", "y_max"): put("COL-QTRB", P + "/domain/" + a, 500)
for seed in (50, 51):
    run(f"COL-QTRB_qtr_s{seed}", "COL-QTRB", seed); quarter_rows(f"COL-QTRB_qtr_s{seed}")
# COL-QTR: the quadrant x >= 0, y >= 0 of the pre-grown colony (symmetry planes at x = 0 and y = 0: no-flux there, Dirichlet on the outer faces only);
# needs the binary with the quarter_colony option (code branch feature/colony-quarter-disc)
QTR = "COL-QTR_dt0.0005_s50"; run(QTR, "COL-QTR")
put("COL-QTR", "+" + P + "/user_parameters", '<quarter_colony type="bool" units="" description="only the quadrant x >= 0, y >= 0 of the initial disc (quarter-domain runs)">true</quarter_colony>')
put("COL-QTR", P + "/domain/x_min", 0); put("COL-QTR", P + "/domain/y_min", 0)
for v in ("oxygen", "glucose", "acetate", "lactate", "CO2"):
    for face in ("xmin", "ymin"): put("COL-QTR", P + f"/microenvironment_setup/variable[@name='{v}']/Dirichlet_options/boundary_value[@ID='{face}']/@enabled", "False")
put(QTR, P + "/overall/max_time", 6)
put(QTR, P + "/overall/dt_diffusion", 0.0005); put(QTR, CD + "/settings/intracellular_dt", 0.0005)
put(QTR, P + "/overall/dt_mechanics", 0.001); put(QTR, P + "/overall/dt_phenotype", 0.001)
put(QTR, P + "/cell_definitions/cell_definition/phenotype/mechanics/cell_cell_repulsion_strength", 1000.0)
for dx in (1, 4):
    rid = f"COL-DX_dx{dx}_s50"; run(rid, "COL-DX"); put(rid, P + "/domain/dx", dx); put(rid, P + "/domain/dy", dx)
run("COL-NEC_x1_s50", "COL-NEC")
with open(f"{UC}/runs_def.csv", "w", newline="") as f: w = csv.writer(f); w.writerow(["run_id", "experiment_id", "seed"]); w.writerows(runs)
with open(f"{UC}/overrides_def.csv", "w", newline="") as f: w = csv.writer(f); w.writerow(["scope", "xpath", "value"]); w.writerows(ov)
print(len(runs), "runs,", len(ov), "override rows")
