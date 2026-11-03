"""Write runs_def.csv and overrides_def.csv: the unit-test experiments (see experiments.csv). Values are absolute.
usage: design_def.py"""
import csv, os, random
UC = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); P = "/PhysiCell_settings"
runs, ov = [], []
run = lambda rid, eid, seed=50: runs.append((rid, eid, seed))
put = lambda scope, xp, v: ov.append((scope, xp, v))
SB = P + "/cell_definitions/cell_definition/phenotype/intracellular/settings"
SBML1 = P + "/cell_definitions/cell_definition/phenotype/intracellular/settings/sbml_filename"
LEDGER = ("+" + SB, "<dfba_conservation_ledger>true</dfba_conservation_ledger>")
GLC = P + "/microenvironment_setup/variable[@name='glucose']/"
DX = [0.5, 1, 2, 4]; DT = [0.001, 0.005, 0.01, 0.02]
def grid(eid, T):
    put(eid, SBML1, "config/shipped/Ecoli_core.xml"); put(eid, *LEDGER)
    put(eid, P + "/overall/max_time", T); put(eid, P + "/parallel/omp_num_threads", 1)
    put(eid, P + "/save/full_data/interval", 30); put(eid, P + "/save/SVG/interval", 30); put(eid, P + "/save/SVG/enable", "false")
def point(rid, dx, dt):
    for a in ("dx", "dy", "dz"): put(rid, P + f"/domain/{a}", dx)
    for a in ("dt_diffusion", "dt_mechanics", "dt_phenotype"): put(rid, P + f"/overall/{a}", dt)   # mechanics/phenotype never below diffusion
    put(rid, SB + "/intracellular_dt", dt)
# UT-BASE: default test (240 min), same config 5 times
put("UT-BASE", SBML1, "config/shipped/Ecoli_core.xml"); put("UT-BASE", *LEDGER); put("UT-BASE", P + "/parallel/omp_num_threads", 1)
for k, rid in enumerate(["UT-BASE_s50"] + [f"UT-BASE_s50_rep{i}" for i in range(1, 5)]): run(rid, "UT-BASE")
# UT-ATPM / UT-IML1515: the default test (UT-BASE settings) with another model: the core model with ATPM at 8.39 (models/, made by scripts/make_models.py) and iML1515
OBJ = "/PhysiCell_settings/cell_definitions/cell_definition/phenotype/intracellular/growth_model/objective_reaction"
for e, f, o in (("UT-ATPM", "models/Ecoli_core_ATPM8.39.xml", None), ("UT-IML1515", "models/iML1515.xml", "R_BIOMASS_Ec_iML1515_core_75p37M")):
    put(e, SBML1, f); put(e, *LEDGER); put(e, P + "/parallel/omp_num_threads", 1)
    if o: put(e, OBJ, o)
    run(f"{e}_s50", e)
# UT-DXDT-A / B: one cell, closed 16 um box, 120 min
grid("UT-DXDT-A", 120); grid("UT-DXDT-B", 120)
put("UT-DXDT-B", GLC + "physical_parameter_set/diffusion_coefficient", 0); put("UT-DXDT-B", GLC + "initial_condition", 1)
for e, tag in (("UT-DXDT-A", "A"), ("UT-DXDT-B", "B")):
    for dx in DX:
        for dt in DT: rid = f"{e}_dx{dx}_dt{dt}"; run(rid, e); point(rid, dx, dt)
# UT-CONV-BOX: 8 cells (created at the origin, pushed apart by mechanics) in the same box, growth and division, 120 min
grid("UT-CONV-BOX", 120); put("UT-CONV-BOX", P + "/user_parameters/number_of_cells", 8)
put("UT-CONV-BOX", GLC + "initial_condition", 50)   # 10 mM (default) is used up by 90 min
for dx in DX:
    for dt in DT: rid = f"UT-CONV-BOX_dx{dx}_dt{dt}"; run(rid, "UT-CONV-BOX"); point(rid, dx, dt)
# UT-THREADS: the CONV-BOX reference point (dx 1, dt 0.01) at 1, 4, 16 threads
grid("UT-THREADS", 120); put("UT-THREADS", P + "/user_parameters/number_of_cells", 8)
put("UT-THREADS", GLC + "initial_condition", 50)   # same box as UT-CONV-BOX
for th in (1, 4, 16):
    rid = f"UT-THREADS_t{th}"; run(rid, "UT-THREADS"); point(rid, 1, 0.01); put(rid, P + "/parallel/omp_num_threads", th)
# UT-VOX-*: the shipped voxel-sharing configs unchanged (0.1 min), only the SBML path
for c in ("n1", "n4", "n8"):
    e = f"UT-VOX-{c.upper()}"; put(e, SBML1, "config/shipped/Ecoli_core.xml"); run(f"{e}_s50", e)
e = "UT-VOX-MIXED"; run(f"{e}_s50", e)
for i in (1, 2): put(e, f"(/PhysiCell_settings/cell_definitions/cell_definition/phenotype/intracellular/settings/sbml_filename)[{i}]", "config/shipped/Ecoli_core.xml")
# UT-ORDER: 12 cells (types A B C N, three each) in one voxel created in 5 orders, at 1 and 4 threads
rows = [t for t in "ABCN" for _ in range(3)]; rng = random.Random(1)
orders = {"identity": rows[:], "reversed": rows[::-1]}
for k in range(3): p = rows[:]; rng.shuffle(p); orders[f"rand{k}"] = p
for i in (1, 2, 3, 4): put("UT-ORDER", f"(/PhysiCell_settings/cell_definitions/cell_definition/phenotype/intracellular/settings/sbml_filename)[{i}]", "config/shipped/Ecoli_core.xml")
for name, p in orders.items():
    for th in (1, 4):
        rid = f"UT-ORDER_{name}_t{th}"; run(rid, "UT-ORDER"); put(rid, "@cells.order", "".join(p)); put(rid, P + "/parallel/omp_num_threads", th)
with open(f"{UC}/runs_def.csv", "w", newline="") as f: w = csv.writer(f); w.writerow(["run_id", "experiment_id", "seed"]); w.writerows(runs)
with open(f"{UC}/overrides_def.csv", "w", newline="") as f: w = csv.writer(f); w.writerow(["scope", "xpath", "value"]); w.writerows(ov)
print(len(runs), "runs,", len(ov), "override rows")
