"""Write runs_def.csv and overrides_def.csv: the definitive tumour experiments (TUM-*) on top of config/shipped/PhysiCell_settings.xml.
Shipped base: 360 x 120 x 30 um box (x from 0, vessel at x = 0), dx = dy = dz = 15 um (two layers), 240 cells (cells_column.csv), 72 h, dt_diffusion = intracellular_dt 0.01 min, dt_mechanics 0.1, dt_phenotype 6, 56 threads, snapshots every 60 min.
Necrosis: the death rate of a cell whose LP is infeasible grows by death_rate_increase at each dFBA update, so the ramp is increase/dt per min2 and depends on dt. In the dt series the increase is scaled with dt (increase = 1.67e-7 x dt / 0.01, as in COL-NECS)
so that the ramp is the shipped 1.67e-5 /min2 in every run. The pre-simulation diffusion (user parameter diffusion_steps, 300 steps of dt_diffusion = 3 min) is scaled the same way (steps = 300 x 0.01 / dt) so that every run starts from the same physical time of diffusion.
Voxel: the cell volume is 2494 um3 and the voxel dx x dx x 15 um3 must not be smaller: dx 15 (3375 um3) is the finest; dx 20 (6000) and 30 (13500) are run; finer voxels are not.
usage: design_def.py"""
import csv, os
UC = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); P = "/PhysiCell_settings"
runs, ov = [], []
run = lambda rid, eid, seed=50: runs.append((rid, eid, seed)); put = lambda s, x, v: ov.append((s, x, v))
CD = P + "/cell_definitions/cell_definition/phenotype/intracellular"
# all TUM experiments: files from config/shipped, ledger on, no SVG, 56 threads (2 runs per MN5 node)
put("*", CD + "/settings/sbml_filename", "config/shipped/MCF7_Zielinski_2017_core_trimmed.sbml")
put("*", P + "/initial_conditions/cell_positions/folder", "./config/shipped")
put("*", "+" + CD + "/settings", "<dfba_conservation_ledger>true</dfba_conservation_ledger>")
put("*", P + "/save/SVG/enable", "false"); put("*", P + "/parallel/omp_num_threads", 56)
INC0, DT0, STEPS0 = 1.67e-7, 0.01, 300
def setdt(rid, dt):
    put(rid, P + "/overall/dt_diffusion", dt); put(rid, CD + "/settings/intracellular_dt", dt)
    put(rid, CD + "/death_model/death_rate_increase", float("%.4g" % (INC0 * dt / DT0))); put(rid, P + "/user_parameters/diffusion_steps", int(round(STEPS0 * DT0 / dt)))
run("TUM-REF_s50", "TUM-REF"); run("TUM-SEED_s51", "TUM-SEED", 51); run("TUM-SEED_s52", "TUM-SEED", 52); run("TUM-REP_s50", "TUM-REP")
run("TUM-PFBA_s50", "TUM-PFBA")
put("TUM-PFBA", "+" + CD + "/growth_model", '<secondary_objective type="parsimonious" enabled="true" fraction_of_optimum="1.0"/>')
for dx in (20, 30):
    rid = f"TUM-DX_dx{dx}_s50"; run(rid, "TUM-DX"); put(rid, P + "/domain/dx", dx); put(rid, P + "/domain/dy", dx)
for dt in (0.005, 0.02):
    rid = f"TUM-DT_dt{dt}_s50"; run(rid, "TUM-DT"); setdt(rid, dt)
with open(f"{UC}/runs_def.csv", "w", newline="") as f: w = csv.writer(f, lineterminator="\n"); w.writerow(["run_id", "experiment_id", "seed"]); w.writerows(runs)
with open(f"{UC}/overrides_def.csv", "w", newline="") as f: w = csv.writer(f, lineterminator="\n"); w.writerow(["scope", "xpath", "value"]); w.writerows(ov)
print(len(runs), "runs,", len(ov), "override rows")
