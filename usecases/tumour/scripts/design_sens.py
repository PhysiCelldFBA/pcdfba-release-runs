"""Write runs_sens.csv, overrides_sens.csv and results/sens_design.csv: the tumour sensitivity block (TUM-GRID, TUM-SOB, TUM-NTH, TUM-LONG) at the second reference (dt 0.005 min, TUM5-*), 24 h.
All runs: seed 50, dt_diffusion = intracellular_dt 0.005 with the necrosis increase (8.35e-8) and diffusion_steps (600) scaled as in TUM-DT, 56 threads, ledger on, no SVG, snapshots every 240 min (24 h = 6 steps).
  TUM-GRID  vessel O2 x vessel glucose, 5 x 5 (the centre, 0.5 / 0.5 mM, is TUM5-REF and is not rerun: 24 runs). O2 0.05 0.1 0.25 0.5 1.0 mM; glucose 0.1 0.25 0.5 1.0 5.0 mM.
  TUM-SOB   16-point scrambled Sobol (scipy, seed 1) on 6 factors, each multiplying the reference value log-uniformly in [0.5, 2]: O2 Vmax, O2 Km, glucose Vmax, glucose Km, D(O2), D(glucose); at vessel O2 0.5 (reference) and 0.1 mM (32 runs).
  TUM-NTH   death_flux_threshold of R_ATPM (reference 1.07) at 0.80, 0.96, 1.18, 1.34 (+/-12% and +/-25%; 4 runs).
  TUM-LONG  72 h runs at four corners of the grid (O2, glucose) = (0.05, 0.1), (0.05, 5), (1.0, 0.1), (1.0, 5) to check that 24 h is long enough (4 runs).
Total 64 runs. usage: design_sens.py"""
import csv, os
import numpy as np
from scipy.stats import qmc
UC = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); P = "/PhysiCell_settings"
CD = P + "/cell_definitions/cell_definition/phenotype/intracellular"
VAR = P + "/microenvironment_setup/variable[@name='%s']"
runs, ov = [], []; design = []
run = lambda rid, eid, seed=50: runs.append((rid, eid, seed)); put = lambda s, x, v: ov.append((s, x, v))
INC0, DT0, STEPS0, DT = 1.67e-7, 0.01, 300, 0.005
# all sensitivity runs: files from config/shipped, ledger on, no SVG, 56 threads (as design_def.py), 24 h, snapshots every 240 min, dt 0.005
put("TUM-GRID", P + "/overall/max_time", 1440); put("TUM-SOB", P + "/overall/max_time", 1440); put("TUM-NTH", P + "/overall/max_time", 1440)
put("TUM-LONG", P + "/overall/max_time", 4320)
def base(rid):
    put(rid, P + "/overall/dt_diffusion", DT); put(rid, CD + "/settings/intracellular_dt", DT)
    put(rid, CD + "/death_model/death_rate_increase", float("%.4g" % (INC0 * DT / DT0))); put(rid, P + "/user_parameters/diffusion_steps", int(round(STEPS0 * DT0 / DT)))
    put(rid, P + "/save/full_data/interval", 240)
def vessel(rid, sub, c):
    put(rid, VAR % sub + "/Dirichlet_boundary_condition", c); put(rid, VAR % sub + "/Dirichlet_options/boundary_value[@ID='xmin']", c)
tm = lambda sub, k: CD + f"/transport_model/exchange[@substrate='{sub}']/{k}"
O2S, GLS = (0.05, 0.1, 0.25, 0.5, 1.0), (0.1, 0.25, 0.5, 1.0, 5.0)
for o in O2S:
    for g in GLS:
        if (o, g) == (0.5, 0.5): continue
        rid = f"TUM-GRID_o{o}_g{g}_s50"; run(rid, "TUM-GRID"); base(rid); vessel(rid, "O2", o); vessel(rid, "D-Glucose", g)
REF = dict(vmax_o2=1.0, km_o2=0.5, vmax_glc=1.0, km_glc=5.0, d_o2=100000, d_glc=8921)
u = qmc.Sobol(d=6, scramble=True, seed=1).random(16); f = 2.0 ** (2 * u - 1)
for o in (0.5, 0.1):
    for i in range(16):
        rid = f"TUM-SOB_o{o}_p{i + 1:02d}_s50"; run(rid, "TUM-SOB"); base(rid)
        if o != 0.5: vessel(rid, "O2", o)
        v = {k: float("%.5g" % (REF[k] * f[i][j])) for j, k in enumerate(REF)}
        put(rid, tm("O2", "Vmax"), v["vmax_o2"]); put(rid, tm("O2", "Km"), v["km_o2"]); put(rid, tm("D-Glucose", "Vmax"), v["vmax_glc"]); put(rid, tm("D-Glucose", "Km"), v["km_glc"])
        put(rid, VAR % "O2" + "/physical_parameter_set/diffusion_coefficient", v["d_o2"]); put(rid, VAR % "D-Glucose" + "/physical_parameter_set/diffusion_coefficient", v["d_glc"])
        design.append(dict(run=rid, o2_vessel=o, point=i + 1, **{k: v[k] for k in REF}, **{"f_" + k: float("%.4f" % f[i][j]) for j, k in enumerate(REF)}))
for th in (0.80, 0.96, 1.18, 1.34):
    rid = f"TUM-NTH_th{th}_s50"; run(rid, "TUM-NTH"); base(rid); put(rid, CD + "/death_model/death_flux_threshold", th)
for o, g in ((0.05, 0.1), (0.05, 5.0), (1.0, 0.1), (1.0, 5.0)):
    rid = f"TUM-LONG_o{o}_g{g}_s50"; run(rid, "TUM-LONG"); base(rid); vessel(rid, "O2", o); vessel(rid, "D-Glucose", g)
with open(f"{UC}/runs_sens.csv", "w", newline="") as fh: w = csv.writer(fh, lineterminator="\n"); w.writerow(["run_id", "experiment_id", "seed"]); w.writerows(runs)
with open(f"{UC}/overrides_sens.csv", "w", newline="") as fh: w = csv.writer(fh, lineterminator="\n"); w.writerow(["scope", "xpath", "value"]); w.writerows(ov)
os.makedirs(f"{UC}/results", exist_ok=True)
with open(f"{UC}/results/sens_design.csv", "w", newline="") as fh: w = csv.DictWriter(fh, fieldnames=list(design[0]), lineterminator="\n"); w.writeheader(); w.writerows(design)
print(len(runs), "runs,", len(ov), "override rows")
