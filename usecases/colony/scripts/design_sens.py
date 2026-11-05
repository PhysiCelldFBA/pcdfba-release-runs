"""Write runs_sens.csv and overrides_sens.csv: the reduced colony sensitivity set (COL-BND, COL-KIN, COL-PFBA), 17 runs, designed 2026-11-05 before any of them was run.
Base = the pre-grown colony of COL-DT_dt0.0005_s50 (R0 200 um disc of 87,685 cells, +-500 um box, dx 2, 2D, 4 mM glucose / 5 mM O2 on the four faces, necrosis x10 = 1.67e-4 per update, dt_diffusion = intracellular_dt 0.0005,
dt_mechanics = dt_phenotype 0.001, ledger on, 112 threads, seed 50, plain pinned binary) run for 6 min with snapshots every 2 min; the reference is that run at 6 min (not rerun).
  COL-BND  boundary grid, glucose {2, 4, 16} mM x O2 {2.5, 5, 10} mM on the four faces and as initial condition (9 points; the centre is the reference: 8 runs). 10 mM O2 with 16 mM glucose is close to the shipped 10 / 11 mM.
  COL-KIN  8-point scrambled Sobol (scipy, seed 1) on 6 factors at the reference boundary: Vmax of glucose, of oxygen and of acetate, Km of glucose, D of glucose (each x0.5-2, log-uniform) and the necrosis threshold of R_ATPM
           (8.39 x 0.75-1.25, uniform). Km of oxygen (0.5 uM) and D of oxygen are left out: both are far from limiting here (cap-free runs, oxygen diffusion length >> colony). 8 points are a screen, not variance indices.
  COL-PFBA parsimonious secondary objective at the reference settings (1 run).
usage: design_sens.py   (also writes results/sens_design.csv)"""
import csv, os
import numpy as np
from lxml import etree
from scipy.stats import qmc
UC = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); P = "/PhysiCell_settings"; CD = P + "/cell_definitions/cell_definition/phenotype/intracellular"
base = etree.parse(f"{UC}/config/shipped/PhysiCell_settings.xml").getroot(); g = lambda xp: float(base.xpath(P + xp)[0].text)
runs, ov, design = [], [], []
run = lambda rid, eid, seed=50: runs.append((rid, eid, seed)); put = lambda s, x, v: ov.append((s, x, v))
SBML = CD + "/settings/sbml_filename"; LEDGER = ("+" + CD + "/settings", "<dfba_conservation_ledger>true</dfba_conservation_ledger>")
EXPS = ["COL-BND", "COL-KIN", "COL-PFBA"]
for e in EXPS:   # = COL-DT_dt0.0005_s50 (design_def.py), 6 min instead of 10
    put(e, SBML, "./config/shipped/Ecoli_core.xml"); put(e, *LEDGER)
    put(e, P + "/overall/max_time", 6); put(e, P + "/save/full_data/interval", 2); put(e, P + "/save/SVG/enable", "false")
    put(e, P + "/user_parameters/colony_radius", 200); put(e, CD + "/death_model/death_rate_increase", 1.67e-4)
    put(e, P + "/overall/dt_diffusion", 0.0005); put(e, CD + "/settings/intracellular_dt", 0.0005)
    put(e, P + "/overall/dt_mechanics", 0.001); put(e, P + "/overall/dt_phenotype", 0.001)
    put(e, P + "/cell_definitions/cell_definition/phenotype/mechanics/cell_cell_repulsion_strength", 1000.0)
VAR = P + "/microenvironment_setup/variable[@name='%s']"
def boundary(rid, sub, c):   # initial condition and the four faces
    put(rid, VAR % sub + "/initial_condition", c)
    for face in ("xmin", "xmax", "ymin", "ymax"): put(rid, VAR % sub + f"/Dirichlet_options/boundary_value[@ID='{face}']", c)
for gl in (2, 4, 16):
    for o2 in (2.5, 5, 10):
        if (gl, o2) == (4, 5): continue
        rid = f"COL-BND_g{gl}_o{o2}_s50"; run(rid, "COL-BND"); boundary(rid, "glucose", gl); boundary(rid, "oxygen", o2)
tm = lambda sub, k: CD + f"/transport_model/exchange[@substrate='{sub}']/{k}"
REF = dict(vmax_glc=float(base.xpath(tm("glucose", "Vmax"))[0].text), km_glc=float(base.xpath(tm("glucose", "Km"))[0].text), vmax_o2=float(base.xpath(tm("oxygen", "Vmax"))[0].text),
           vmax_ac=float(base.xpath(tm("acetate", "Vmax"))[0].text), d_glc=g("/microenvironment_setup/variable[@name='glucose']/physical_parameter_set/diffusion_coefficient"),
           thr=float(base.xpath(CD + "/death_model/death_flux_threshold")[0].text))
u = qmc.Sobol(d=6, scramble=True, seed=1).random(8); f = 2.0 ** (2 * u[:, :5] - 1); fth = 0.75 + 0.5 * u[:, 5]
keys = ["vmax_glc", "km_glc", "vmax_o2", "vmax_ac", "d_glc"]
for i in range(8):
    rid = f"COL-KIN_p{i + 1:02d}_s50"; run(rid, "COL-KIN"); v = {k: float("%.5g" % (REF[k] * f[i][j])) for j, k in enumerate(keys)}; thr = float("%.4g" % (REF["thr"] * fth[i]))
    put(rid, tm("glucose", "Vmax"), v["vmax_glc"]); put(rid, tm("glucose", "Km"), v["km_glc"]); put(rid, tm("oxygen", "Vmax"), v["vmax_o2"]); put(rid, tm("acetate", "Vmax"), v["vmax_ac"])
    put(rid, VAR % "glucose" + "/physical_parameter_set/diffusion_coefficient", v["d_glc"]); put(rid, CD + "/death_model/death_flux_threshold", thr)
    design.append(dict(run=rid, **v, necrosis_threshold=thr, **{"f_" + k: float("%.4f" % f[i][j]) for j, k in enumerate(keys)}, f_threshold=float("%.4f" % fth[i])))
run("COL-PFBA_s50", "COL-PFBA"); put("COL-PFBA", "+" + CD + "/growth_model", '<secondary_objective type="parsimonious" enabled="true" fraction_of_optimum="1.0"/>')
with open(f"{UC}/runs_sens.csv", "w", newline="") as fh: w = csv.writer(fh, lineterminator="\n"); w.writerow(["run_id", "experiment_id", "seed"]); w.writerows(runs)
with open(f"{UC}/overrides_sens.csv", "w", newline="") as fh: w = csv.writer(fh, lineterminator="\n"); w.writerow(["scope", "xpath", "value"]); w.writerows(ov)
os.makedirs(f"{UC}/results", exist_ok=True)
with open(f"{UC}/results/sens_design.csv", "w", newline="") as fh: w = csv.DictWriter(fh, fieldnames=list(design[0]), lineterminator="\n"); w.writeheader(); w.writerows(design)
print(len(runs), "runs,", len(ov), "override rows")
