"""Write runs_def.csv and overrides_def.csv: the definitive acetate experiments (ACE-REF ... ACE-PULSE) on top of overrides.csv (the `*` common settings).
All values are absolute so the tables can be read without the base; sensitivity values are the base value x0.5 and x2.
usage: design_def.py   (run from anywhere; rewrites the two files)"""
import csv, os
from lxml import etree
UC = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
base = etree.parse(f"{UC}/config/shipped/PhysiCell_settings.xml").getroot(); P = "/PhysiCell_settings"
g = lambda xp: float(base.xpath(P + xp)[0].text)
fmt = lambda v: f"{v:.6g}"
runs, ov = [], []
def run(rid, eid, seed=50): runs.append((rid, eid, seed))
def put(scope, xp, val): ov.append((scope, xp, val))
LEDGER = ("+" + P + "/cell_definitions/cell_definition/phenotype/intracellular/settings", "<dfba_conservation_ledger>true</dfba_conservation_ledger>")
DEF = ["ACE-REF", "ACE-SEED", "ACE-REP", "ACE-CONV", "ACE-DTM", "ACE-PFBA", "ACE-SENS", "ACE-PULSE", "ACE-FVA1", "ACE-FVA6S"]
for e in DEF: put(e, *LEDGER)
run("ACE-REF_s50", "ACE-REF"); run("ACE-SEED_s51", "ACE-SEED", 51); run("ACE-SEED_s52", "ACE-SEED", 52); run("ACE-REP_s50", "ACE-REP")
IC = P + "/cell_definitions/cell_definition/phenotype/intracellular/settings/intracellular_dt"
for dx in (2, 6):
    r = f"ACE-CONV_dx{dx}_s50"; run(r, "ACE-CONV"); [put(r, P + f"/domain/{a}", dx) for a in ("dx", "dy", "dz")]
for dt in (0.001, 0.05):
    r = f"ACE-CONV_dt{dt}_s50"; run(r, "ACE-CONV"); put(r, P + "/overall/dt_diffusion", dt); put(r, IC, dt)
run("ACE-DTM_m0.001_s50", "ACE-DTM"); put("ACE-DTM", P + "/overall/dt_mechanics", 0.001); put("ACE-DTM", P + "/overall/dt_phenotype", 0.001)
run("ACE-PFBA_s50", "ACE-PFBA")
put("ACE-PFBA", "+" + P + "/cell_definitions/cell_definition/phenotype/intracellular/growth_model", '<secondary_objective type="parsimonious" enabled="true" fraction_of_optimum="1.0"/>')
T = P + "/cell_definitions/cell_definition/phenotype/intracellular/transport_model/exchange[@substrate='%s']/%s"
V = P + "/microenvironment_setup/variable[@name='%s']/%s"
sens = [(f"{k}_{s}", T % (s, k)) for s in ("glucose", "acetate", "oxygen") for k in ("Km", "Vmax")]
sens += [(f"D_{s}", V % (s, "physical_parameter_set/diffusion_coefficient")) for s in ("glucose", "acetate", "oxygen")]
sens += [("glc0", V % ("glucose", "initial_condition")), ("pulse", P + "/user_parameters/o2_concentration"), ("tpulse", P + "/user_parameters/reinject_oxygen_time")]
for name, xp in sens:
    b = float(base.xpath(xp)[0].text)
    for f in (0.5, 2):
        r = f"ACE-SENS_{name}_x{f}_s50"; run(r, "ACE-SENS"); put(r, xp, fmt(b * f))
for f in (0.5, 2):
    r = f"ACE-SENS_ncells_x{f}_s50"; run(r, "ACE-SENS"); put(r, "@cells.n", int(20 * f))
# same cell density in a box of different volume (z extent changed, x and y kept at 36 um; dx 4 stays): 40 cells in 36 x 36 x 72 (x2) and 9 cells in 36 x 36 x 16 (x0.44; 18 um is not a multiple of dx 4), density 4.29e-4 and 4.34e-4 per um3 against 4.29e-4 in the reference (20 cells in 36^3)
for f, n, hz in ((2, 40, 36), (0.44, 9, 8)):
    r = f"ACE-SENS_vol_x{f}_s50"; run(r, "ACE-SENS"); put(r, "@cells.n", n); put(r, P + "/domain/z_min", -hz); put(r, P + "/domain/z_max", hz)
put("ACE-PULSE", P + "/overall/max_time", 900)
for s in (50, 51, 52): r = f"ACE-PULSE_o2_10_s{s}"; run(r, "ACE-PULSE", s); put(r, P + "/user_parameters/o2_concentration", 10)
for o in (20, 40, 80): r = f"ACE-PULSE_o2_{o}_s50"; run(r, "ACE-PULSE"); put(r, P + "/user_parameters/o2_concentration", o)
run("ACE-FVA1_s50", "ACE-FVA1"); put("ACE-FVA1", P + "/overall/max_time", 150); put("ACE-FVA1", P + "/save/full_data/interval", 1)
run("ACE-FVA6S_s50", "ACE-FVA6S"); put("ACE-FVA6S", P + "/overall/max_time", 80); put("ACE-FVA6S", P + "/save/full_data/interval", 0.1)
with open(f"{UC}/runs_def.csv", "w", newline="") as f: w = csv.writer(f); w.writerow(["run_id", "experiment_id", "seed"]); w.writerows(runs)
with open(f"{UC}/overrides_def.csv", "w", newline="") as f: w = csv.writer(f); w.writerow(["scope", "xpath", "value"]); w.writerows(ov)
print(len(runs), "runs,", len(ov), "override rows")
