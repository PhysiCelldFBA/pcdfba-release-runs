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
# Block 2 (2026-11-07), per variant at seed 50: convergence in dt and dx, pFBA, and the sensitivity block. Criteria fixed before the runs (experiments.csv).
INTRA = P + "/cell_definitions/cell_definition/phenotype/intracellular"; GLU = P + "/microenvironment_setup/variable[@name='glucose']"
KMX = INTRA + "/transport_model/exchange[@substrate='glucose']/Km"; XMAX = GLU + "/Dirichlet_options/boundary_value[@ID='xmax']"; DGLU = GLU + "/physical_parameter_set/diffusion_coefficient"; COST = P + "/user_parameters/motility_cost_at_vmax"
VAR = {"REF": dict(km=0.5, cost=48.9), "GR": dict(km=0.02, cost=66.4)}; DG, BND = 46800.0, 0.65
def setdt(rid, dt): put(rid, P + "/overall/dt_diffusion", dt); put(rid, INTRA + "/settings/intracellular_dt", dt)
def setdx(rid, dx): put(rid, P + "/domain/dx", dx); put(rid, P + "/domain/dy", dx)
for v, d in VAR.items():
    for dt in (0.002, 0.0005): rid = f"MOT-DT-{v}_dt{dt}_s50"; run(rid, f"MOT-DT-{v}"); setdt(rid, dt)
    for dx in (2.5, 10): rid = f"MOT-DX-{v}_dx{dx}_s50"; run(rid, f"MOT-DX-{v}"); setdx(rid, dx)
    rid = f"MOT-PFBA-{v}_s50"; run(rid, f"MOT-PFBA-{v}"); put(rid, "+" + INTRA + "/growth_model", '<secondary_objective type="parsimonious" enabled="true" fraction_of_optimum="1.0"/>')
    for fk in (0.5, 1, 2):
        for fb in (0.5, 1, 2):
            if fk == 1 and fb == 1: continue
            rid = f"MOT-SENS-{v}_Km{fk:g}_B{fb:g}_s50"; run(rid, f"MOT-SENS-{v}"); put(rid, KMX, round(d["km"] * fk, 6)); put(rid, XMAX, round(BND * fb, 6))
    for f in (0.5, 2): rid = f"MOT-SENS-{v}_cost{f:g}_s50"; run(rid, f"MOT-SENS-{v}"); put(rid, COST, round(d["cost"] * f, 4))
    for f in (0.5, 2): rid = f"MOT-SENS-{v}_D{f:g}_s50"; run(rid, f"MOT-SENS-{v}"); put(rid, DGLU, DG * f)
# Block 3 (2026-11-07): diagnosis of the MOT-REF dependence on dt and dx. At dt 0.0005 the founder (x = -140, on a voxel face) never moves: the glucose of its voxel is 0.021 mM, the ATP-optimised
# solution cannot meet the basal ATP flux (8.39) at that uptake, and the motility flux is 0. The glucose near the x_min wall (Dirichlet 0) is above the steady profile by about 19 x dt mM. Test: other start positions
# (voxel centres -142.5 and -137.5, and -130, -120, where the steady glucose is 0.016, 0.027, 0.039, 0.065 mM... see RESULTS.md) at three dt, two dx at -120, and intermediate dt at the shipped position.
START = P + "/initial_conditions/cell_positions/folder"
for x0 in ("-142.5", "-137.5", "-130", "-120"):
    for dt in (0.0005, 0.001, 0.002): rid = f"MOT-START-REF_x{x0}_dt{dt}_s50"; run(rid, "MOT-START-REF"); put(rid, START, f"./config/start/x{x0}"); setdt(rid, dt)
for dx in (2.5, 10): rid = f"MOT-START-REF_x-120_dx{dx}_s50"; run(rid, "MOT-START-REF"); put(rid, START, "./config/start/x-120"); setdx(rid, dx)
for dt in (0.00075, 0.0014): rid = f"MOT-DT-REF_dt{dt}_s50"; run(rid, "MOT-DT-REF"); setdt(rid, dt)
# Block 4 (2026-11-07), MOT2-*: the full set again with the founder cell starting at x = -120 um (config/start/x-120/cells.csv) instead of -140 (on a voxel face, 0.01-0.02 mM of glucose in the steady
# profile, below the 0.032 mM that the basal ATP flux needs): reference and seeds, convergence in dt and dx, pFBA and the sensitivity block, same overrides as MOT-* otherwise. The MOT-* runs (start -140) are kept as the record of the shipped start.
OLD = [r for r in runs if r[0].startswith(("MOT-REF", "MOT-GR", "MOT-DT-", "MOT-DX-", "MOT-PFBA-", "MOT-SENS-")) and "dt0.00075" not in r[0] and "dt0.0014" not in r[0]]
new = lambda t: t.replace("MOT-", "MOT2-", 1)
for rid, eid, seed in OLD:
    run(new(rid), new(eid), seed); put(new(rid), START, "./config/start/x-120")
    for sc, xp, v in [o for o in ov if o[0] == rid]: put(new(rid), xp, v)
# Block 5 (2026-11-07): one run of the motility-only variant with Km(glucose) 0.8 mM (the value reported in the manuscript, S5.6), founder at x = -120 um, to see whether the regime without division within 10 h is reproduced.
rid = "MOT2-KM-REF_Km0.8_s50"; run(rid, "MOT2-KM-REF"); put(rid, START, "./config/start/x-120"); put(rid, KMX, 0.8)
# Block 6 (2026-11-07), MOT3-*: regimes selected by the glucose supply. Both behaviours use the growth-recovery kinetics (Km 0.02 mM, threshold 0.35, cost 66.4 = the file motility_proliferation); the founder starts at x = -120 um.
#   high supply = boundary glucose 0.65 mM (the MOT2-GR reference, not rerun); low supply = 0.325 mM, below the switching threshold (0.35 mM): the cell is expected to reach the boundary and stay in motility mode.
# MOT3-LOW: seeds 50-52 and repeat, dt 0.0005 / 0.002, dx 2.5 / 10, pFBA, sensitivity (Km x boundary grid, cost, D) around the low-supply reference; MOT3-KM-LOW / MOT3-KM-HIGH: Km scan (speed against Km) at the two supplies.
LOWB = 0.325
def low(rid, eid, seed=50, b=LOWB): run(rid, eid, seed); put(rid, START, "./config/start/x-120"); put(rid, XMAX, b)
for sd in (50, 51, 52): low(f"MOT3-LOW_s{sd}", "MOT3-LOW", sd)
low("MOT3-LOW-REP_s50", "MOT3-LOW")
for dt in (0.0005, 0.002): rid = f"MOT3-DT-LOW_dt{dt}_s50"; low(rid, "MOT3-DT-LOW"); setdt(rid, dt)
for dx in (2.5, 10): rid = f"MOT3-DX-LOW_dx{dx}_s50"; low(rid, "MOT3-DX-LOW"); setdx(rid, dx)
rid = "MOT3-PFBA-LOW_s50"; low(rid, "MOT3-PFBA-LOW"); put(rid, "+" + INTRA + "/growth_model", '<secondary_objective type="parsimonious" enabled="true" fraction_of_optimum="1.0"/>')
for fk in (0.5, 1, 2):
    for fb in (0.5, 1, 2):
        if fk == 1 and fb == 1: continue
        rid = f"MOT3-SENS-LOW_Km{fk:g}_B{fb:g}_s50"; low(rid, "MOT3-SENS-LOW", b=round(LOWB * fb, 6)); put(rid, KMX, round(0.02 * fk, 6))
for f in (0.5, 2): rid = f"MOT3-SENS-LOW_cost{f:g}_s50"; low(rid, "MOT3-SENS-LOW"); put(rid, COST, round(66.4 * f, 4))
for f in (0.5, 2): rid = f"MOT3-SENS-LOW_D{f:g}_s50"; low(rid, "MOT3-SENS-LOW"); put(rid, DGLU, DG * f)
for km in (0.005, 0.08, 0.16, 0.32, 0.64, 1.28):
    rid = f"MOT3-KM-LOW_Km{km:g}_s50"; low(rid, "MOT3-KM-LOW"); put(rid, KMX, km)
    rid = f"MOT3-KM-HIGH_Km{km:g}_s50"; low(rid, "MOT3-KM-HIGH", b=0.65); put(rid, KMX, km)
# Block 7 (2026-11-08), MOT4-HIGH: the high supply moved from 0.65 to 1.0 mM. At 0.65 mM the glucose at the wall region is 0.608 mM, 1.4% above the switching threshold (0.6 mM): seed 50 stayed in ATP mode at the wall (recorded mode 1, growth_rate 0) while seeds 51 and 52 switched.
# At 1.0 mM the wall region is far above the threshold. Same growth-recovery kinetics (Km 0.02 mM, threshold 0.35, cost 66.4), founder at x = -120 um: seeds 50-52 and repeat, dt 0.0005 / 0.002, dx 2.5 / 10, pFBA, and the Km scan (MOT4-KM-HIGH). The low supply (MOT3-LOW) is unchanged.
HIGHB = 1.0
def high(rid, eid, seed=50, b=HIGHB): run(rid, eid, seed); put(rid, START, "./config/start/x-120"); put(rid, XMAX, b)
for sd in (50, 51, 52): high(f"MOT4-HIGH_s{sd}", "MOT4-HIGH", sd)
high("MOT4-HIGH-REP_s50", "MOT4-HIGH")
for dt in (0.0005, 0.002): rid = f"MOT4-DT-HIGH_dt{dt}_s50"; high(rid, "MOT4-DT-HIGH"); setdt(rid, dt)
for dx in (2.5, 10): rid = f"MOT4-DX-HIGH_dx{dx}_s50"; high(rid, "MOT4-DX-HIGH"); setdx(rid, dx)
rid = "MOT4-PFBA-HIGH_s50"; high(rid, "MOT4-PFBA-HIGH"); put(rid, "+" + INTRA + "/growth_model", '<secondary_objective type="parsimonious" enabled="true" fraction_of_optimum="1.0"/>')
for km in (0.005, 0.01, 0.04, 0.08, 0.16, 0.32, 0.64, 1.28, 2.56):
    rid = f"MOT4-KM-HIGH_Km{km:g}_s50"; high(rid, "MOT4-KM-HIGH"); put(rid, KMX, km)
with open(f"{UC}/runs_def.csv", "w", newline="") as f: w = csv.writer(f, lineterminator="\n"); w.writerow(["run_id", "experiment_id", "seed"]); w.writerows(runs)
with open(f"{UC}/overrides_def.csv", "w", newline="") as f: w = csv.writer(f, lineterminator="\n"); w.writerow(["scope", "xpath", "value"]); w.writerows(ov)
print(len(runs), "runs,", len(ov), "override rows")
