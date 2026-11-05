"""Write runs_dt125.csv / overrides_dt125.csv: TUM-DT125, one run at dt_diffusion = intracellular_dt 0.00125 min (half of the smallest dt of TUM-DT), 24 h, seed 50, to test whether the O2 field at 6-24 h still moves with dt (compare_dt_profiles.py). Scaling as in design_def.py (necrosis increase 1.67e-7 x dt/0.01, diffusion_steps 300 x 0.01/dt); snapshots every 60 min as the other TUM-DT runs. usage: design_dt125.py"""
import csv, os
UC = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); P = "/PhysiCell_settings"; CD = P + "/cell_definitions/cell_definition/phenotype/intracellular"; DT = 0.00125; rid = "TUM-DT125_dt0.00125_s50"
ov = [("TUM-DT125", P + "/overall/max_time", 1440), ("TUM-DT125", P + "/overall/dt_diffusion", DT), ("TUM-DT125", CD + "/settings/intracellular_dt", DT),
      ("TUM-DT125", CD + "/death_model/death_rate_increase", float("%.4g" % (1.67e-7 * DT / 0.01))), ("TUM-DT125", P + "/user_parameters/diffusion_steps", int(round(300 * 0.01 / DT)))]
with open(f"{UC}/runs_dt125.csv", "w", newline="") as f: w = csv.writer(f, lineterminator="\n"); w.writerow(["run_id", "experiment_id", "seed"]); w.writerow([rid, "TUM-DT125", 50])
with open(f"{UC}/overrides_dt125.csv", "w", newline="") as f: w = csv.writer(f, lineterminator="\n"); w.writerow(["scope", "xpath", "value"]); w.writerows(ov)
print(rid, len(ov), "override rows")
