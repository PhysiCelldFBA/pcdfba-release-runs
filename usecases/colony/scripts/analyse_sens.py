"""Reduced colony sensitivity set (COL-BND, COL-KIN, COL-PFBA; design_sens.py) against the reference COL-DT_dt0.0005_s50, all at t = 6 min.
usage: analyse_sens.py   (every runs/COL-BND_* COL-KIN_* COL-PFBA_*; reference runs/COL-DT_dt0.0005_s50; uses pg_analyse.analyse for the annuli of 10 um)
Writes results/sens_summary.csv (one row per run at 6 min: dead %, core radius, living-rim width = R_colony - R_core, ring edges, acetate-using fraction of the living cells, capped % from results/cap_check_<run>.csv if present, differences to the reference and the criterion fixed on 2026-11-05:
core radius within one 10 um annulus AND dead fraction within 10% relative) and results/sens_profiles.csv (annulus profiles at 6 min). Does not touch pg_summary.csv / pg_profiles.csv."""
import sys, os, glob, numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); UC = HERE + "/.."; sys.path.insert(0, HERE)
from pg_analyse import analyse
REF = "COL-DT_dt0.0005_s50"; T = 6.0
runs = [REF] + sorted(os.path.basename(p) for pat in ("COL-BND_*", "COL-KIN_*", "COL-PFBA_*") for p in glob.glob(f"{UC}/runs/{pat}") if os.path.isdir(p))
rows, profs = [], []
for r in runs:
    a, b = analyse(f"{UC}/runs/{r}"); a = a[abs(a.t - T) < 1e-6]; b = b[abs(b.t - T) < 1e-6]
    if a.empty: print("no 6 min snapshot:", r); continue
    row = a.iloc[0].to_dict(); row["rim_width"] = row["R_colony"] - row["R_core"]
    cc = f"{UC}/results/cap_check_{r}.csv"
    if os.path.exists(cc):
        c = pd.read_csv(cc); c = c[abs(c.t - T) < 1e-6]
        for s in ("glucose", "oxygen", "acetate"):
            v = c[c.substrate == s].frac_capped_pct; row[f"capped_{s}_pct"] = float(v.iloc[0]) if len(v) else np.nan
    rows.append(row); profs.append(b)
S = pd.DataFrame(rows).set_index("run"); ref = S.loc[REF]
S["d_core_um"] = S.R_core - ref.R_core; S["d_dead_rel_pct"] = (S.dead_pct / ref.dead_pct - 1) * 100; S["d_ringin_um"] = S.ring_in - ref.ring_in; S["d_ringout_um"] = S.ring_out - ref.ring_out
S["d_ac_pts"] = S.ac_pct_alive - ref.ac_pct_alive; S["d_rim_um"] = S.rim_width - ref.rim_width
S["within_criterion"] = (S.d_core_um.abs() <= 10) & (S.d_dead_rel_pct.abs() <= 10)
S.reset_index().to_csv(f"{UC}/results/sens_summary.csv", index=False); pd.concat(profs).to_csv(f"{UC}/results/sens_profiles.csv", index=False)
pd.set_option("display.width", 260); print(S[["dead_pct", "R_core", "rim_width", "ring_in", "ring_out", "ac_pct_alive", "d_core_um", "d_dead_rel_pct", "within_criterion"]].round(2).to_string())
