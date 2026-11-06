"""Undetermined exchanges of the acetate switch over time from the fva_replay tables of ACE-REF_s50 (results/fva/ACE-REF_s50_snap*_cells.csv; made with fva_replay.py ... --snap N --n 200). Writes results/fva_scan_ACE-REF.csv. usage: fva_scan.py"""
import os, re, glob, pandas as pd
UC = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); rows = []
for f in sorted(glob.glob(f"{UC}/results/fva/ACE-REF_s50_snap*_cells.csv"), key=lambda x: int(re.search(r"snap(\d+)", x).group(1))):
    d = pd.read_csv(f); d = d[d.fraction == 0.99999]; r = dict(t_min=d.t.iloc[0], cells=d.cell.nunique(), infeasible=int(d.drop_duplicates("cell").obj_fba.isna().sum()), infeasible_zero_growth=int(d.drop_duplicates("cell").query("obj_fba != obj_fba and growth_rec == 0").shape[0]))
    for ex in ("R_EX_ac_e", "R_EX_co2_e", "R_EX_o2_e", "R_EX_glc__D_e"):
        g = d[(d.exchange == ex) & d.obj_fba.notna()]; k = ex[5:-2]; r[k + "_undet_pct"] = round(100 * (g.width > 1e-2).mean(), 1) if len(g) else None; r[k + "_max_width"] = round(g.width.max(), 4) if len(g) else None
    rows.append(r)
pd.DataFrame(rows).to_csv(f"{UC}/results/fva_scan_ACE-REF.csv", index=False); print(pd.DataFrame(rows).to_string(index=False))
