"""Field profiles and cell state at early times across dt (tests whether the dt dependence of the proliferative edge is in transport or in the cell response). usage: compare_dt_profiles.py [HOURS ...]   (dt 0.00125 = TUM-DT125, 24 h only)  -> results/tum_dt_profiles.csv
Per run and time: mean over y,z of O2 and glucose per x column, their difference to the smallest-dt run at that hour (0.00125 up to 24 h, else 0.0025) (max over x and at 100 um), number of living cells, proliferative edge, mean O2 uptake rate per growing cell."""
import sys, numpy as np, pandas as pd, scipy.io as sio, xml.etree.ElementTree as ET
sys.path.insert(0, "../../common"); from pcxml import labels
R = {0.02: "TUM-DT_dt0.02_s50", 0.01: "TUM-REF_s50", 0.005: "TUM5-REF_s50", 0.0025: "TUM-DT_dt0.0025_s50", 0.00125: "TUM-DT125_dt0.00125_s50"}; Hs = [float(a) for a in sys.argv[1:]] or [6, 12, 24]
def load(run, H):
    xf = f"runs/{run}/output{int(H):08d}.xml"; r = ET.parse(xf).getroot(); lab = labels(r); c = sio.loadmat(xf.replace(".xml", "_cells.mat"))["cells"]
    mm = sio.loadmat(xf.replace(".xml", "_microenvironment0.mat"))["multiscale_microenvironment"]; ids = {v.get("name"): int(v.get("ID")) for v in r.iter("variable") if v.get("ID") is not None and v.get("name")}
    f = pd.DataFrame({"x": mm[0], "O2": mm[4 + ids["O2"]], "glc": mm[4 + ids["D-Glucose"]]}).groupby("x").mean()
    g = c[lab.index("growth_rate")]; dead = c[lab.index("dead")]; x = c[lab.index("position[0]")]; al = dead == 0; gr = al & (g > 0.5 * 0.0137)
    return f, dict(cells=len(g), alive=int(al.sum()), n_grow=int(gr.sum()), prolif_edge=float(np.percentile(x[gr], 95)) if gr.any() else np.nan, mean_growth_living=float(g[al].mean()))
rows = []
for H in Hs:
    smallest = 0.00125 if H <= 24 else 0.0025; ref, _ = load(R[smallest], H)   # differences are against the smallest dt that exists at that hour
    for dt, run in R.items():
        if dt < smallest: continue
        f, s = load(run, H); dO = (f.O2 - ref.O2); dG = (f.glc - ref.glc); x100 = f.index[np.argmin(abs(f.index - 100))]
        rows.append(dict(t_h=H, dt=dt, **s, O2_100um=f.O2[x100], glc_100um=f.glc[x100], dO2_max_mM=dO.abs().max(), dO2_rel_max=(dO.abs() / ref.O2).max(), dglc_max_mM=dG.abs().max(), dglc_rel_max=(dG.abs() / ref.glc).max()))
d = pd.DataFrame(rows); d.to_csv("results/tum_dt_profiles.csv", index=False); pd.set_option("display.width", 250); print(d.round(4).to_string(index=False))
