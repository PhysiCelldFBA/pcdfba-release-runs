"""FBA vs pFBA exchange fluxes of the tumour. usage: compare_pfba.py RUN_FBA RUN_PFBA OUTCSV [HOURS=72]
Per cell net exchange flux (mmol/gDW/h) of every exchange reaction with a nonzero flux, in 10 um bins of x (living cells): mean of each run, difference, and the summary per substrate
(maximum over bins with >= 5 cells in both runs of |difference| / (largest |mean| of the substrate), and the correlation across bins). The two runs share the seed, so the cell positions at t = 0 agree but the trajectories do not (thread scheduling)."""
import sys, os, numpy as np, pandas as pd, scipy.io as sio, xml.etree.ElementTree as ET
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..", "common")); from pcxml import labels
a, b, out = sys.argv[1].rstrip("/"), sys.argv[2].rstrip("/"), sys.argv[3]; H = float(sys.argv[4]) if len(sys.argv) > 4 else 72.0
def load(d):
    xf = f"{d}/output{int(H):08d}.xml"; lab = labels(ET.parse(xf).getroot()); m = sio.loadmat(xf.replace(".xml", "_cells.mat"))["cells"]
    df = pd.DataFrame({k: m[lab.index(k)] for k in lab if k.startswith("R_EX") and not k.endswith("_lb")}); df["x"] = m[lab.index("position[0]")]; df["alive"] = m[lab.index("dead")] == 0
    df["gr"] = m[lab.index("growth_rate")]; return df[df.alive]
A, B = load(a), load(b); ex = [c for c in A if c.startswith("R_EX") and (A[c].abs().max() > 0 or B[c].abs().max() > 0)]; e = np.arange(0, 380, 10); rows = []
for lo, hi in zip(e[:-1], e[1:]):
    sa, sb = A[(A.x >= lo) & (A.x < hi)], B[(B.x >= lo) & (B.x < hi)]
    for c in ex + ["gr"]: rows.append(dict(rxn=c, x0=lo, n_fba=len(sa), n_pfba=len(sb), mean_fba=sa[c].mean() if len(sa) else np.nan, mean_pfba=sb[c].mean() if len(sb) else np.nan))
R = pd.DataFrame(rows); R["diff"] = R.mean_pfba - R.mean_fba; R.to_csv(out + "_bins.csv", index=False); S = []
for c, d in R.groupby("rxn"):
    d = d[(d.n_fba >= 5) & (d.n_pfba >= 5)]; sc = max(d.mean_fba.abs().max(), d.mean_pfba.abs().max())
    S.append(dict(rxn=c, bins=len(d), scale=sc, max_abs_diff=d["diff"].abs().max(), max_rel_to_scale=d["diff"].abs().max() / sc if sc > 0 else np.nan, corr=np.corrcoef(d.mean_fba, d.mean_pfba)[0, 1] if len(d) > 2 and sc > 0 else np.nan,
                  n_cells_fba=len(A), n_cells_pfba=len(B)))
S = pd.DataFrame(S).sort_values("max_rel_to_scale", ascending=False); S.to_csv(out + "_summary.csv", index=False); pd.set_option("display.width", 200); print(S.round(4).to_string(index=False))
