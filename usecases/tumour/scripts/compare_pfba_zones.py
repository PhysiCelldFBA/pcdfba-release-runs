"""FBA vs pFBA exchange fluxes of the tumour per zone at 72 h, against the seed/repeat noise. usage: compare_pfba_zones.py OUTCSV
Mean flux (mmol/gDW/h) of the living cells with x < 100 um (proliferative) and 100 <= x < 190 um (hypoxic) in TUM5-REF_s50, TUM5-PFBA_s50 and the noise runs TUM5-REP_s50, TUM5-SEED_s51, TUM5-SEED_s52;
diff_pct = (pFBA - REF) / |REF|; noise_pct = largest |run - REF| / |REF| over the three noise runs."""
import sys, os, numpy as np, pandas as pd, scipy.io as sio, xml.etree.ElementTree as ET
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..", "common")); from pcxml import labels
def load(r):
    xf = f"runs/TUM5-{r}/output00000072.xml"; lab = labels(ET.parse(xf).getroot()); m = sio.loadmat(xf.replace(".xml", "_cells.mat"))["cells"]
    ex = [k for k in lab if k.startswith("R_EX") and not k.endswith("_lb")]; df = pd.DataFrame({k: m[lab.index(k)] for k in ex}); df["x"] = m[lab.index("position[0]")]; return df[m[lab.index("dead")] == 0], ex
R = {r: load(r) for r in ("REF_s50", "REP_s50", "SEED_s51", "SEED_s52", "PFBA_s50")}; runs = {k: v[0] for k, v in R.items()}; ex = R["REF_s50"][1]; rows = []
for c in ex:
    for z, (a, b) in {"proliferative (<100 um)": (0, 100), "hypoxic (100-190 um)": (100, 190)}.items():
        v = {r: d[(d.x >= a) & (d.x < b)][c].mean() for r, d in runs.items()}; ref = v["REF_s50"]
        if max(abs(x) for x in v.values()) < 1e-6: continue
        rows.append(dict(rxn=c[5:-2], zone=z, ref=ref, pfba=v["PFBA_s50"], diff_pct=100 * (v["PFBA_s50"] - ref) / abs(ref) if ref else np.nan, noise_pct=100 * max(abs(v[r] - ref) for r in ("REP_s50", "SEED_s51", "SEED_s52")) / abs(ref) if ref else np.nan))
T = pd.DataFrame(rows); T.to_csv(sys.argv[1], index=False); pd.set_option("display.width", 200, "display.max_rows", 100); print(T.round(3).to_string(index=False))
k = T.dropna(); print(f"{(k.diff_pct.abs() > k.noise_pct).sum()} of {len(k)} zone means differ by more than the largest seed/repeat difference; {(k.diff_pct.abs() > 20).sum()} by more than 20%")
