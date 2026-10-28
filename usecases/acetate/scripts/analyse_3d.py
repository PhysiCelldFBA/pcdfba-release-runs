"""Regime times and spreads for the acetate 3D batch. Usage: analyse_3d.py [run names]; /opt/anaconda3/bin/python."""
import sys, os, glob, numpy as np, pandas as pd, scipy.io as sio, xml.etree.ElementTree as ET
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..", "common")); from pcxml import labels
HERE = os.path.dirname(os.path.abspath(__file__)); R = os.path.join(HERE, "..", "runs"); OUT = os.path.join(HERE, "..", "results")
names = sys.argv[1:] or sorted(n for n in os.listdir(R) if n.startswith("ACE-"))
rows = []
for r in names:
    xs = [x for x in sorted(glob.glob(f"{R}/{r}/output0*.xml")) if "_" not in os.path.basename(x)]
    T, N, AC, GL, O2 = [], [], [], [], []
    for s in xs:
        root = ET.parse(s).getroot(); L = labels(root); m = sio.loadmat(s.replace(".xml", "_cells.mat"))["cells"]
        alive = m[L.index("dead")] == 0
        T.append(float(root.find(".//current_time").text)); N.append(int(alive.sum()))
        for lst, k in ((AC, "R_EX_ac_e"), (GL, "R_EX_glc__D_e"), (O2, "R_EX_o2_e")):
            k = k if k in L else {"R_EX_ac_e": "acetate_flux", "R_EX_glc__D_e": "glucose_flux", "R_EX_o2_e": "oxygen_flux"}[k]   # batch 1 binary wrote the old names
            lst.append(float(np.median(m[L.index(k)][alive])) if alive.any() else np.nan)
    T, N, AC, GL, O2 = map(np.array, (T, N, AC, GL, O2))
    sec = T[AC > 0.5]; reuse = T[AC < -0.5]
    rows.append(dict(run=r.replace("ac3d_", ""), snaps=len(T), n_end=int(N[-1]), n_half=int(N[len(N)//2]),
                     ferm_start=sec.min() if len(sec) else np.nan, ferm_end=sec.max() if len(sec) else np.nan,
                     reuse_start=reuse.min() if len(reuse) else np.nan, ac_peak=AC.max(), ac_min=AC.min(), glc_min=GL.min(), o2_min=O2.min()))
D = pd.DataFrame(rows); D.to_csv(f"{OUT}/regimes_{names[0].split('_')[0]}.csv" if len({n.split('_')[0] for n in names}) == 1 else f"{OUT}/regimes_3d.csv", index=False)
pd.set_option("display.width", 250); print(D.round(2).to_string())
