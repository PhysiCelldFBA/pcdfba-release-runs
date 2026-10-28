"""Pre-grown colony runs: radial profile per snapshot, necrotic-core radius and acetate-uptake ring.
usage: pg_analyse.py RUNDIR [RUNDIR ...]   (RUNDIR holds output*.xml / *_cells.mat / *_microenvironment0.mat). Run with /opt/anaconda3/bin/python.
Definitions (annuli of width W=10 um around the origin, alive = not dead):
  dead_frac(r)  = dead / all cells in the annulus
  ac_frac(r)    = alive cells with R_EX_ac_e < -1e-6 / alive cells
  core radius   = outer edge of the contiguous dead-dominated (dead_frac >= 0.5) region starting at the centre; 0 if none
  ring          = annuli with ac_frac >= 0.2 and >= 10 alive cells; inner and outer edge
"""
import sys, glob, os, numpy as np, pandas as pd, scipy.io as sio, xml.etree.ElementTree as ET
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..", "common")); from pcxml import labels
W, DEAD_T, AC_T, MINN = 10.0, 0.5, 0.2, 10
def analyse(d):
    rows, prof = [], []
    for s in sorted(glob.glob(f"{d}/output*[0-9].xml")):
        root = ET.parse(s).getroot(); lab = labels(root); t = float(root.find(".//current_time").text)
        m = sio.loadmat(s.replace(".xml", "_cells.mat"))["cells"]
        x, y = m[lab.index("position[0]")], m[lab.index("position[1]")]; r = np.hypot(x, y)
        dead = m[lab.index("dead")] > 0; ac = m[lab.index("R_EX_ac_e")] < -1e-6; gr = m[lab.index("growth_rate")]
        me = sio.loadmat(s.replace(".xml", "_microenvironment0.mat"))["multiscale_microenvironment"]
        mx, my = me[0], me[1]; rv = np.hypot(mx, my); names = ["oxygen", "glucose", "acetate"]
        edges = np.arange(0, r.max() + W, W); k = np.digitize(r, edges) - 1
        nA = len(edges) - 1; nd = np.zeros(nA); na = np.zeros(nA); nac = np.zeros(nA); ng = np.zeros(nA)
        for i in range(nA):
            sel = k == i; a = sel & ~dead
            nd[i] = (sel & dead).sum(); na[i] = a.sum(); nac[i] = (a & ac).sum(); ng[i] = (a & (gr > 1e-6)).sum()
        df = nd / np.maximum(nd + na, 1); af = nac / np.maximum(na, 1)
        # core: contiguous from the centre while dead_frac >= DEAD_T
        rc = 0.0
        for i in range(nA):
            if nd[i] + na[i] >= MINN and df[i] >= DEAD_T: rc = edges[i + 1]
            else: break
        ringi = [i for i in range(nA) if na[i] >= MINN and af[i] >= AC_T]
        r_in, r_out = (edges[ringi[0]], edges[ringi[-1] + 1]) if ringi else (np.nan, np.nan)
        row = dict(run=os.path.basename(d), t=t, n=len(r), dead_pct=100 * dead.mean(), ac_pct_alive=100 * ac[~dead].mean() if (~dead).any() else np.nan,
                   growing_pct_alive=100 * (gr[~dead] > 1e-6).mean() if (~dead).any() else np.nan, R_colony=np.percentile(r, 99.5), R_core=rc, ring_in=r_in, ring_out=r_out)
        for j, nm in enumerate(names):
            c = me[4 + j]; row[nm + "_centre"] = float(c[rv < 3 * 2].mean()); row[nm + "_edge"] = float(c[(rv > row["R_colony"] + 20) & (rv < row["R_colony"] + 40)].mean()) if row["R_colony"] + 40 < 480 else np.nan
        rows.append(row)
        for i in range(nA): prof.append(dict(run=os.path.basename(d), t=t, r_mid=edges[i] + W / 2, n=int(nd[i] + na[i]), dead_frac=df[i], ac_frac=af[i], growing_frac=ng[i] / max(na[i], 1)))
    return pd.DataFrame(rows), pd.DataFrame(prof)
if __name__ == "__main__":
    R, P = [], []
    for d in sys.argv[1:]:
        a, b = analyse(d.rstrip("/")); R.append(a); P.append(b)
    R, P = pd.concat(R), pd.concat(P)
    out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "results"); os.makedirs(out, exist_ok=True)
    R.to_csv(f"{out}/pg_summary.csv", index=False); P.to_csv(f"{out}/pg_profiles.csv", index=False)
    pd.set_option("display.width", 250); print(R.round(3).to_string(index=False))
