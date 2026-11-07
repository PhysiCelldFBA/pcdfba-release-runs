"""H2 and acetate near the border of the box, per snapshot, with the living cell counts of each species.
usage: border_profile.py OUTPREFIX RUNDIR [RUNDIR ...]    writes OUTPREFIX_border.csv
Rings: ring k = voxels at k voxels from the edge of the box (ring 0 = the boundary voxels, fixed at the Dirichlet value when that is on); columns c_<sub>_ring0 / ring1 / ring2 / ring5 / ring10 (mean over the ring, mM), band5 = mean of rings 1-5, box = box mean, cmax = maximum, plus the
mean concentration in the voxels of the cells of each species (cell_<sub>_<species>) and the living cell counts nCB, nMB."""
import sys, glob, os, numpy as np, pandas as pd, scipy.io as sio, xml.etree.ElementTree as ET
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..", "common")); from pcxml import labels
VAR = {"H2": 7, "acetate": 11}; rows = []
for d in (x.rstrip("/") for x in sys.argv[2:]):
    for f in sorted(glob.glob(f"{d}/output*[0-9].xml")):
        try:
            root = ET.parse(f).getroot(); lab = labels(root); t = float(root.find(".//current_time").text)
            c = sio.loadmat(f.replace(".xml", "_cells.mat"))["cells"]; me = sio.loadmat(f.replace(".xml", "_microenvironment0.mat"))["multiscale_microenvironment"]
        except Exception as e: print("skip", f, e); continue
        g = lambda k: c[lab.index(k)]; ty, dead = g("cell_type"), g("dead") > 0; x, y = g("position[0]"), g("position[1]")
        xs, ys = np.unique(me[0]), np.unique(me[1]); ix, iy = np.searchsorted(xs, me[0]), np.searchsorted(ys, me[1]); ring = np.minimum.reduce([ix, len(xs) - 1 - ix, iy, len(ys) - 1 - iy])
        r = dict(run=os.path.basename(d), t=t, nCB=int(((ty == 0) & ~dead).sum()), nMB=int(((ty == 1) & ~dead).sum()))
        vx, vy = np.clip(np.searchsorted(xs, x), 0, len(xs) - 1), np.clip(np.searchsorted(ys, y), 0, len(ys) - 1)
        for sub, row in VAR.items():
            v = me[row]; r[f"c_{sub}_box"] = v.mean(); r[f"c_{sub}_max"] = v.max()
            for k in (0, 1, 2, 5, 10): r[f"c_{sub}_ring{k}"] = v[ring == k].mean()
            r[f"c_{sub}_band5"] = v[(ring >= 1) & (ring <= 5)].mean()
            grid = np.zeros((len(ys), len(xs))); grid[iy, ix] = v
            for k, nm in ((0, "CB"), (1, "MB")):
                s = (ty == k) & ~dead; r[f"cell_{sub}_{nm}"] = grid[vy[s], vx[s]].mean() if s.any() else np.nan
        rows.append(r)
df = pd.DataFrame(rows); df.to_csv(sys.argv[1] + "_border.csv", index=False); print(len(df), "rows ->", sys.argv[1] + "_border.csv")
