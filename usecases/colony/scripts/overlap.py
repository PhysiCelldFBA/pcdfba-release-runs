"""Cell overlap per snapshot: spheres of radius (3V/4pi)^(1/3) at the cell positions; two cells overlap when d < r_i + r_j.
usage: overlap.py RUNDIR   ()   writes results/overlap_<run>.csv
Reports n, % of cells with >=1 overlap, mean and max relative depth (r_i+r_j-d)/(r_i+r_j) over overlapping pairs, mean neighbours within 2(r_i+r_j)."""
import sys, glob, os, numpy as np, pandas as pd, scipy.io as sio, xml.etree.ElementTree as ET
from scipy.spatial import cKDTree
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..", "common")); from pcxml import labels
d = sys.argv[1].rstrip("/"); rows = []
for s in sorted(glob.glob(f"{d}/output*[0-9].xml")):
    root = ET.parse(s).getroot(); lab = labels(root); t = float(root.find(".//current_time").text)
    m = sio.loadmat(s.replace(".xml", "_cells.mat"))["cells"]
    p = np.c_[m[lab.index("position[0]")], m[lab.index("position[1]")], m[lab.index("position[2]")]]
    r = (3 * m[lab.index("total_volume")] / (4 * np.pi)) ** (1 / 3); n = len(r)
    tree = cKDTree(p); pairs = tree.query_pairs(2 * r.max() * 1.0, output_type="ndarray")
    if len(pairs):
        dist = np.linalg.norm(p[pairs[:, 0]] - p[pairs[:, 1]], axis=1); rs = r[pairs[:, 0]] + r[pairs[:, 1]]
        ov = dist < rs; depth = (rs - dist)[ov] / rs[ov]; ids = np.unique(pairs[ov].ravel())
        near = np.bincount(pairs[dist < 2 * rs].ravel(), minlength=n).mean()
    else: ov = np.array([], bool); depth = np.array([]); ids = np.array([], int); near = 0
    rows.append(dict(t=t, n=n, r_mean_um=r.mean(), pct_cells_overlapping=100 * len(ids) / n, mean_depth=depth.mean() if len(depth) else 0,
                     max_depth=depth.max() if len(depth) else 0, pairs_overlapping=int(ov.sum()), mean_neighbours_2x=near))
out = pd.DataFrame(rows); R = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "results"); os.makedirs(R, exist_ok=True); out.to_csv(os.path.join(R, f"overlap_{os.path.basename(d)}.csv"), index=False); pd.set_option("display.width", 220); print(out.round(3).to_string(index=False))
