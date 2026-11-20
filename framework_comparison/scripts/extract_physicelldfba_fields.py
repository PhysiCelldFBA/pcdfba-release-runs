"""Snapshot fields and cells of a PhysiCelldFBA consortium run for the microenvironment figures of S12 (make_microenv.py).
usage: python3 extract_physicelldfba_fields.py <run dir> <out.npz> [T_MIN ...]      (default 600 1800 3000)
For each requested time, the snapshot output*.xml closest to it. Writes
  t_h                        (T,)          time of each snapshot
  xs, ys                     (nx,), (ny,)  voxel centres (um)
  glucose, H2, acetate, methane (T, ny, nx) concentration per voxel (mM; rows of microenvironment0.mat 4, 7, 11, 5 as in
                                           usecases/crossfeeding/scripts/make_supp_snapshots.py), indexed [t, y, x]
  cells_<i>                  (n_i, 5)      x, y (um), cell_type (0 C. beijerinckii, 1 M. barkeri), dead, total volume (um3)
The domain is one voxel thick, so each field is the 2D slab."""
import sys, glob, os, numpy as np, scipy.io as sio, xml.etree.ElementTree as ET
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "common")); from pcxml import labels

run, out = sys.argv[1].rstrip("/"), sys.argv[2]
TIMES = [float(v) for v in sys.argv[3:]] or [600.0, 1800.0, 3000.0]
VAR = {"glucose": 4, "methane": 5, "H2": 7, "acetate": 11}
files = sorted(glob.glob(f"{run}/output*[0-9].xml"))
tof = lambda f: float(ET.parse(f).getroot().find(".//current_time").text)
tt = np.array([tof(f) for f in files])
res = {k: [] for k in VAR}; res["t_h"] = []
for i, T in enumerate(TIMES):
    f = files[int(np.argmin(abs(tt - T)))]; root = ET.parse(f).getroot(); lab = labels(root)
    me = sio.loadmat(f.replace(".xml", "_microenvironment0.mat"))["multiscale_microenvironment"]
    xs, ys = np.unique(me[0]), np.unique(me[1]); assert len(np.unique(me[2])) == 1, "expected a one-voxel-thick domain"
    ix, iy = np.searchsorted(xs, me[0]), np.searchsorted(ys, me[1])
    for k, r in VAR.items():
        g = np.zeros((len(ys), len(xs))); g[iy, ix] = me[r]; res[k].append(g)
    c = sio.loadmat(f.replace(".xml", "_cells.mat"))["cells"]; g = lambda k: c[lab.index(k)]
    res[f"cells_{i}"] = np.c_[g("position[0]"), g("position[1]"), g("cell_type"), g("dead"), g("total_volume")]
    res["t_h"].append(tof(f) / 60.0); print(os.path.basename(f), tof(f), "min", c.shape[1], "cells")
np.savez_compressed(out, xs=xs, ys=ys, **{k: np.array(v) for k, v in res.items() if not k.startswith("cells_")},
                    **{k: v for k, v in res.items() if k.startswith("cells_")})
