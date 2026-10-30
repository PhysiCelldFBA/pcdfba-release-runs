"""Full vs quarter colony growth: per snapshot cell count (quarter x4), dead %, R99, axis/diagonal radius. usage: qtr_compare.py OUT.csv RUNDIR ..."""
import sys, glob, os, numpy as np, pandas as pd, scipy.io as sio, xml.etree.ElementTree as ET
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..", "common")); from pcxml import labels
rows = []
for d in sys.argv[2:]:
    qt = "_qtr_" in d
    for s in sorted(glob.glob(f"{d}/output*[0-9].xml")):
        root = ET.parse(s).getroot(); lab = labels(root); t = float(root.find(".//current_time").text)
        m = sio.loadmat(s.replace(".xml", "_cells.mat"))["cells"]
        x, y = m[lab.index("position[0]")], m[lab.index("position[1]")]; r = np.hypot(x, y); dead = m[lab.index("dead")] > 0
        ang = np.degrees(np.arctan2(np.abs(y), np.abs(x))); ax = (ang < 10) | (ang > 80); dg = (ang > 35) & (ang < 55)
        q = lambda v, p: float(np.percentile(v, p)) if len(v) else np.nan
        rows.append(dict(run=os.path.basename(d), t=t, n=len(r), n_quadrant=int(((x >= 0) & (y >= 0)).sum()), dead_pct=100 * dead.mean(), R99=q(r, 99), R_axis=q(r[ax], 99), R_diag=q(r[dg], 99)))
pd.DataFrame(rows).to_csv(sys.argv[1], index=False)
