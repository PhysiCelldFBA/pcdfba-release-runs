"""Summaries for UT-BASE, UT-CONV-BOX and UT-THREADS (the existing check_*.py scripts cover DXDT, VOX and ORDER).
usage: analyse_ut.py   (from the use-case directory; writes results/ut_convbox.csv, results/ut_base_threads.csv)
Observables at the last snapshot (CONV-BOX: at 30, 60, 90 and 120 min): total cell volume (um3), cell count, mean growth rate, and the remaining amount (conc x voxel volume) of glucose, oxygen, acetate.
CONV-BOX is compared at each time with the plateau run (dx 4, dt 0.001, where the uptake cap cannot act); the criterion is 5% in total volume and cell count at 30, 60 and 90 min (fixed in advance)."""
import glob, os, re, sys, numpy as np, pandas as pd, scipy.io as sio, xml.etree.ElementTree as ET
def last(d, t=None):
    xs = sorted(glob.glob(f"{d}/output*[0-9].xml"))
    if t is not None: xs = [x for x in xs if abs(float(ET.parse(x).getroot().find(".//current_time").text) - t) < 1e-6]
    if not xs: return None
    xf = xs[-1]; r = ET.parse(xf).getroot()
    lab = {l.text: (int(l.get("index")), int(l.get("size"))) for l in r.findall(".//cell_population/custom/simplified_data/labels/label")}
    m = sio.loadmat(os.path.join(d, r.find(".//cell_population/custom/simplified_data/filename").text))["cells"]
    names = [x.get("name") for x in r.findall(".//microenvironment/domain/variables/variable")]
    me = sio.loadmat(os.path.join(d, r.find(".//microenvironment/domain/data/filename").text))["multiscale_microenvironment"]
    o = dict(t=float(r.find(".//current_time").text), n=m.shape[1], vol=float(m[lab["total_volume"][0]].sum()), mu=float(m[lab["growth_rate"][0]].mean()))
    for i, s in enumerate(names):
        if s in ("glucose", "oxygen", "acetate"): o[s] = float((me[3] * me[4 + i]).sum())
    return o
rows = []
for d in sorted(glob.glob("runs/UT-*")):
    n = os.path.basename(d); m = re.match(r"UT-(THREADS|BASE)_(.*)", n)
    if m and glob.glob(d + "/output*.xml"): rows.append(dict(run=n, exp=m.group(1), tag=m.group(2), **last(d)))
df = pd.DataFrame(rows); os.makedirs("results", exist_ok=True); pd.set_option("display.width", 250); pd.set_option("display.float_format", lambda x: f"{x:.6g}")
TIMES = (30, 60, 90, 120); PLATEAU = (4.0, 0.001); CRIT = 5.0
cr = []
for d in sorted(glob.glob("runs/UT-CONV-BOX_*")):
    n = os.path.basename(d); tag = n.split("_", 1)[1]
    for t in TIMES:
        o = last(d, t)
        if o: cr.append(dict(run=n, tag=tag, **o))
cb = pd.DataFrame(cr)
if len(cb):
    cb["dx"] = cb.tag.str.extract(r"dx([\d.]+)")[0].astype(float); cb["dtv"] = cb.tag.str.extract(r"dt([\d.]+)")[0].astype(float)
    pl = cb[(cb.dx == PLATEAU[0]) & (cb.dtv == PLATEAU[1])].set_index("t")
    for q in ("vol", "n", "mu", "glucose", "oxygen", "acetate"): cb["d_" + q + "_%"] = [100 * (r[q] / pl.loc[r.t, q] - 1) if r.t in pl.index and pl.loc[r.t, q] else np.nan for _, r in cb.iterrows()]
    cb = cb.sort_values(["t", "dx", "dtv"]); cb.to_csv("results/ut_convbox.csv", index=False)
    print(f"CONV-BOX plateau dx {PLATEAU[0]} dt {PLATEAU[1]}; criterion {CRIT}% in volume and cell count at 30, 60, 90 min")
    for t in TIMES:
        c = cb[cb.t == t]; print(f"\nt = {t} min"); print(c[["dx", "dtv", "n", "vol", "d_vol_%", "d_n_%", "d_glucose_%"]].to_string(index=False))
    c = cb[cb.t.isin((30, 60, 90))]; bad = c[(c["d_vol_%"].abs() > CRIT) | (c["d_n_%"].abs() > CRIT)]
    print(f"\n{len(c)} comparisons at 30/60/90 min, {len(bad)} outside {CRIT}%:"); print(bad[["dx", "dtv", "t", "d_vol_%", "d_n_%"]].to_string(index=False) if len(bad) else "none")
if df.empty: df = pd.DataFrame(columns=["run", "exp", "tag", "t", "n", "vol", "mu", "glucose", "oxygen", "acetate"])
bt = df[df.exp.isin(["BASE", "THREADS"])]; bt.to_csv("results/ut_base_threads.csv", index=False) if len(bt) else None
print("\nBASE (5 identical runs):"); b = bt[bt.exp == "BASE"]; print(b[["run", "n", "vol", "mu", "glucose", "oxygen"]].to_string(index=False))
print("spread (max-min)/mean %:", {q: round(100 * (b[q].max() - b[q].min()) / b[q].mean(), 6) for q in ("vol", "mu", "glucose", "oxygen")})
print("\nTHREADS:"); print(bt[bt.exp == "THREADS"][["run", "n", "vol", "mu", "glucose", "oxygen"]].to_string(index=False))
