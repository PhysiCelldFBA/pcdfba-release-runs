"""Summaries for UT-BASE, UT-CONV-BOX and UT-THREADS (the existing check_*.py scripts cover DXDT, VOX and ORDER).
usage: analyse_ut.py   (from the use-case directory; writes results/ut_convbox.csv, results/ut_base_threads.csv)
Observables at the last snapshot: total cell volume (um3), cell count, mean growth rate, and the remaining amount (conc x voxel volume) of glucose, oxygen, acetate.
CONV-BOX is compared with its finest run (smallest dx, smallest dt); the criterion is 5% (fixed in advance)."""
import glob, os, re, sys, numpy as np, pandas as pd, scipy.io as sio, xml.etree.ElementTree as ET
def last(d):
    xf = sorted(glob.glob(f"{d}/output*[0-9].xml"))[-1]; r = ET.parse(xf).getroot()
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
    n = os.path.basename(d); m = re.match(r"UT-(CONV-BOX|THREADS|BASE)_(.*)", n)
    if m and glob.glob(d + "/output*.xml"): rows.append(dict(run=n, exp=m.group(1), tag=m.group(2), **last(d)))
df = pd.DataFrame(rows); os.makedirs("results", exist_ok=True); pd.set_option("display.width", 250); pd.set_option("display.float_format", lambda x: f"{x:.6g}")
cb = df[df.exp == "CONV-BOX"].copy(); cb["dx"] = cb.tag.str.extract(r"dx([\d.]+)")[0].astype(float); cb["dtv"] = cb.tag.str.extract(r"dt([\d.]+)")[0].astype(float)
ref = cb.sort_values(["dx", "dtv"]).iloc[0]
for q in ("vol", "n", "mu", "glucose", "oxygen", "acetate"): cb["d_" + q + "_%"] = 100 * (cb[q] / ref[q] - 1)
cb = cb.sort_values(["dx", "dtv"]); cb.to_csv("results/ut_convbox.csv", index=False)
print(f"CONV-BOX reference dx {ref.dx} dt {ref.dtv}: vol {ref.vol:.5g}, n {int(ref.n)}, mu {ref.mu:.5g}")
print(cb[["dx", "dtv", "n", "vol", "d_vol_%", "d_n_%", "d_mu_%", "d_glucose_%", "d_oxygen_%"]].to_string(index=False))
bt = df[df.exp.isin(["BASE", "THREADS"])]; bt.to_csv("results/ut_base_threads.csv", index=False)
print("\nBASE (5 identical runs):"); b = bt[bt.exp == "BASE"]; print(b[["run", "n", "vol", "mu", "glucose", "oxygen"]].to_string(index=False))
print("spread (max-min)/mean %:", {q: round(100 * (b[q].max() - b[q].min()) / b[q].mean(), 6) for q in ("vol", "mu", "glucose", "oxygen")})
print("\nTHREADS:"); print(bt[bt.exp == "THREADS"][["run", "n", "vol", "mu", "glucose", "oxygen"]].to_string(index=False))
