"""Compare the dx/dt runs with their plateau and relate the deviation to a cap number.
usage: cap_rule.py   (from the use-case directory; /opt/anaconda3/bin/python)   writes results/ut_vs_plateau.csv, results/fig_cap_rule.png

Plateau = dx 4 um, dt 0.001 min of the same experiment (largest voxel, smallest step: the uptake cap cannot act there).
Cap number  R = n * q * m_dw * dt / (c * V_vox)  =  q * X * dt / c,
  q    specific uptake of the limiting substrate at the start (Michaelis-Menten, mmol/gDW/min),
  m_dw dry mass of one cell (g) = V_cell*(1-fluid)*density,   X = n*m_dw/V_vox the local dry biomass concentration (g/L),
  c    concentration of the substrate (mM),  V_vox = dx^3 (um3; 1 um3 = 1e-15 L),  n = cells sharing the voxel at the start.
R is the fraction of the voxel's substrate that the cells in it ask for in one time step (the time to empty the voxel is c/(q X))."""
import glob, os, re, sys, numpy as np, pandas as pd, scipy.io as sio, xml.etree.ElementTree as ET
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from lxml import etree
cfg = etree.parse("config/shipped/PhysiCell_settings.xml").getroot(); g = lambda p: float(cfg.xpath(p)[0].text)
Vmax = g("//exchange[@substrate='glucose']/Vmax") / 60.0; Km = g("//exchange[@substrate='glucose']/Km")      # mmol/gDW/min, mM
c0 = g("//variable[@name='glucose']/initial_condition"); Vc = g("//cell_definition//phenotype/volume/total"); ff = g("//phenotype/volume/fluid_fraction")
rho = g("//growth_model/cell_density") * 1e-12                                                              # g/um3
mdw = Vc * (1 - ff) * rho; q = Vmax * c0 / (Km + c0)
def R(n, dx, dt, c=c0): return n * q * mdw * dt / (c * dx ** 3 * 1e-15)
def snap(d, k=-1):
    xf = sorted(glob.glob(f"{d}/output*[0-9].xml"))[k]; r = ET.parse(xf).getroot()
    lab = {l.text: (int(l.get("index")), int(l.get("size"))) for l in r.findall(".//cell_population/custom/simplified_data/labels/label")}
    m = sio.loadmat(os.path.join(d, r.find(".//cell_population/custom/simplified_data/filename").text))["cells"]
    L = open(f"{d}/conservation_ledger.csv").read().strip().split("\n"); h = L[0].split(","); lv = dict(zip(h, [float(x) for x in L[-1].split(",")]))
    return dict(t=float(r.find(".//current_time").text), n=m.shape[1], vol=float(m[lab["total_volume"][0]].sum()), glc=-lv["glucose_cells_dF"])
rows = []
for exp, n0, ts in (("UT-DXDT-A", 1, [None]), ("UT-CONV-BOX", 8, [1, 2, 3])):
    for d in sorted(glob.glob(f"runs/{exp}_dx*")):
        dx, dt = map(float, re.findall(r"dx([\d.]+)_dt([\d.]+)", d)[0])
        for k in ts:
            s = snap(d, -1 if k is None else k); rows.append(dict(exp=exp, dx=dx, dt=dt, R=R(n0, dx, dt), **s))
df = pd.DataFrame(rows); out = []
for (exp, t), grp in df.groupby(["exp", df.t.round(0)]):
    ref = grp[(grp.dx == 4) & (grp.dt == 0.001)].iloc[0]
    for _, r in grp.iterrows(): out.append(dict(r, d_vol=100 * (r.vol / ref.vol - 1), d_glc=100 * (r.glc / ref.glc - 1) if ref.glc else np.nan, d_n=int(r.n - ref.n)))
res = pd.DataFrame(out); os.makedirs("results", exist_ok=True); res.to_csv("results/ut_vs_plateau.csv", index=False)
pd.set_option("display.width", 220)
a = res[res.exp == "UT-DXDT-A"].sort_values(["dx", "dt"])
print(f"q = {q:.4f} mmol/gDW/min, m_dw = {mdw:.3e} g, c0 = {c0} mM  ->  R = {R(1,1,1):.3f} * n * dt / dx^3 (dt in min, dx in um)")
print("\nUT-DXDT-A against the plateau (dx 4, dt 0.001), criterion 0.5% on volume and glucose taken:")
print(a[["dx", "dt", "R", "d_vol", "d_glc", "d_n"]].round(4).to_string(index=False))
b = res[(res.exp == "UT-CONV-BOX")]
print("\nUT-CONV-BOX against the plateau, worst of the volume deviations over t = 30, 60, 90 min:")
w = b.groupby(["dx", "dt"]).agg(R=("R", "first"), worst_vol=("d_vol", lambda x: x.abs().max()), worst_n=("d_n", lambda x: x.abs().max())).reset_index()
print(w.round(3).to_string(index=False))
fig, ax = plt.subplots(1, 2, figsize=(9, 3.6))
for exp, col, lab, tol in (("UT-DXDT-A", "d_vol", "1 cell, volume at 120 min", 0.5), ("UT-CONV-BOX", "d_vol", "8 cells, volume (worst of 30-90 min)", 5)):
    s = (a if exp == "UT-DXDT-A" else b.groupby(["dx", "dt"]).agg(R=("R", "first"), d_vol=("d_vol", lambda x: x.abs().max())).reset_index())
    s = s.assign(y=s[col].abs().clip(lower=1e-4)); i = 0 if exp == "UT-DXDT-A" else 1
    for dx, m in zip((0.5, 1, 2, 4), "svo^"): z = s[s.dx == dx]; ax[i].scatter(z.R, z.y, marker=m, label=f"dx {dx}")
    ax[i].axhline(tol, color="k", ls="--", lw=.8); ax[i].set_xscale("log"); ax[i].set_yscale("log"); ax[i].set_xlabel("cap number R"); ax[i].set_ylabel("|deviation from plateau| (%)"); ax[i].set_title(lab, fontsize=9)
ax[0].legend(fontsize=7); fig.tight_layout(); fig.savefig("results/fig_cap_rule.png", dpi=150)
