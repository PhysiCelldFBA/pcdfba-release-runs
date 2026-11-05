"""Tumour main figure (Fig. 5) from one run. usage: make_figure.py RUNDIR OUTPREFIX [HOURS=72]
Panel A: cells at HOURS coloured by growth rate (dead cells dark, cells with a positive necrosis rate grey); the dashed lines are the proliferative edge and the edge of the living cells of the run (analyse_tum.py definitions: 95th percentiles of x).
Panel B: net exchange rate of 15 metabolites per cell (fmol/min) against distance to the vessel, mean +- SD of the cells in 20 bins of 16 um (3-bin running mean), the same panels as the submitted figure.
Uses the plotting modules of the sample project (PhysiCelldFBA/scripts: plot_cells_spatial.py, plot_exchange_profiles.py); the cell table is built from the *_cells.mat file because pctk is not needed for that."""
import sys, os, glob, io, numpy as np, pandas as pd, scipy.io as sio, xml.etree.ElementTree as ET
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..", "common")); from pcxml import labels
SC = "/home/mponce/projects/devel/PhysiCelldFBA/scripts"; sys.path.insert(0, SC)
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.image import imread
from plot_cells_spatial import plot_cells_spatial
from plot_exchange_profiles import plot_exchange_profile_grid, DEFAULT_EXCHANGES
d, out = sys.argv[1].rstrip("/"), sys.argv[2]; H = float(sys.argv[3]) if len(sys.argv) > 3 else 72.0
xf = f"{d}/output{int(round(H * 60 / 60)):08d}.xml"
r = ET.parse(xf).getroot(); lab = labels(r); m = sio.loadmat(xf.replace(".xml", "_cells.mat"))["cells"]; t = float(r.find(".//current_time").text) / 60
assert abs(t - H) < 1e-6, (t, H)
df = pd.DataFrame({k: m[lab.index(k)] for k in lab if k.startswith("R_EX") and not k.endswith("_lb")})
for k, n in (("position[0]", "x_position"), ("position[1]", "y_position"), ("position[2]", "z_position"), ("total_volume", "total_volume"), ("growth_rate", "growth_rate"), ("dead", "dead"), ("death_rates[1]", "death_rates_1"), ("fluid_fraction", "fluid_fraction")): df[n] = m[lab.index(k)]
df["time"] = t; alive = df.dead == 0; grow = alive & (df.growth_rate > 0.5 * 0.0137)
pe, le = float(np.percentile(df.x_position[grow], 95)), float(np.percentile(df.x_position[alive], 95))
print(f"t = {t:.0f} h, cells {len(df)}, dead {int((~alive).sum())}, proliferative edge {pe:.1f} um, living edge {le:.1f} um, mean fluid fraction of living cells {df.fluid_fraction[alive].mean():.3f}")
sty = dict(X_MIN=0.0, X_MAX=320.0, Y_MIN=-55.0, Y_MAX=55.0, skip_dead=False, vmax=0.026803, vmin=0.0, vessel_style="bar", vessel_label=True, show_zones=True, zone_x=(pe, le), zone_labels=True, zone_bands=False,
           show_nuclei=True, nucleus_scale=1.0, cell_shape="blob", roughness=0.05, shape_seed=0, perspective=0.1, depth_fade=0.05, show_x_ruler=True, show_scale_bar=False, x_major=50.0, x_minor=10.0,
           xlabel="Distance from vessel (µm)", fontsize=20, annotation_fontsize=25, figsize=(16, 5.6))
fa = plot_cells_spatial(df, time_point=t, color_config={"value_fn": lambda x: x.growth_rate, "cmap": plt.cm.Greens, "label": "Growth rate (1/h)"}, **sty)
ba = io.BytesIO(); fa.savefig(ba, dpi=200, bbox_inches="tight", pad_inches=0.05, format="png"); plt.close(fa)
fb = plot_exchange_profile_grid(mode="flux", time_point=t, df_cells=df, zone_x=(pe, le), solid_fraction=float(1 - df.fluid_fraction[alive].mean()), dpi=200)
for tx in fb.texts:
    if "export" in tx.get_text(): tx.set_text("Cell net exchange rate (fmol/min)")
bb = io.BytesIO(); fb.savefig(bb, dpi=200, bbox_inches="tight", pad_inches=0.05, format="png"); plt.close(fb)
ia, ib = imread(io.BytesIO(ba.getvalue())), imread(io.BytesIO(bb.getvalue())); W = max(ia.shape[1], ib.shape[1])
fig = plt.figure(figsize=(W / 200, (ia.shape[0] + ib.shape[0]) / 200 * 1.02), dpi=200); h = ia.shape[0] + ib.shape[0]
a1 = fig.add_axes([0, ib.shape[0] / h, ia.shape[1] / W, ia.shape[0] / h]); a1.imshow(ia); a1.axis("off")
a2 = fig.add_axes([0, 0, ib.shape[1] / W, ib.shape[0] / h]); a2.imshow(ib); a2.axis("off")
fig.text(0.005, 0.995, "A", fontsize=40, fontweight="bold", va="top"); fig.text(0.005, ib.shape[0] / h + 0.002, "B", fontsize=40, fontweight="bold", va="bottom")
for ext in ("png", "pdf"): fig.savefig(f"{out}.{ext}", dpi=200)
pd.DataFrame([dict(run=os.path.basename(d), t_h=t, cells=len(df), dead=int((~alive).sum()), prolif_edge=pe, living_edge=le)]).to_csv(out + "_values.csv", index=False)
print("wrote", out)
