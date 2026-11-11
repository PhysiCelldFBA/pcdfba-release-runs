"""Supplementary Figs. S8 (concentration profiles) and S9 (uptake stack) from one run. usage: make_supp_figs.py RUNDIR OUTDIR [HOURS=72]
Same plotting modules as make_figure.py. Corrections against the submitted S9: per-cell rates use dry mass = (1 - mean fluid fraction of the living cells) = 0.30 of the volume instead of 0.75, and the O2 row is multiplied by the cell volume like the other rows (the submitted O2 row omitted it). Zone lines = proliferative and living edges of the run."""
import sys, os, numpy as np, pandas as pd, scipy.io as sio, xml.etree.ElementTree as ET
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..", "common")); from pcxml import labels
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))   # plot_cells_spatial.py and plot_exchange_profiles.py are in this folder (copies of the PhysiCelldFBA scripts/ modules with the font-size options used here)
import types, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from PIL import Image
sys.modules["pctk"] = types.ModuleType("pctk"); sys.modules["pctk"].multicellds = types.SimpleNamespace()   # not needed: no reader is used
import figure_paths; figure_paths.ensure_paths = lambda *a, **k: None
import plot_cancer_tissue_spatial as P
from plot_cancer_tissue_spatial_gif import _render_uptake_stack
from plot_exchange_profiles import plot_exchange_profile_grid
d, out = sys.argv[1].rstrip("/"), sys.argv[2]; H = float(sys.argv[3]) if len(sys.argv) > 3 else 72.0; os.makedirs(out, exist_ok=True)
xf = f"{d}/output{int(round(H)):08d}.xml"; r = ET.parse(xf).getroot(); lab = labels(r); m = sio.loadmat(xf.replace(".xml", "_cells.mat"))["cells"]; t = float(r.find(".//current_time").text) / 60; assert abs(t - H) < 1e-6
df = pd.DataFrame({k: m[lab.index(k)] for k in lab if k.startswith("R_EX") and not k.endswith("_lb")})
for k, n in (("position[0]", "x_position"), ("position[1]", "y_position"), ("position[2]", "z_position"), ("total_volume", "total_volume"), ("growth_rate", "growth_rate"), ("dead", "dead"), ("death_rates[1]", "death_rates_1"), ("fluid_fraction", "fluid_fraction")): df[n] = m[lab.index(k)]
df["time"] = t; alive = df.dead == 0; grow = alive & (df.growth_rate > 0.5 * 0.0137)
pe, le = float(np.percentile(df.x_position[grow], 95)), float(np.percentile(df.x_position[alive], 95)); dry = float(1 - df.fluid_fraction[alive].mean())
print(f"t={t:.0f} h cells {len(df)} prolif edge {pe:.1f} living edge {le:.1f} dry fraction {dry:.3f}")
# S9: patch the module constants and the O2 row
P._SOLID = dry; P._o2_uptake = lambda row: -float(row.R_EX_o2_e) * float(row.total_volume) * P._RHO * P._SOLID / P._MIN_PER_H
specs = {s["stem"]: s for s in P._figure_specs()}; specs["oxygen_uptake"]["required"] = ["R_EX_o2_e", "total_volume"]
sty = dict(P._DEFAULT_STYLE); sty["zone_x"] = (pe, le); sty["fontsize"] = 27; sty["annotation_fontsize"] = 28   # S9 is drawn on a 20 in canvas (about 17.9 in after trimming) and printed at 4.7 in (0.75 text width): about 7 pt ticks and 7.4 pt labels as printed; the colorbar labels are long and set the limit
fr = _render_uptake_stack(df, [t], [specs[s] for s in P._DEFAULT_STACK_STEMS], sty, show_colorbar=True)
Image.fromarray(fr[0]).save(f"{out}/supp_fig_S9_uptake_stack_{int(H)}h.png")
# values for the caption
rows = []
for s in P._DEFAULT_STACK_STEMS:
    v = np.array([specs[s]["color_config"]["value_fn"](x) for x in df[alive].itertuples()]); rows.append(dict(stem=s, max=v.max(), p98=np.percentile(v, 98), mean=v.mean()))
pd.DataFrame(rows).to_csv(f"{out}/supp_fig_S9_values.csv", index=False)
# S8: concentration grid, mean over y (both z layers)
mm = sio.loadmat(xf.replace(".xml", "_microenvironment0.mat"))["multiscale_microenvironment"]
names = {v.get("name"): int(v.get("ID")) for v in r.iter("variable") if v.get("ID") is not None and v.get("name")}
dm = pd.DataFrame({"x": mm[0], "y": mm[1]})
for n, i in names.items(): dm[n] = mm[4 + i]
fb = plot_exchange_profile_grid(mode="concentration", time_point=t, df_micro=dm, zone_x=(pe, le), dpi=300, figsize=(6.0, 4.2), title_fontsize=8, tick_fontsize=7, label_fontsize=8, line_width=1.0)   # S8 printed at 0.95 text width (6.0 in)
fb.savefig(f"{out}/supp_fig_S8_concentration_grid_{int(H)}h.png", dpi=300, bbox_inches="tight"); print("wrote", out)
