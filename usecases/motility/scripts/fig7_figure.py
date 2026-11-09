"""Fig. 7: A-B trajectories of the cell over the glucose field (two colorbars above), C-D glucose at the position of the cell against time with the two switching thresholds (0.6 mM to grow, 0.35 mM to return to motility), E-F position (left) and number of cells (right) over the 10 h; one legend box at the bottom. The time series are the population mean over the whole run (the first division is at about 285-295 min).
usage: fig7_figure.py [--high runs/MOT4-HIGH_s51] [--low runs/MOT3-LOW_s50] [--out results/fig7/fig7] [--scale 3]     (from the use-case directory; PYTHONPATH must contain pctk/src)
The Km scan that was panel G is now a supplementary figure (fig7_panel_g.py -> results/fig7/panel_g.png). The glucose at the cell is the linear profile between 0 at x = -150 um and the supply at x = +150 um (checked against the last snapshot)."""
import sys, os, argparse
from pathlib import Path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, matplotlib as mpl
mpl.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize
import fig7_composite as fc
ap = argparse.ArgumentParser(); ap.add_argument("--high", default="runs/MOT4-HIGH_s51"); ap.add_argument("--low", default="runs/MOT3-LOW_s50"); ap.add_argument("--out", default="results/fig7/fig7"); ap.add_argument("--scale", type=float, default=3.0); a = ap.parse_args(); S = a.scale
fc.apply_style(S)
sc = [fc.make_scenario(Path(a.high), "1.0 mM at the right wall", 1e-3, 1e-3), fc.make_scenario(Path(a.low), "0.325 mM at the right wall", 1e-3, 1e-3)]
color_col, cell_label, cmap_name = fc.metric_for_color("growth_rate"); cell_cmap = fc.truncated_cmap(cmap_name, low=0.20, high=0.92); glc_cmap = fc.resolve_glucose_cmap("yellow")
comp, gv, cv = [], [], []
for s in sc:
    snaps = fc.load_cell_snapshots(s.reader, n_frames=7, color_col=color_col, frame_spacing="linear"); Gi, vals, _, _ = fc.glucose_grid(s.reader, "last", s.domain_x, s.domain_y)
    comp.append((snaps, Gi)); gv.append(vals); cv.append(np.concatenate([sn[color_col] for sn in snaps]))
glc_norm = Normalize(float(np.nanmin(np.concatenate(gv))), float(np.nanmax(np.concatenate(gv)))); cell_norm = Normalize(0.0, max(float(np.nanmax(np.concatenate(cv))), 1e-9))
SUP = (1.0, 0.325); TH_G, TH_M = 0.6, 0.35; XL, XR = -150.0, 150.0
fig = plt.figure(figsize=(7.6 * S, 5.6 * S)); W = 0.30
gs = fig.add_gridspec(4, 1, height_ratios=[0.20, 0.95, 1.55, 0.24], left=0.125, right=0.915, bottom=0.03, top=0.975, hspace=0.34)
cb_gs = gs[0].subgridspec(1, 5, width_ratios=[1.1, 1.7, 0.5, 1.7, 1.1], wspace=0.0)
cax_g, cax_c = fig.add_subplot(cb_gs[0, 1]), fig.add_subplot(cb_gs[0, 3])
top = gs[1].subgridspec(1, 2, wspace=W); mid = gs[2].subgridspec(2, 2, wspace=W, hspace=0.16)
ax_comp = [fig.add_subplot(top[0, c]) for c in range(2)]; ax_glc = [fig.add_subplot(mid[0, c]) for c in range(2)]
ax_pos = [fig.add_subplot(mid[1, c], sharex=ax_glc[c]) for c in range(2)]; ax_leg = fig.add_subplot(gs[3]); ax_leg.axis("off")
for c, s in enumerate(sc):
    ax = ax_comp[c]; fc.draw_composite_panel(ax, s, comp[c][0], comp[c][1], s.domain_x, s.domain_y, glc_norm, glc_cmap, cell_norm, cell_cmap, color_col, cell_radius=4.0, x_scale="linear", show_time_labels=True, scale=S)
    ax.set_xlabel("Position along gradient axis (µm)")
    if c == 0: ax.set_ylabel("y position (µm)")
    else: ax.set_yticklabels([])
for cax, norm, cm, lab in ((cax_g, glc_norm, glc_cmap, "Glucose (mM)"), (cax_c, cell_norm, cell_cmap, cell_label)):
    cb = fig.colorbar(ScalarMappable(norm=norm, cmap=cm), cax=cax, orientation="horizontal"); cb.set_label(lab, labelpad=4 * S); cb.ax.xaxis.set_ticks_position("top"); cb.ax.xaxis.set_label_position("top")
    cb.ax.tick_params(labelsize=5 * S, length=2 * S, pad=2 * S); cb.outline.set_linewidth(0.5 * S)
tmax = max(float(s.df["time_h"].max()) for s in sc); ncmax = max(float(s.df["n_cells"].max()) for s in sc); twins = []; H = {}
for c, s in enumerate(sc):
    df, t = s.df, s.df["time_h"]; gloc = SUP[c] * (df["x_position"] - XL) / (XR - XL); ax = ax_glc[c]
    H["glc"], = ax.plot(t, gloc, color="#f1c232", lw=2.2 * S / 3, label="Glucose at the cell"); H["g"] = ax.axhline(TH_G, color="0.25", ls="--", lw=0.9 * S / 3, label="Switch to growth (%.2f mM)" % TH_G); H["m"] = ax.axhline(TH_M, color="0.25", ls=":", lw=1.1 * S / 3, label="Return to motility (%.2f mM)" % TH_M)
    ax.set_ylim(0, 1.02); ax.set_xlim(0, tmax); ax.tick_params(labelbottom=False); ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    if c == 0: ax.set_ylabel("Glucose (mM)", fontsize=6 * S, labelpad=6 * S)
    else: ax.tick_params(labelleft=False)
    ax = ax_pos[c]; ax2 = ax.twinx(); twins.append(ax2)
    H["x"], = ax.plot(t, df["x_position"], color="black", lw=1.6 * S / 3, label="Position of the colony"); H["n"], = ax2.plot(t, df["n_cells"], color="#1b7837", lw=1.6 * S / 3, ls="--", label="Number of cells")
    ax.set_ylim(-150, 150); ax2.set_ylim(0, max(ncmax * 1.15, 2)); ax.set_xlim(0, tmax); ax.set_xlabel("Time (h)"); ax.spines["top"].set_visible(False); ax2.spines["top"].set_visible(False)
    if c == 0: ax.set_ylabel("Position x (µm)", fontsize=6 * S, labelpad=6 * S); ax2.tick_params(labelright=False)
    else: ax.tick_params(labelleft=False); ax2.set_ylabel("Number of cells", labelpad=8 * S, fontsize=6 * S)
for ax, ax2 in zip(ax_pos, twins): ax2.set_position(ax.get_position())
for cx in (cax_g, cax_c): q = cx.get_position(); cx.set_position([q.x0, q.y0 - 0.045, q.width, q.height])
hh = [H["glc"], H["g"], H["m"], H["x"], H["n"]]; leg = ax_leg.legend(hh, [h.get_label() for h in hh], loc="center", ncol=5, frameon=True, fancybox=False, edgecolor="black", handlelength=1.6, handletextpad=0.5, columnspacing=1.5, fontsize=5.4 * S, borderpad=0.6 * S); leg.get_frame().set_linewidth(0.5 * S)
for ax, L_ in ((ax_comp[0], "A"), (ax_comp[1], "B"), (ax_glc[0], "C"), (ax_glc[1], "D"), (ax_pos[0], "E"), (ax_pos[1], "F")): ax.text({"A": -0.16, "C": -0.16, "E": -0.16}.get(L_, -0.06), 1.12 if L_ in "AB" else 1.05, L_, transform=ax.transAxes, fontsize=7 * S, fontweight="bold", va="bottom")
Path(a.out).parent.mkdir(parents=True, exist_ok=True); fig.savefig(a.out + ".png", dpi=300); fig.savefig(a.out + ".svg"); print("wrote", a.out)
