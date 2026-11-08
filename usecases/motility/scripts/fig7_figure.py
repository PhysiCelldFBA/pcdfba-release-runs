"""Fig. 7 in one layout (panels A-G): A-B trajectories over the glucose field with the two colorbars above them (horizontal), C-D migration speed and growth rate with their legend centered above the pair, E-F fluxes with their legend centered above the pair, G time to reach x = 25 um against Km.
usage: fig7_figure.py [--high runs/MOT4-HIGH_s51] [--low runs/MOT3-LOW_s50] [--tmax-min 280] [--out results/fig7/fig7] [--scale 3]     (from the use-case directory; PYTHONPATH must contain pctk/src)
The time series (C-F) stop at --tmax-min (before the first division, so that they are those of the founder cell); the trajectories (A-B) keep all snapshots."""
import sys, os, argparse
from pathlib import Path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, matplotlib as mpl
mpl.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize
import fig7_composite as fc, fig7_panel_g as pg
ap = argparse.ArgumentParser(); ap.add_argument("--high", default="runs/MOT4-HIGH_s51"); ap.add_argument("--low", default="runs/MOT3-LOW_s50"); ap.add_argument("--tmax-min", type=float, default=280.0)
ap.add_argument("--out", default="results/fig7/fig7"); ap.add_argument("--scale", type=float, default=3.0); a = ap.parse_args(); S = a.scale
fc.apply_style(S)
sc = [fc.make_scenario(Path(a.high), "1.0 mM at the right wall", 1e-3, 1e-3), fc.make_scenario(Path(a.low), "0.325 mM at the right wall", 1e-3, 1e-3)]
for s in sc: s.df = s.df[s.df["time_h"] * 60.0 <= a.tmax_min].reset_index(drop=True)
color_col, cell_label, cmap_name = fc.metric_for_color("growth_rate"); cell_cmap = fc.truncated_cmap(cmap_name, low=0.20, high=0.92); glc_cmap = fc.resolve_glucose_cmap("yellow")
comp, gv, cv = [], [], []
for s in sc:
    snaps = fc.load_cell_snapshots(s.reader, n_frames=7, color_col=color_col, frame_spacing="linear"); Gi, vals, _, _ = fc.glucose_grid(s.reader, "last", s.domain_x, s.domain_y)
    comp.append((snaps, Gi)); gv.append(vals); cv.append(np.concatenate([sn[color_col] for sn in snaps]))
glc_norm = Normalize(float(np.nanmin(np.concatenate(gv))), float(np.nanmax(np.concatenate(gv)))); cell_norm = Normalize(0.0, max(float(np.nanmax(np.concatenate(cv))), 1e-9))
fig = plt.figure(figsize=(7.6 * S, 9.0 * S))
gs = fig.add_gridspec(8, 2, height_ratios=[0.20, 0.95, 0.26, 1.0, 0.26, 1.0, 0.22, 1.25], left=0.11, right=0.93, bottom=0.06, top=0.97, wspace=0.55, hspace=0.34)
cb_gs = gs[0, :].subgridspec(1, 5, width_ratios=[1.1, 1.7, 0.5, 1.7, 1.1], wspace=0.0)
cax_g, cax_c = fig.add_subplot(cb_gs[0, 1]), fig.add_subplot(cb_gs[0, 3])
ax_comp = [fig.add_subplot(gs[1, c]) for c in range(2)]; ax_leg1 = fig.add_subplot(gs[2, :]); ax_line = [fig.add_subplot(gs[3, c]) for c in range(2)]
ax_leg2 = fig.add_subplot(gs[4, :]); ax_flux = [fig.add_subplot(gs[5, c]) for c in range(2)]; ax_g = fig.add_subplot(gs[7, :].subgridspec(1, 3, width_ratios=[0.55, 1.9, 0.55])[0, 1])
for ax in (ax_leg1, ax_leg2): ax.axis("off")
for c, s in enumerate(sc):
    ax = ax_comp[c]; fc.draw_composite_panel(ax, s, comp[c][0], comp[c][1], s.domain_x, s.domain_y, glc_norm, glc_cmap, cell_norm, cell_cmap, color_col, cell_radius=4.0, x_scale="linear", show_time_labels=True, scale=S)
    ax.set_xlabel("Position along gradient axis (µm)")
    if c == 0: ax.set_ylabel("y position (µm)")
    else: ax.set_yticklabels([])
for cax, norm, cm, lab in ((cax_g, glc_norm, glc_cmap, "Glucose (mM)"), (cax_c, cell_norm, cell_cmap, cell_label)):
    cb = fig.colorbar(ScalarMappable(norm=norm, cmap=cm), cax=cax, orientation="horizontal"); cb.set_label(lab, labelpad=4 * S); cb.ax.xaxis.set_ticks_position("top"); cb.ax.xaxis.set_label_position("top")
    cb.ax.tick_params(labelsize=5 * S, length=2 * S, pad=2 * S); cb.outline.set_linewidth(0.5 * S)
ms = max(float(s.df["migration_speed"].max()) for s in sc); mg = max(float(s.df["growth_rate"].max()) for s in sc); mf = max(float(max(s.df["atp_flux"].max(), s.df["motility_atp_flux"].max(), s.df["glucose_flux"].abs().max())) for s in sc); mt = max(float(s.df["time_h"].max()) for s in sc)
fcol = {"atp": "#3366AA", "mot": "#88CCEE", "glc": "#f1c232"}; twins = []; h1 = None
for c, s in enumerate(sc):
    df, t = s.df, s.df["time_h"]; ax = ax_line[c]; ax2 = ax.twinx(); twins.append(ax2)
    l1, = ax.plot(t, df["migration_speed"], color="black", ls="-", label="Migration speed"); l2, = ax2.plot(t, df["growth_rate"], color="black", ls="--", label="Growth rate")
    ax.set_ylim(0, max(0.05, ms * 1.08)); ax2.set_ylim(0, max(0.05, mg * 1.12)); ax.set_xlim(0, mt); ax.set_xticklabels([]); ax.set_ylabel("Migration speed\n(µm/min)"); ax2.set_ylabel("Growth rate\n(h$^{-1}$)", labelpad=10 * S)
    ax.spines["top"].set_visible(False); ax2.spines["top"].set_visible(False); h1 = h1 or [l1, l2]
    ax = ax_flux[c]; ax.plot(t, df["atp_flux"], color=fcol["atp"], label="Total ATP flux"); ax.plot(t, df["motility_atp_flux"], color=fcol["mot"], label="Motility ATP flux"); ax.plot(t, df["glucose_flux"].abs(), color=fcol["glc"], ls="--", label="|Glucose flux|")
    ax.set_ylim(0, max(1.0, mf * 1.12)); ax.set_xlim(0, mt); ax.set_xlabel("Time (h)"); ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    if c == 0: ax.set_ylabel("Flux\n(mmol·gDW$^{-1}$·h$^{-1}$)")
    else: ax.set_yticklabels([])
for ax, ax2 in zip(ax_line, twins): ax2.set_position(ax.get_position())
ax_leg1.legend(h1, [h.get_label() for h in h1], loc="center", ncol=2, frameon=False, handlelength=2.4 * S, columnspacing=2.0 * S)
hf, lf = ax_flux[-1].get_legend_handles_labels(); ax_leg2.legend(hf, lf, loc="center", ncol=3, frameon=False, handlelength=2.4 * S, columnspacing=2.0 * S)
pg.draw_panel_g(ax_g, k=2.2 * S / 3)
for ax, L in ((ax_comp[0], "A"), (ax_comp[1], "B"), (ax_line[0], "C"), (ax_line[1], "D"), (ax_flux[0], "E"), (ax_flux[1], "F"), (ax_g, "G")): ax.text(-0.16 if L in "ACEG" else -0.20, 1.12 if L in "AB" else 1.05, L, transform=ax.transAxes, fontsize=7 * S, fontweight="bold", va="bottom")
Path(a.out).parent.mkdir(parents=True, exist_ok=True); fig.savefig(a.out + ".png", dpi=300); fig.savefig(a.out + ".svg"); print("wrote", a.out)
