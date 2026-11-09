"""Fig. 7: A-B trajectories of the cell over the glucose field (vertical glucose colorbar to the right of B), C-D glucose at the position of the colony against time with the two switching thresholds
(0.6 mM to grow, 0.35 mM to return to motility), E-F position (left) and number of cells (right) over the 10 h; one legend box at the bottom (cell states and lines).
Cells are drawn in two colours: moving (growth rate 0) and proliferating (growth rate > 0). FOR DISPLAY ONLY the positions of overlapping cells are pushed apart (minimum centre distance 2 x radius; the centre of mass is kept),
because the cells of a young colony overlap in the simulation (mechanics time step), and a cell at a wall is drawn inside the domain; the simulation data and all curves are not changed.
usage: fig7_figure.py [--high runs/MOT4-HIGH_s51] [--low runs/MOT3-LOW_s50] [--out results/fig7/fig7]     (from the use-case directory; PYTHONPATH must contain pctk/src)
The Km scan that was panel G is a supplementary figure (fig7_panel_g.py -> results/fig7/panel_g.png). The glucose at the cell is the linear profile between 0 at x = -150 um and the supply at x = +150 um (checked against the last snapshot).
Type sizes: common/plot_style.py (printed sizes; the figure is drawn at its printed width, 6.3 in)."""
import sys, os, argparse
from pathlib import Path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..", "common"))
import numpy as np, matplotlib as mpl
mpl.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize
from matplotlib.lines import Line2D
from matplotlib.patches import Circle
import fig7_composite as fc
from plot_style import apply_style, FS, SC
ap = argparse.ArgumentParser(); ap.add_argument("--high", default="runs/MOT4-HIGH_s51"); ap.add_argument("--low", default="runs/MOT3-LOW_s50"); ap.add_argument("--out", default="results/fig7/fig7"); a = ap.parse_args()
FW, FH = 6.3, 4.6; apply_style(FW); CMOV, CPRO = "#D62728", "#1b7837"; R = 4.0
sc = [fc.make_scenario(Path(a.high), "1.0 mM at the right wall", 1e-3, 1e-3), fc.make_scenario(Path(a.low), "0.325 mM at the right wall", 1e-3, 1e-3)]
glc_cmap = fc.resolve_glucose_cmap("yellow"); comp, gv = [], []
for s in sc:
    snaps = fc.load_cell_snapshots(s.reader, n_frames=7, color_col="growth_rate", frame_spacing="linear"); Gi, vals, _, _ = fc.glucose_grid(s.reader, "last", s.domain_x, s.domain_y)
    comp.append((snaps, Gi)); gv.append(vals)
glc_norm = Normalize(float(np.nanmin(np.concatenate(gv))), float(np.nanmax(np.concatenate(gv))))
SUP = (1.0, 0.325); TH_G, TH_M = 0.6, 0.35; XL, XR = -150.0, 150.0
print("snapshot keys:", list(comp[0][0][0].keys()))

def spread(x, y, r, it=400, seed=1):
    """display only: push overlapping cells apart to a centre distance of 2 r, keep the centre of mass"""
    P = np.c_[x, y].astype(float); n = len(P)
    if n < 2: return x, y
    rng = np.random.default_rng(seed); P += rng.normal(0, 0.05, P.shape); c0 = P.mean(0)
    for _ in range(it):
        moved = False
        for i in range(n):
            for j in range(i + 1, n):
                d = P[j] - P[i]; dist = np.hypot(*d)
                if dist < 2.05 * r:
                    u = d / max(dist, 1e-9); push = (2.05 * r - dist) / 2; P[i] -= u * push; P[j] += u * push; moved = True
        if not moved: break
    P += c0 - P.mean(0); return P[:, 0], P[:, 1]

def tlabel(sn):
    for k in ("time_h", "t_h"):
        if k in sn: return "%.1f h" % sn[k]
    for k in ("time_min", "t_min", "time"):
        if k in sn: return "%.1f h" % (sn[k] / 60)
    return ""

def draw_cells(ax, snaps, Gi, s):
    ax.imshow(Gi, origin="lower", extent=[*s.domain_x, *s.domain_y], norm=glc_norm, cmap=glc_cmap, interpolation="bilinear", aspect="auto", zorder=1, alpha=0.92)
    kept = []; sep = (s.domain_x[1] - s.domain_x[0]) / len(snaps) * 0.6; trans = ax.get_xaxis_transform()
    for k, sn in enumerate(snaps):
        x, y = spread(sn["x"], sn["y"], R) if k == len(snaps) - 1 else (np.asarray(sn["x"]), np.asarray(sn["y"]))
        x = np.clip(x, s.domain_x[0] + 1.3 * R, s.domain_x[1] - 1.3 * R)   # display only: a cell at the wall is drawn inside the domain
        pro = np.asarray(sn["growth_rate"]) > 1e-6
        for i in range(len(x)):
            ax.add_patch(Circle((x[i], y[i]), R * 1.18, facecolor="white", edgecolor="none", alpha=0.92, zorder=2 + k))
            ax.add_patch(Circle((x[i], y[i]), R, facecolor=CPRO if pro[i] else CMOV, edgecolor="0.15", linewidth=0.3, zorder=3 + k))
    for k in range(len(snaps) - 1, -1, -1):
        cx = float(np.mean(snaps[k]["x"]))
        if any(abs(cx - kx) < sep for kx in kept): continue
        kept.append(cx); ax.plot([cx, cx], [1.0, 1.03], transform=trans, color="0.35", lw=0.5, clip_on=False)
        ax.text(cx, 1.05, tlabel(snaps[k]), transform=trans, ha="center", va="bottom", fontsize=FS["annot"], color="0.2")
    ax.set_xlim(*s.domain_x); ax.set_ylim(*s.domain_y); ax.set_aspect("equal", adjustable="box")

fig = plt.figure(figsize=(FW, FH)); W = 0.30
gs = fig.add_gridspec(4, 1, height_ratios=[0.95, 1.75, 0.10, 0.27], left=0.14, right=0.875, bottom=0.05, top=0.99, hspace=0.20)
top = gs[0].subgridspec(1, 2, wspace=W); mid = gs[1].subgridspec(2, 2, wspace=W, hspace=0.32)
ax_comp = [fig.add_subplot(top[0, c]) for c in range(2)]; ax_glc = [fig.add_subplot(mid[0, c]) for c in range(2)]
ax_pos = [fig.add_subplot(mid[1, c], sharex=ax_glc[c]) for c in range(2)]; ax_leg = fig.add_subplot(gs[3]); ax_leg.axis("off")
for c, s in enumerate(sc):
    ax = ax_comp[c]; draw_cells(ax, comp[c][0], comp[c][1], s); ax.set_xlabel("x (µm)")
    if c == 0: ax.set_ylabel("y (µm)")
    else: ax.tick_params(labelleft=False)
cb = fig.colorbar(ScalarMappable(norm=glc_norm, cmap=glc_cmap), cax=ax_comp[1].inset_axes([1.05, 0.0, 0.03, 1.0]), orientation="vertical"); cb.set_label("Glucose (mM)"); cb.outline.set_linewidth(0.5)
sh = [Line2D([], [], marker="o", ls="", markersize=6, markerfacecolor=col, markeredgecolor="0.15", markeredgewidth=0.3, label=lab) for col, lab in ((CMOV, "Moving cell"), (CPRO, "Proliferating cell"))]
ax2_top = max(ncmax_ := max(float(s.df["n_cells"].max()) for s in sc) * 1.15, 2); tmax = max(float(s.df["time_h"].max()) for s in sc); ncmax = max(float(s.df["n_cells"].max()) for s in sc); twins = []; H = {}
for c, s in enumerate(sc):
    df, t = s.df, s.df["time_h"]; gloc = SUP[c] * (df["x_position"] - XL) / (XR - XL); ax = ax_glc[c]
    H["glc"], = ax.plot(t, gloc, color=SC["glucose"], lw=FS["lw"] * 1.3, label="Glucose at the cell"); H["g"] = ax.axhline(TH_G, color="0.25", ls="--", lw=0.8, label="Switch to growth (%.1f mM)" % TH_G); H["m"] = ax.axhline(TH_M, color="0.25", ls=":", lw=1.0, label="Return to motility (%.2f mM)" % TH_M)
    ax.set_ylim(0, 1.02); ax.set_xlim(0, tmax); ax.tick_params(labelbottom=False)
    if c == 0: ax.set_ylabel("Glucose (mM)")
    else: ax.tick_params(labelleft=False)
    ax = ax_pos[c]; ax2 = ax.twinx(); twins.append(ax2); ax2.spines["right"].set_visible(True)
    H["x"], = ax.plot(t, df["x_position"], color="black", label="Position of the founder cell"); H["n"], = ax2.plot(t, df["n_cells"], color=CPRO, ls="--", label="Number of cells")
    ax.set_ylim(-150, 150); ax2.set_ylim(0, max(ncmax * 1.15, 2)); ax2.set_yticks(range(1, int(ax2_top) + 1, 2)); ax.set_xlim(0, tmax); ax.set_xlabel("Time (h)")
    if c == 0: ax.set_ylabel("Position x (µm)"); ax2.tick_params(labelright=False)
    else: ax.tick_params(labelleft=False); ax2.set_ylabel("Number of cells")
for ax, ax2 in zip(ax_pos, twins): ax2.set_position(ax.get_position())
hh = sh + [H["glc"], H["x"], H["g"], H["n"], H["m"]]; leg = ax_leg.legend(hh, [h.get_label() for h in hh], loc="upper center", ncol=4, frameon=True, fancybox=False, edgecolor="black", handlelength=1.8, handletextpad=0.5, columnspacing=1.4, borderpad=0.5); leg.get_frame().set_linewidth(0.5)
for ax, L_ in ((ax_comp[0], "A"), (ax_comp[1], "B"), (ax_glc[0], "C"), (ax_glc[1], "D"), (ax_pos[0], "E"), (ax_pos[1], "F")): ax.text(-0.31 if L_ in "ACE" else -0.12, 1.04, L_, transform=ax.transAxes, fontsize=FS["letter"], fontweight="bold", va="top", ha="left")
Path(a.out).parent.mkdir(parents=True, exist_ok=True); fig.savefig(a.out + ".png", dpi=600); fig.savefig(a.out + ".svg"); print("wrote", a.out)
