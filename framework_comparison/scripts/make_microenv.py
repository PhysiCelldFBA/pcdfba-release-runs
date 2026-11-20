"""Microenvironment of the consortium of Fig. 6 in PhysiCelldFBA, COMETS and BacArena, 10 to 30 h (Supplementary Figure of S12).

  results/supp_fig_benchmark_microenv.{png,pdf}
    columns: the three tools, each at 10, 20 and 30 h; rows: the cells of the two species, then glucose, H2, acetate and methane in
    the 200 um domain. Each row has one colour scale for all its panels (glucose linear, the others logarithmic).

Inputs (all in data/):
  fig6_physicelldfba/physicelldfba_CF3-NC-CTL_s50_fields.npz   extract_physicelldfba_fields.py on CF3-NC-CTL_s50 (MN5,
      pcdfba-review-runs-cfnc/uc/runs, snapshots 600, 1200, 1800 min): the configuration of Fig. 6 rerun to 30 h, used because the
      raw output of CF3-NC50_s50 is not available (PC_FIELDS overrides the path)
  fig6_comets/comets_fig6_conv_fields.npz                       extract_comets_fields.py on BM-CF-COMETS
  fig6_bacarena/bacarena_fig6_closed_nopool_{media,cells}.csv.gz extract_bacarena_fields.R on BM-CF-BACARENA
Run from framework_comparison/:  python3 scripts/make_microenv.py
"""
import os, sys
from pathlib import Path
import numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, Normalize
from matplotlib.lines import Line2D

HERE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(HERE.parent / "common")); from plot_style import substrate_cmap, apply_style
D, OUT = HERE / "data", HERE / "results"
TIMES = (10.0, 20.0, 30.0)
L = 100.0                                                   # half width of the domain (um)
TOOLS = {"PhysiCelldFBA": "#1b6ca8", "COMETS": "#d1495b", "BacArena": "#3c8d2f"}
SPC = {"CB": "#2b6cb0", "MB": "#b8460b"}                    # C. beijerinckii, M. barkeri
PACKED = {"CB": 2.13e-5, "MB": 2.37e-5}                     # g/cm2, packing density of the COMETS run (comets_cfref_10um.py)
MET = {"glucose": ("Glucose", substrate_cmap("glucose"), Normalize(0.5, 1.0)),
       "H2": ("H$_2$", substrate_cmap("h2"), LogNorm(1e-3, 1.0)),
       "acetate": ("Acetate", substrate_cmap("acetate"), LogNorm(1e-3, 0.5)),
       "methane": ("CH$_4$", substrate_cmap("methane"), LogNorm(1e-4, 0.06))}


# ----------------------------------------------------------------------------- loading: {tool: {time: dict}}
def load_pc():
    z = np.load(os.environ.get("PC_FIELDS", D / "fig6_physicelldfba/physicelldfba_CF3-NC-CTL_s50_fields.npz"))
    xs = z["xs"]; h = xs[1] - xs[0]; out = {}
    for i, t in enumerate(z["t_h"]):
        c = z[f"cells_{i}"]; live = c[:, 3] == 0
        out[round(float(t), 2)] = dict(**{k: z[k][i] for k in MET}, ext=[xs[0] - h / 2, xs[-1] + h / 2] * 2,
                                       cells={"CB": c[live & (c[:, 2] == 0), :2], "MB": c[live & (c[:, 2] == 1), :2]})
    return out


def load_comets():
    z = np.load(D / "fig6_comets/comets_fig6_conv_fields.npz"); g = float(z["grid_um"]); N = int(z["N"]); out = {}
    cen = (np.arange(N) - 1) * g + g / 2 - L; sl = slice(1, N - 1)     # index 0 and N-1: static ring outside the domain
    X, Y = np.meshgrid(cen[sl], cen[sl])
    for t in TIMES:
        i = int(np.argmin(abs(z["t_h"] - t))); boxes = {}
        for s, nm in (("CB", "c_beijerinckii"), ("MB", "m_barkeri")):
            f = np.clip(z[nm][i][sl, sl].T / (PACKED[s] * g ** 2 * 1e-8), 0, 1); ok = f >= 0.02
            boxes[s] = (X[ok], Y[ok], f[ok])
        out[t] = dict(**{k: z[k][i][sl, sl].T for k in MET}, ext=[-L, L, -L, L], boxes=boxes)
    return out


def load_bacarena():
    m = pd.read_csv(D / "fig6_bacarena/bacarena_fig6_closed_nopool_media.csv.gz")
    c = pd.read_csv(D / "fig6_bacarena/bacarena_fig6_closed_nopool_cells.csv.gz")
    n = int(m.x.max()); w = 2 * L / n; out = {}
    for t in TIMES:
        q = m[m.step_h == int(t)]; f = {}
        for k in MET:
            a = np.zeros((n, n)); a[q.y.values - 1, q.x.values - 1] = q[k].values; f[k] = a
        p = c[c.step_h == int(t)]
        cells = {s: np.c_[(p[p.species == sp].x - 0.5) * w - L, (p[p.species == sp].y - 0.5) * w - L]
                 for s, sp in (("CB", "c_beijerinckii"), ("MB", "m_barkeri"))}
        out[t] = dict(**f, ext=[-L, L, -L, L], cells=cells)
    return out


DATA = {"PhysiCelldFBA": load_pc(), "COMETS": load_comets(), "BacArena": load_bacarena()}
for tool, d in DATA.items():
    missing = [t for t in TIMES if t not in d]; assert not missing, (tool, missing)

# ----------------------------------------------------------------------------- figure
k = apply_style(6.3)
rows = ["cells"] + list(MET)
fig = plt.figure(figsize=(6.3, 4.1))
gs_l, gs_r, gs_t, gs_b, gap = 0.075, 0.885, 0.875, 0.075, 0.012
W = (gs_r - gs_l - 2 * gap * 3) / 9; H = (gs_t - gs_b) / len(rows)
for j, tool in enumerate(TOOLS):
    x0 = gs_l + j * (3 * W + 3 * gap)
    fig.text(x0 + 1.5 * W, gs_t + 0.07, tool, ha="center", color="black", fontweight="bold")
    for c_, t in enumerate(TIMES):
        d = DATA[tool][t]
        fig.text(x0 + (c_ + 0.5) * W, gs_t + 0.015, f"{t:g} h", ha="center", fontsize=6.5 * k, color="0.3")
        for r, row in enumerate(rows):
            ax = fig.add_axes([x0 + c_ * W, gs_t - (r + 1) * H, W * 0.97, H * 0.97])
            ax.set_xlim(-L, L); ax.set_ylim(-L, L); ax.set_xticks([]); ax.set_yticks([]); ax.set_aspect("equal")
            for sp in ax.spines.values(): sp.set_visible(True); sp.set_linewidth(0.4 * k); sp.set_color("0.55")
            if row == "cells":
                if "cells" in d:
                    for s, xy in d["cells"].items():
                        ax.scatter(xy[:, 0], xy[:, 1], s=(0.08 if tool == "PhysiCelldFBA" else 0.15) * k, c=SPC[s], lw=0,
                                   marker="o" if tool == "PhysiCelldFBA" else "s", rasterized=True)
                else:
                    for s, (bx, by, f) in d["boxes"].items():
                        ax.scatter(bx + (-2.5 if s == "CB" else 2.5), by, s=(0.3 + 2.2 * f) * k, c=SPC[s], marker="s", lw=0)
            else:
                lab, cmap, norm = MET[row]
                ax.imshow(np.clip(d[row], norm.vmin, None), origin="lower", extent=d["ext"], cmap=cmap, norm=norm,
                          interpolation="bilinear" if tool == "PhysiCelldFBA" else "nearest", rasterized=True)
            if j == 0 and c_ == 0:
                ax.set_ylabel({"cells": "Cells"}.get(row, MET.get(row, ("",))[0]), fontsize=7 * k, rotation=0, ha="right", va="center", labelpad=3)
for r, row in enumerate(rows[1:], start=1):
    lab, cmap, norm = MET[row]
    cax = fig.add_axes([gs_r + 0.012, gs_t - (r + 1) * H + 0.12 * H, 0.011, 0.76 * H])
    cb = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), cax=cax); cb.outline.set_linewidth(0.4 * k)
    cb.ax.tick_params(labelsize=5.5 * k, length=1.5 * k, width=0.4 * k, pad=1)
    if isinstance(norm, LogNorm): cb.ax.minorticks_off()
fig.text(gs_r + 0.075, gs_t - 3 * H, "mM", rotation=90, va="center", fontsize=6.5 * k)
fig.legend(handles=[Line2D([], [], ls="", marker="o", ms=3.5 * k, color=SPC["CB"], label="$C.~beijerinckii$"),
                    Line2D([], [], ls="", marker="o", ms=3.5 * k, color=SPC["MB"], label="$M.~barkeri$")],
           loc="lower center", ncol=2, bbox_to_anchor=(0.48, -0.005), frameon=False, handletextpad=0.2, columnspacing=1.5, fontsize=6.5 * k)
for ext in ("png", "pdf"): fig.savefig(OUT / f"supp_fig_benchmark_microenv.{ext}", dpi=400, facecolor="white")
print("written", OUT / "supp_fig_benchmark_microenv.png")
