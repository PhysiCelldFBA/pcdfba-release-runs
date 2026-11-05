"""Supplementary figure: reduced colony sensitivity set at 6 min (COL-BND 3x3 boundary grid, COL-KIN 8 Sobol sets, COL-PFBA) against the reference COL-DT_dt0.0005_s50, from results/sens_profiles.csv and results/sens_summary.csv.
A: annulus profiles (10 um) of the dead fraction (black) and of the fraction of living cells taking up acetate (red) against the distance from the centre; the reference in grey; boundary grid (glucose rows, oxygen columns; the centre is the reference).
B: the same for the 8 kinetic sets and pFBA. C: differences to the reference (core radius, dead fraction, ring edges, acetate-using fraction, rim width) and the capped fraction of glucose; pale = within the criterion fixed on 2026-11-05 (core radius within 10 um and dead fraction within 10% relative).
usage: make_figure_sens.py   -> results/fig_colony_sens.{png,pdf}"""
import os, numpy as np, pandas as pd, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.colors import TwoSlopeNorm
UC = os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."); R = UC + "/results"
P = pd.read_csv(f"{R}/sens_profiles.csv"); S = pd.read_csv(f"{R}/sens_summary.csv").set_index("run"); REF = "COL-DT_dt0.0005_s50"
D = pd.read_csv(f"{R}/sens_design.csv").set_index("run")
COL = dict(dead="#222222", ac="#cc0000", ref="#9a9a9a", ok="#2e8b57")
prof = lambda run: P[P.run == run].sort_values("r_mid")
def panel(ax, run, title, ok=None):
    r0 = prof(REF); ax.plot(r0.r_mid, r0.dead_frac, color=COL["ref"], lw=4.5, alpha=0.55); ax.plot(r0.r_mid, r0.ac_frac, color=COL["ref"], lw=4.5, alpha=0.55)
    if run != REF and run in set(P.run):
        p = prof(run); ax.plot(p.r_mid, p.dead_frac, color=COL["dead"], lw=1.6); ax.plot(p.r_mid, p.ac_frac, color=COL["ac"], lw=1.6)
    elif run == REF:
        ax.plot(r0.r_mid, r0.dead_frac, color=COL["dead"], lw=1.6); ax.plot(r0.r_mid, r0.ac_frac, color=COL["ac"], lw=1.6)
    ax.set_xlim(0, 210); ax.set_ylim(-0.03, 1.03); ax.set_xticks([0, 100, 200]); ax.set_yticks([0, 0.5, 1]); ax.tick_params(labelsize=8.5)
    [ax.spines[s].set_visible(False) for s in ("top", "right")]; ax.set_title(title, fontsize=9.5, loc="left")
    if ok: [ax.spines[s].set_color(COL["ok"]) or ax.spines[s].set_linewidth(2.2) for s in ("left", "bottom")]
fig = plt.figure(figsize=(11, 17.5)); gs = GridSpec(3, 1, height_ratios=[3.1, 3.1, 3.5], hspace=0.32, figure=fig)
ga = gs[0].subgridspec(3, 3, hspace=0.55, wspace=0.28); GL, OX = (2, 4, 16), (2.5, 5, 10)
for i, gl in enumerate(GL):
    for j, ox in enumerate(OX):
        run = REF if (gl, ox) == (4, 5) else f"COL-BND_g{gl}_o{ox}_s50"; ax = fig.add_subplot(ga[i, j])
        panel(ax, run, f"glucose {gl} mM, O$_2$ {ox} mM" + (" (reference)" if run == REF else ""), bool(S.within_criterion.get(run, False)) and run != REF)
        if j == 0: ax.set_ylabel("fraction", fontsize=9)
        if i == 2: ax.set_xlabel("distance from the centre ($\\mu$m)", fontsize=9)
fig.text(0.03, 0.965, "A", fontsize=22, fontweight="bold")
gb = gs[1].subgridspec(3, 3, hspace=0.55, wspace=0.28); kin = sorted(r for r in S.index if r.startswith("COL-KIN"))
labs = lambda r: "set " + r.split("_p")[1][:2] + "  (" + ", ".join(f"{D.loc[r, 'f_' + k]:.1f}" for k in ("vmax_glc", "vmax_o2", "vmax_ac", "km_glc", "d_glc")) + f", {D.loc[r, 'necrosis_threshold']:.1f})"
for k, run in enumerate(kin + (["COL-PFBA_s50"] if "COL-PFBA_s50" in S.index else [])):
    ax = fig.add_subplot(gb[k // 3, k % 3]); panel(ax, run, labs(run) if run.startswith("COL-KIN") else "pFBA (reference settings)", bool(S.within_criterion.get(run, False)))
    if k % 3 == 0: ax.set_ylabel("fraction", fontsize=9)
    if k >= 6: ax.set_xlabel("distance from the centre ($\\mu$m)", fontsize=9)
from matplotlib.lines import Line2D
fig.legend(handles=[Line2D([], [], color=COL["ref"], lw=4.5, alpha=0.55, label="reference"), Line2D([], [], color=COL["dead"], lw=1.6, label="dead fraction"), Line2D([], [], color=COL["ac"], lw=1.6, label="living cells taking up acetate"),
    Line2D([], [], color=COL["ok"], lw=2.2, label="within the criterion")], loc="upper center", bbox_to_anchor=(0.5, 0.688), ncol=4, fontsize=9.5, frameon=False)
fig.text(0.03, 0.655, "B", fontsize=22, fontweight="bold")
fig.text(0.5, 0.653, "kinetic sets: factors of Vmax glucose, Vmax O$_2$, Vmax acetate, Km glucose, D glucose, and the necrosis threshold value", ha="center", fontsize=8.5, color="#444")
axc = fig.add_subplot(gs[2]); order = [r for r in S.index if r != REF]
cols = [("d_core_um", "core radius\n(µm)", 10), ("d_dead_rel_pct", "dead fraction\n(% relative)", 10), ("d_rim_um", "living-rim width\n(µm)", None), ("d_ringin_um", "ring inner edge\n(µm)", None), ("d_ringout_um", "ring outer edge\n(µm)", None), ("d_ac_pts", "acetate-using\nliving cells (points)", None), ("capped_glucose_pct", "glucose capped\n(% of cells)", None), ("capped_oxygen_pct", "O$_2$ capped\n(% of cells)", None)]
CAPC = ("capped_glucose_pct", "capped_oxygen_pct")
M = np.array([[pd.to_numeric(S.loc[r, c], errors="coerce") for c, _, _ in cols] for r in order], float)
for j, (c, _, tol) in enumerate(cols):
    if c in CAPC: continue
    v = M[:, j]; sc = tol if tol else (np.nanmax(np.abs(v)) or 1)
    axc.imshow(np.where(np.isnan(v), np.nan, np.clip(v / sc, -3, 3))[:, None], extent=(j - .5, j + .5, len(order) - .5, -.5), cmap=plt.get_cmap("RdBu_r").copy(), norm=TwoSlopeNorm(0, -3, 3), aspect="auto")
for j, (c, _, _) in enumerate(cols):
    if c in CAPC: axc.imshow(np.clip(M[:, j], 0, 100)[:, None], extent=(j - .5, j + .5, len(order) - .5, -.5), cmap="Greys", vmin=0, vmax=140, aspect="auto")
for i in range(len(order)):
    for j, (c, _, tol) in enumerate(cols):
        v = M[i, j]
        if np.isnan(v): axc.text(j, i, "none", ha="center", va="center", fontsize=8, color="#555"); continue
        inside = (abs(v) <= tol) if tol else False
        axc.text(j, i, f"{v:+.0f}" if c not in CAPC else f"{v:.0f}", ha="center", va="center", fontsize=8.5, fontweight="normal" if inside else "bold", color=("#111" if v < 70 else "white") if c in CAPC else ("#111" if abs(v) < 2.2 * (tol or 1e9) else "white"))
for i, r in enumerate(order):
    if S.within_criterion[r]: axc.add_patch(plt.Rectangle((-.5, i - .5), 2, 1, fill=False, ec=COL["ok"], lw=1.8))
axc.set_xlim(-.5, len(cols) - .5); axc.set_ylim(len(order) - .5, -.5); axc.set_xticks(range(len(cols))); axc.xaxis.tick_top(); axc.set_xticklabels([l for _, l, _ in cols], fontsize=8.5)
axc.set_yticks(range(len(order))); axc.set_yticklabels([o.replace("_s50", "").replace("COL-", "") for o in order], fontsize=9); axc.tick_params(length=0); [s.set_visible(False) for s in axc.spines.values()]
axc.set_title("Difference to the reference at 6 min (green frame: both criterion columns within tolerance)", fontsize=10.5, loc="left", pad=48); fig.text(0.03, 0.34, "C", fontsize=22, fontweight="bold")
fig.subplots_adjust(left=0.11, right=0.97, top=0.95, bottom=0.03)
for e in ("png", "pdf"): fig.savefig(f"{R}/fig_colony_sens.{e}", dpi=200, bbox_inches="tight")
print("within criterion:", [r for r in order if S.within_criterion[r]])
