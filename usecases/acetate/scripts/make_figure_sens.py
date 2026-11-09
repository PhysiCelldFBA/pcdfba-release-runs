"""Supplementary figure: the 26 ACE-SENS runs (each parameter x0.5 and x2, seed 50) against ACE-REF, 600 min, from results/ace_series.csv and results/ace_summary.csv.
A: acetate and glucose box-mean concentrations against time, one panel per parameter (reference grey, x0.5 dashed, x2 solid; a green frame marks a parameter whose two runs are within 5% of ACE-REF on all six observables).
B: relative change of the six observables of the criterion (fermentation start and end, acetate-reuse start, peak acetate, final biomass, acetate floor) for every run; cells within 5% are pale, 'none' = the event did not occur.
usage: make_figure_sens.py   -> results/fig_acetate_sens.{png,pdf}, results/ace_sens_rel.csv"""
import os, numpy as np, pandas as pd, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.colors import TwoSlopeNorm
import sys; sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..", "common")); from plot_style import SC, apply_style, FS
apply_style(6.3)
UC = os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."); R = UC + "/results"
S = pd.read_csv(f"{R}/ace_series.csv"); Q = pd.read_csv(f"{R}/ace_summary.csv").set_index("run")
OBS = [("ferm_start", "Fermentation\nstart"), ("ferm_end", "Fermentation\nend"), ("reuse_start", "Acetate-reuse\nstart"), ("ac_peak", "Peak\nacetate"), ("biomass_end_gL", "Final\nbiomass"), ("ac_floor", "Acetate\nfloor")]
PAR = [("D_glucose", "$D$ glucose"), ("D_oxygen", "$D$ oxygen"), ("D_acetate", "$D$ acetate"), ("Km_glucose", "$K_m$ glucose"), ("Km_oxygen", "$K_m$ oxygen"), ("Km_acetate", "$K_m$ acetate"),
       ("Vmax_glucose", "$V_{max}$ glucose"), ("Vmax_oxygen", "$V_{max}$ oxygen"), ("Vmax_acetate", "$V_{max}$ acetate"), ("glc0", "Initial glucose"), ("ncells", "Cell number (fixed box)"), ("vol", "Cells and volume\n(x0.44, x2)"), ("pulse", "O$_2$ pulse size"), ("tpulse", "O$_2$ pulse time")]
FAC = lambda p: ("0.44", "2") if p == "vol" else ("0.5", "2")   # the volume pair is x0.44 (9 cells in 36x36x16 um) and x2 (40 cells in 36x36x72 um)
rid = lambda p, f: f"ACE-SENS_{p}_x{f}_s50"; ref = Q.loc["ACE-REF_s50", [o for o, _ in OBS]]
rel = pd.DataFrame({(p, i): (Q.loc[rid(p, f), [o for o, _ in OBS]] / ref - 1) * 100 for p, _ in PAR for i, f in zip(("lo", "hi"), FAC(p))}).T
rel.index = [f"{p}_{i}" for p, i in rel.index]; rel.to_csv(f"{R}/ace_sens_rel.csv", float_format="%.3f")
ok = {p: bool((rel.loc[f"{p}_lo"].abs() <= 5).all() and (rel.loc[f"{p}_hi"].abs() <= 5).all()) for p, _ in PAR}
print("parameters within 5% on all six observables (both factors):", [p for p, v in ok.items() if v], "-> runs", 2 * sum(ok.values()))
COL = dict(glc=SC["glucose"], ace=SC["acetate"], ref="#9a9a9a", fg="#222222", ok="#2e8b57")
def get(run): d = S[S.run == run].sort_values("t"); return d.t.values / 60, d.glucose_mean.values, d.acetate_mean.values
fig = plt.figure(figsize=(6.30, 9.74)); gs = GridSpec(2, 1, height_ratios=[4.1, 3.3], hspace=0.27, figure=fig)
ga = gs[0].subgridspec(4, 4, hspace=0.95, wspace=0.42); t0, g0, a0 = get("ACE-REF_s50"); ymax = 17.0
for i, (p, lab) in enumerate(PAR):
    ax = fig.add_subplot(ga[i // 4, i % 4]); ax.plot(t0, g0, color=COL["ref"], lw=2.58, alpha=0.55, solid_capstyle="round"); ax.plot(t0, a0, color=COL["ref"], lw=2.58, alpha=0.55, solid_capstyle="round")
    for f, ls, lw in ((FAC(p)[0], (0, (4, 2)), 1.5), (FAC(p)[1], "-", 1.5)):
        t, g, a = get(rid(p, f)); ax.plot(t, g, color=COL["glc"], ls=ls, lw=lw); ax.plot(t, a, color=COL["ace"], ls=ls, lw=lw)
    ax.set_xlim(0, 10); ax.set_ylim(0, ymax); ax.set_xticks([0, 5, 10]); ax.set_yticks([0, 5, 10, 15]); ax.set_title(lab, fontsize=FS["title"], color=COL["fg"], loc="left")
    ax.tick_params(labelsize=FS["label"]); [ax.spines[s].set_visible(False) for s in ("top", "right")]
    if ok[p]: [ax.spines[s].set_color(COL["ok"]) or ax.spines[s].set_linewidth(1.3) for s in ("left", "bottom")]
    if i % 4 == 0: ax.set_ylabel("Mean conc. (mM)", fontsize=FS["label"])
    if i // 4 == 3 or i >= 9: ax.set_xlabel("Time (h)", fontsize=FS["label"])
axl = fig.add_subplot(ga[3, 2:]); axl.axis("off")
from matplotlib.lines import Line2D
h = [Line2D([], [], color=COL["ref"], lw=2.58, alpha=0.55, label="reference (ACE-REF)"), Line2D([], [], color="k", lw=0.86, ls=(0, (4, 2)), label="parameter $\\times$0.5\n($\\times$0.44, volume pair)"), Line2D([], [], color="k", lw=0.86, label="parameter $\\times$2"),
     Line2D([], [], color=COL["glc"], lw=1.15, label="glucose"), Line2D([], [], color=COL["ace"], lw=1.15, label="acetate"), Line2D([], [], color=COL["ok"], lw=1.26, label="both runs within 5% on\nall six observables (B)")]
axl.legend(handles=h, loc="center left", ncol=2, frameon=False, fontsize=FS["label"], labelspacing=0.9, columnspacing=1.2)
fig.text(0.03, 0.965, "A", fontsize=FS["letter"], fontweight="bold")
axb = fig.add_subplot(gs[1]); M = np.full((len(PAR), 12), np.nan)
for i, (p, _) in enumerate(PAR):
    for j, (o, _) in enumerate(OBS): M[i, 2 * j] = rel.loc[f"{p}_lo", o]; M[i, 2 * j + 1] = rel.loc[f"{p}_hi", o]
cm = plt.get_cmap("RdBu_r").copy(); cm.set_bad("#d9d9d9"); axb.imshow(np.clip(M, -100, 100), cmap=cm, norm=TwoSlopeNorm(0, -100, 100), aspect="auto", alpha=1)
for i in range(M.shape[0]):
    for j in range(12):
        v = M[i, j]
        if np.isnan(v): axb.text(j, i, "none", ha="center", va="center", fontsize=FS["tick"], color="#444"); continue
        if abs(v) <= 5: axb.add_patch(plt.Rectangle((j - .5, i - .5), 1, 1, color="white", alpha=0.9, lw=0.00))
        axb.text(j, i, f"{v:+.0f}" if abs(v) >= 1 else ("0" if abs(v) < 0.05 else f"{v:+.1f}"), ha="center", va="center", fontsize=FS["tick"], color="#777" if abs(v) <= 5 else ("white" if abs(v) > 60 else "#111"), fontweight="normal" if abs(v) <= 5 else "bold")
axb.set_xticks(range(12)); axb.set_xticklabels(["0.5", "2"] * 6, fontsize=FS["label"]); axb.set_yticks(range(len(PAR))); axb.set_yticklabels([l for _, l in PAR], fontsize=FS["label"])
for j, (_, lab) in enumerate(OBS):
    axb.text(2 * j + .5, -1.15, lab, ha="center", va="bottom", fontsize=FS["label"])

    if j: axb.axvline(2 * j - .5, color="white", lw=1.72)
axb.set_xlabel("parameter factor (0.44 and 2 for the volume pair)", fontsize=FS["label"]); axb.tick_params(length=0); [s.set_visible(False) for s in axb.spines.values()]
axb.set_title("Change relative to the reference (%)", fontsize=FS["title"], loc="left", pad=46); fig.text(0.03, 0.425, "B", fontsize=FS["letter"], fontweight="bold")
fig.subplots_adjust(left=0.1, right=0.97, top=0.95, bottom=0.05)
for e in ("png", "pdf"): fig.savefig(f"{R}/fig_acetate_sens.{e}", dpi=200, bbox_inches="tight")
