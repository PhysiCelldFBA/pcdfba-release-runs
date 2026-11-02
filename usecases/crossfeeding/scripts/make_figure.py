"""Cross-feeding figure (manuscript Fig. 6 layout) from one run: A schematic (cropped from the current manuscript figure), B biomass volume per species, C box-mean concentrations,
D net export of each species (fmol/s), E growth rate and H2 / acetate fields with streamlines of the concentration gradient (+grad c, as in the manuscript: arrows point to higher concentration) at an early and a late time.
usage: make_figure.py RUNDIR OUT_PREFIX SCHEMATIC.png [T_EARLY_MIN=30] [T_LATE_MIN=last snapshot]
Net export = sum over living cells of flux (mmol/gDW/h) x dry weight (volume x (1 - fluid fraction) x 1.04e-12 g/um3) x 1e12 fmol/mmol / 3600 s/h; positive = export."""
import sys, glob, os, numpy as np, scipy.io as sio, xml.etree.ElementTree as ET, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt, matplotlib.image as mpimg
from matplotlib.colors import LinearSegmentedColormap
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..", "common")); from pcxml import labels
run, out, schem = sys.argv[1].rstrip("/"), sys.argv[2], sys.argv[3]
TE = float(sys.argv[4]) if len(sys.argv) > 4 else 30.0
files = sorted(glob.glob(f"{run}/output*[0-9].xml")); VAR = {"glucose": 4, "methane": 5, "CO2": 6, "H2": 7, "acetate": 11}
EXP = {"Glucose": "R_EX_glc_D_e", "H$_2$": "R_EX_h2_e", "CO$_2$": "R_EX_co2_e", "Acetate": "R_EX_ac_e", "CH$_4$": "R_EX_ch4_e"}
ts, vol, conc, exp, snaps = [], [], [], [], {}
for f in files:
    root = ET.parse(f).getroot(); lab = labels(root); t = float(root.find(".//current_time").text)
    c = sio.loadmat(f.replace(".xml", "_cells.mat"))["cells"]; g = lambda k: c[lab.index(k)]
    ty, v, ff, dead = g("cell_type"), g("total_volume"), g("fluid_fraction"), g("dead") > 0; dw = v * (1 - ff) * 1.04e-12; alive = ~dead
    me = sio.loadmat(f.replace(".xml", "_microenvironment0.mat"))["multiscale_microenvironment"]
    ts.append(t / 60); vol.append([v[(ty == k) & alive].sum() for k in (0, 1)]); conc.append({k: me[r].mean() for k, r in VAR.items()})
    exp.append({(k, nm): (g(col) * dw * 1e12 / 3600)[(ty == k) & alive].sum() for k in (0, 1) for nm, col in EXP.items()})
    snaps[t] = (me, g("position[0]"), g("position[1]"), g("growth_rate"), ty, dead)
ts, vol = np.array(ts), np.array(vol); TL = float(sys.argv[5]) if len(sys.argv) > 5 else max(snaps); TE = min(snaps, key=lambda t: abs(t - TE)); TL = min(snaps, key=lambda t: abs(t - TL))
print("cells CB/MB at the end:", [int(((snaps[max(snaps)][4] == k) & ~snaps[max(snaps)][5]).sum()) for k in (0, 1)], "t_end %.1f h" % ts[-1], "early %g late %g min" % (TE, TL))
plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False})
fig = plt.figure(figsize=(14.5, 8.4), dpi=200); og = fig.add_gridspec(2, 3, height_ratios=[0.95, 1.8], width_ratios=[1, 1, 1.15], left=0.04, right=0.98, top=0.97, bottom=0.07, wspace=0.28, hspace=0.18)
a = fig.add_subplot(og[0, 0:2]); im = mpimg.imread(schem); H, W = im.shape[:2]; a.imshow(im[: int(0.375 * H), : int(0.56 * W)]); a.axis("off")
lb = og[1, 0].subgridspec(2, 1, hspace=0.5); b = fig.add_subplot(lb[0]); c_ = fig.add_subplot(lb[1])
b.plot(ts, vol[:, 0], color="#1f77b4", label="C. beijerinckii"); b.plot(ts, vol[:, 1], color="#c0551a", label="M. barkeri"); b.set_ylabel("Biomass volume (µm³)"); b.set_xlabel("Time (h)"); b.legend(frameon=False, loc="upper left", fontsize=8)
cc = {"glucose": "#f0c030", "methane": "#18a081", "H2": "#3a6ea5", "acetate": "#c00000"}
for k, col in cc.items(): c_.plot(ts, [d[k] for d in conc], color=col, label={"glucose": "Glucose", "methane": "CH$_4$", "H2": "H$_2$", "acetate": "Acetate"}[k])
c_.set_ylabel("Concentration (mM)"); c_.set_xlabel("Time (h)"); c_.legend(frameon=False, ncol=4, fontsize=8, loc="upper center", bbox_to_anchor=(0.5, -0.38))
for q, lab_ in ((b, "B"), (c_, "C")): q.text(-0.2, 1.08, lab_, transform=q.transAxes, fontsize=20, fontweight="bold")
dg = og[1, 1].subgridspec(2, 1, hspace=0.45); col = {"Glucose": "#f0c030", "H$_2$": "#3a6ea5", "CO$_2$": "#8e44ad", "Acetate": "#c00000", "CH$_4$": "#18a081"}
for i, (k, names, nm_) in enumerate(((0, ["Glucose", "H$_2$", "CO$_2$", "Acetate"], "C. beijerinckii"), (1, ["H$_2$", "Acetate", "CO$_2$", "CH$_4$"], "M. barkeri"))):
    d = fig.add_subplot(dg[i])
    for nm in names: d.plot(ts, [e[(k, nm)] for e in exp], color=col[nm], label=nm)
    d.axhline(0, color="grey", ls="--", lw=0.8); d.set_ylabel("Net export rate (fmol/s)\n" + nm_); d.legend(frameon=False, fontsize=8, loc="center left", bbox_to_anchor=(1.0, 0.5))
    if i == 1: d.set_xlabel("Time (h)")
    if i == 0: d.text(-0.25, 1.08, "D", transform=d.transAxes, fontsize=20, fontweight="bold")
eg = og[:, 2].subgridspec(3, 2, hspace=0.12, wspace=0.08); cb_ = plt.cm.Blues; mb_ = LinearSegmentedColormap.from_list("o", ["#fdf1e6", "#8b2d00"])
for j, T in enumerate((TE, TL)):
    me, x, y, gr, ty, dead = snaps[T]; xs, ys = np.unique(me[0]), np.unique(me[1]); ext = [xs[0] - (xs[1] - xs[0]) / 2, xs[-1] + (xs[1] - xs[0]) / 2] * 2
    ax = fig.add_subplot(eg[0, j]); ax.set_title("t = %.1f h" % (T / 60), fontsize=11)
    for k, cm in ((0, cb_), (1, mb_)):
        s = (ty == k) & ~dead; ax.scatter(x[s], y[s], c=gr[s], cmap=cm, vmin=0, vmax=max(0.12 if k == 0 else 0.02, 1e-9), s=2.5, lw=0)
    lim = ext[1]; ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim); ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values(): sp.set_visible(True)
    if j == 0: ax.set_ylabel("Growth rate", fontsize=11); ax.text(-0.2, 1.1, "E", transform=ax.transAxes, fontsize=20, fontweight="bold")
    for i, (nm, row, cm, lab_) in enumerate((("H2", 7, "Blues", "H$_2$ (mM)"), ("acetate", 11, "Reds", "Acetate (mM)"))):
        ax = fig.add_subplot(eg[1 + i, j]); g = np.zeros((len(ys), len(xs))); g[np.searchsorted(ys, me[1]), np.searchsorted(xs, me[0])] = me[row]
        ax.imshow(g, origin="lower", extent=ext, cmap=cm, vmin=0, vmax=max(g.max(), 1e-12)); gy, gx = np.gradient(g, ys, xs)
        ax.streamplot(xs, ys, gx, gy, color="#333333", density=0.9, linewidth=0.7, arrowsize=0.8); ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim); ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([])
        for sp in ax.spines.values(): sp.set_visible(True)
        if j == 0: ax.set_ylabel(lab_, fontsize=11)
for k, (cm, vm, lab_, y0) in enumerate(((cb_, 0.12, "CB (1/h)", 0.85), (mb_, 0.02, "MB (1/h)", 0.71))):
    cax = fig.add_axes([0.985, y0 - 0.02, 0.007, 0.07]); fig.colorbar(plt.cm.ScalarMappable(norm=plt.Normalize(0, vm), cmap=cm), cax=cax, label=lab_)
fig.savefig(out + ".png", bbox_inches="tight", facecolor="white"); fig.savefig(out + ".pdf", bbox_inches="tight", facecolor="white")
