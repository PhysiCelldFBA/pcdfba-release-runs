"""Cross-feeding figure (manuscript Fig. 6) from one run (open box, no caps, 50 h: CF3-NC50_s50): A biomass volume per species, B box-mean concentrations, D mean net export rate per living cell (amol/s) of C. beijerinckii and M. barkeri (shared y-axis),
B mean growth rate of the living cells of each species (the closed-box value of M. barkeri, CF3-NF_s50, is only printed), E growth rate and H2, acetate and CH4 fields at an early and a late time (colours of the substrates from common/plot_style.py; five contour lines at equal steps of the maximum, without the outer voxel layer (the boundary value)).
usage: make_figure.py RUNDIR OUT_PREFIX [T_EARLY_MIN=30] [T_LATE_MIN=last snapshot] [CLOSED_RUNDIR=runs/CF3-NF_s50]
Net export per cell = mean over living cells of flux (mmol/gDW/h) x dry weight (volume x (1 - fluid fraction) x 1.04e-12 g/um3) x 1e15 amol/mmol / 3600 s/h; positive = export.
The schematic of the submitted Fig. 6 is now a supplementary figure (make_schematic.py)."""
import sys, glob, os, numpy as np, scipy.io as sio, xml.etree.ElementTree as ET, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..", "common")); from pcxml import labels; from plot_style import SC, substrate_cmap
run, out = sys.argv[1].rstrip("/"), sys.argv[2]
TE = float(sys.argv[3]) if len(sys.argv) > 3 else 30.0
files = sorted(glob.glob(f"{run}/output*[0-9].xml")); VAR = {"glucose": 4, "methane": 5, "CO2": 6, "H2": 7, "acetate": 11}
EXP = {"Glucose": "R_EX_glc_D_e", "H$_2$": "R_EX_h2_e", "CO$_2$": "R_EX_co2_e", "Acetate": "R_EX_ac_e", "CH$_4$": "R_EX_ch4_e"}
def tof(f): return float(ET.parse(f).getroot().find(".//current_time").text)
tl = [tof(files[0]), tof(files[-1])]; TL = float(sys.argv[4]) if len(sys.argv) > 4 else tl[1]
CLOSED = sys.argv[5] if len(sys.argv) > 5 else os.path.join(os.path.dirname(run), "CF3-NF_s50")
ts, vol, conc, exp, grm, snaps = [], [], [], [], [], {}
for f in files:
    root = ET.parse(f).getroot(); lab = labels(root); t = float(root.find(".//current_time").text)
    c = sio.loadmat(f.replace(".xml", "_cells.mat"))["cells"]; g = lambda k: c[lab.index(k)]
    ty, v, ff, dead = g("cell_type"), g("total_volume"), g("fluid_fraction"), g("dead") > 0; dw = v * (1 - ff) * 1.04e-12; alive = ~dead
    ts.append(t / 60); vol.append([v[(ty == k) & alive].sum() for k in (0, 1)]); grm.append([g("growth_rate")[(ty == k) & alive].mean() for k in (0, 1)])
    exp.append({(k, nm): (g(col) * dw * 1e15 / 3600)[(ty == k) & alive].mean() for k in (0, 1) for nm, col in EXP.items()})
    me = sio.loadmat(f.replace(".xml", "_microenvironment0.mat"))["multiscale_microenvironment"]; conc.append({k: me[r].mean() for k, r in VAR.items()})
    if abs(t - TE) < 1 or abs(t - TL) < 1: snaps[t] = (me, g("position[0]"), g("position[1]"), g("growth_rate"), ty, dead)
ts, vol, grm = np.array(ts), np.array(vol), np.array(grm); TE = min(snaps, key=lambda t: abs(t - TE)); TL = min(snaps, key=lambda t: abs(t - TL))
fc = sorted(glob.glob(f"{CLOSED}/output*[0-9].xml"))[-1]; rc = ET.parse(fc).getroot(); lc = labels(rc); cc_ = sio.loadmat(fc.replace(".xml", "_cells.mat"))["cells"]
closed_mb = float(cc_[lc.index("growth_rate")][(cc_[lc.index("cell_type")] == 1) & (cc_[lc.index("dead")] == 0)].mean())
print("cells CB/MB at the end:", [int(((snaps[TL][4] == k) & ~snaps[TL][5]).sum()) for k in (0, 1)], "t_end %.1f h" % ts[-1], "early %g late %g min" % (TE, TL), "closed-box MB growth %.4f" % closed_mb)
plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False})
fig = plt.figure(figsize=(15.5, 9.4), dpi=200); og = fig.add_gridspec(2, 1, height_ratios=[0.55, 1.5], left=0.06, right=0.985, top=0.965, bottom=0.04, hspace=0.20)
top = og[0].subgridspec(1, 4, wspace=0.34, width_ratios=[1, 1, 1, 2.1]); LET = lambda q, s, x=-0.2: q.text(x, 1.06, s, transform=q.transAxes, fontsize=20, fontweight="bold")
cb_c, mb_c = "#1f77b4", "#c0551a"
a = fig.add_subplot(top[0]); a.plot(ts, vol[:, 0], color=cb_c, label="C. beijerinckii"); a.plot(ts, vol[:, 1], color=mb_c, label="M. barkeri"); a.set_yscale("log"); a.set_ylabel("Biomass volume (µm³)"); a.set_xlabel("Time (h)"); a.legend(frameon=False, loc="upper left", fontsize=8.5); LET(a, "A")
e = fig.add_subplot(top[1]); e.plot(ts[1:], grm[1:, 0], color=cb_c, label="C. beijerinckii"); e.plot(ts[1:], grm[1:, 1], color=mb_c, label="M. barkeri")
e.set_ylabel("Mean growth rate (1/h)"); e.set_xlabel("Time (h)"); e.set_ylim(0, None); e.legend(frameon=False, fontsize=8.5, loc="center left", bbox_to_anchor=(0.02, 0.61)); LET(e, "B", -0.24)
b = fig.add_subplot(top[2]); cc = {"glucose": SC["glucose"], "methane": SC["methane"], "H2": SC["h2"], "acetate": SC["acetate"]}
for k, col in cc.items(): b.plot(ts[ts >= 1 / 60], np.array([d[k] for d in conc])[ts >= 1 / 60], color=col, label={"glucose": "Glucose", "methane": "CH$_4$", "H2": "H$_2$", "acetate": "Acetate"}[k])
b.set_ylabel("Concentration (mM)"); b.set_xlabel("Time (h)"); b.legend(frameon=False, ncol=1, fontsize=8.5, loc="center left", bbox_to_anchor=(0.02, 0.40)); LET(b, "C")
col = {"Glucose": SC["glucose"], "H$_2$": SC["h2"], "CO$_2$": SC["co2"], "Acetate": SC["acetate"], "CH$_4$": SC["methane"]}
xp = top[3].subgridspec(1, 2, wspace=0.06); dax = None
for i, (k, names, nm_) in enumerate(((0, ["Glucose", "H$_2$", "CO$_2$", "Acetate"], "C. beijerinckii"), (1, ["H$_2$", "Acetate", "CO$_2$", "CH$_4$"], "M. barkeri"))):
    d = fig.add_subplot(xp[i], sharey=dax); dax = dax or d
    for nm in names: d.plot(ts, [e_[(k, nm)] for e_ in exp], color=col[nm], label=nm)
    d.axhline(0, color="grey", ls="--", lw=0.8); d.set_xlabel("Time (h)"); d.set_title(nm_, fontsize=10, style="italic"); d.legend(frameon=False, fontsize=8.5, loc="upper left", ncol=2)
    d.set_ylim(-3.4, 3.9)
    if i == 0: d.set_ylabel("Net export rate per cell (amol/s)"); LET(d, "D", -0.24)
    else: d.tick_params(labelleft=False)
eg = og[1].subgridspec(2, 4, hspace=0.10, wspace=0.42); cb_ = plt.cm.Blues; mb_ = LinearSegmentedColormap.from_list("o", ["#fdf1e6", "#8b2d00"])
for j, T in enumerate((TE, TL)):
    me, x, y, gr, ty, dead = snaps[T]; xs, ys = np.unique(me[0]), np.unique(me[1]); ext = [xs[0] - (xs[1] - xs[0]) / 2, xs[-1] + (xs[1] - xs[0]) / 2] * 2; lim = ext[1]
    ax = fig.add_subplot(eg[j, 0]); ax.set_ylabel("t = %.1f h" % (T / 60), fontsize=12)
    for k, cm, vm in ((0, cb_, 0.12), (1, mb_, 0.02)):
        s = (ty == k) & ~dead; ax.scatter(x[s], y[s], c=gr[s], cmap=cm, vmin=0, vmax=vm, s=2.5, lw=0)
    ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim); ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values(): sp.set_visible(True)
    for q, (cm, vm, lab_, y0) in enumerate(((cb_, 0.12, "CB (1/h)", 0.53), (mb_, 0.02, "MB (1/h)", 0.0))):
        cax = ax.inset_axes([1.04, y0, 0.05, 0.47]); cbx = fig.colorbar(plt.cm.ScalarMappable(norm=plt.Normalize(0, vm), cmap=cm), cax=cax); cbx.set_label(lab_, fontsize=8, labelpad=2); cbx.ax.tick_params(labelsize=7)
    if j == 0: ax.text(-0.2, 1.12, "E", transform=ax.transAxes, fontsize=20, fontweight="bold"); ax.set_title("Growth rate", fontsize=11, pad=6)
    for i, (nm, row, cm, lab_) in enumerate((("H2", 7, substrate_cmap("h2"), "H$_2$ (mM)"), ("acetate", 11, substrate_cmap("acetate"), "Acetate (mM)"), ("methane", 5, substrate_cmap("methane"), "CH$_4$ (mM)"))):
        ax = fig.add_subplot(eg[j, 1 + i]); g = np.zeros((len(ys), len(xs))); g[np.searchsorted(ys, me[1]), np.searchsorted(xs, me[0])] = me[row]
        imh = ax.imshow(g, origin="lower", extent=ext, cmap=cm, vmin=0, vmax=max(g.max(), 1e-12))
        cbx = fig.colorbar(imh, cax=ax.inset_axes([1.04, 0, 0.05, 1])); cbx.ax.tick_params(labelsize=7); cbx.formatter.set_powerlimits((-2, 3)); cbx.update_ticks(); cbx.ax.yaxis.get_offset_text().set_fontsize(7)
        ax.contour(xs[1:-1], ys[1:-1], g[1:-1, 1:-1], levels=np.linspace(0, g.max(), 7)[1:-1], colors="k", linewidths=0.6, alpha=0.55)
        ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim); ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([])
        for sp in ax.spines.values(): sp.set_visible(True)
        if j == 0: ax.set_title(lab_, fontsize=11, pad=6)
fig.savefig(out + ".png", bbox_inches="tight", pad_inches=0.2, facecolor="white"); fig.savefig(out + ".pdf", bbox_inches="tight", pad_inches=0.2, facecolor="white")
