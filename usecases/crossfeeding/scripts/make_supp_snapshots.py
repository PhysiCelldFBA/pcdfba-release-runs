"""Supplementary Figures S10 (growth rate of the cells) and S11 (glucose, CH4, H2 and acetate fields) of the cross-feeding use case from one run (CF3-NC50_s50: open box, no caps, 50 h)
at 0.5, 12.5, 25, 37.5 and 50 h. Colours from common/plot_style.py; the colour scale of each field is the same at all times (zero to the maximum over the five times).
usage: make_supp_snapshots.py RUNDIR OUT_PREFIX    ->  OUT_PREFIX_growth_rate.png, OUT_PREFIX_fields.png"""
import sys, glob, os, numpy as np, scipy.io as sio, xml.etree.ElementTree as ET, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..", "common")); from pcxml import labels; from plot_style import substrate_cmap, apply_style, FS
run, out = sys.argv[1].rstrip("/"), sys.argv[2]
TIMES = (30, 750, 1500, 2250, 3000); VAR = {"glucose": 4, "methane": 5, "H2": 7, "acetate": 11}
files = sorted(glob.glob(f"{run}/output*[0-9].xml")); tof = lambda f: float(ET.parse(f).getroot().find(".//current_time").text); tt = np.array([tof(f) for f in files])
snaps = {}
for T in TIMES:
    f = files[int(np.argmin(abs(tt - T)))]; root = ET.parse(f).getroot(); lab = labels(root)
    c = sio.loadmat(f.replace(".xml", "_cells.mat"))["cells"]; g = lambda k: c[lab.index(k)]
    me = sio.loadmat(f.replace(".xml", "_microenvironment0.mat"))["multiscale_microenvironment"]
    snaps[T] = (tof(f), me, g("position[0]"), g("position[1]"), g("growth_rate"), g("cell_type"), g("dead") > 0)
me0 = snaps[TIMES[0]][1]; xs, ys = np.unique(me0[0]), np.unique(me0[1]); h = xs[1] - xs[0]; ext = [xs[0] - h / 2, xs[-1] + h / 2] * 2; lim = ext[1]
apply_style(6.3); cb_ = plt.cm.Blues; mb_ = LinearSegmentedColormap.from_list("o", ["#fdf1e6", "#8b2d00"])
def frame(ax):
    ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim); ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([]); [sp.set_visible(True) for sp in ax.spines.values()]
# S10: growth rate
fig, axs = plt.subplots(1, 5, figsize=(6.3, 1.55), gridspec_kw=dict(wspace=0.04, right=0.87))
for ax, T in zip(axs, TIMES):
    t, me, x, y, gr, ty, dead = snaps[T]
    for k, cm, vm in ((0, cb_, 0.12), (1, mb_, 0.02)):
        s = (ty == k) & ~dead; ax.scatter(x[s], y[s], c=gr[s], cmap=cm, vmin=0, vmax=vm, s=0.5, lw=0)
    frame(ax); ax.set_title("t = %.1f h" % (t / 60))
for cm, vm, lab, y0 in ((cb_, 0.12, "CB (1/h)", 0.57), (mb_, 0.02, "MB (1/h)", 0.13)):
    cax = fig.add_axes([0.885, y0 * 0.85 + 0.03, 0.014, 0.31]); cbx = fig.colorbar(plt.cm.ScalarMappable(norm=plt.Normalize(0, vm), cmap=cm), cax=cax); cbx.set_label(lab)
fig.savefig(out + "_growth_rate.png", dpi=300, bbox_inches="tight", pad_inches=0.03, facecolor="white"); plt.close(fig)
# S11: fields
rows = (("glucose", "Glucose (mM)"), ("methane", "CH$_4$ (mM)"), ("H2", "H$_2$ (mM)"), ("acetate", "Acetate (mM)"))
fig, axs = plt.subplots(4, 5, figsize=(6.3, 5.4), gridspec_kw=dict(wspace=0.06, hspace=0.06, right=0.90))
for r, (nm, lab) in enumerate(rows):
    fields = []
    for T in TIMES:
        me = snaps[T][1]; g = np.zeros((len(ys), len(xs))); g[np.searchsorted(ys, me[1]), np.searchsorted(xs, me[0])] = me[VAR[nm]]; fields.append(g)
    vm = max(f.max() for f in fields); cm = substrate_cmap({"H2": "h2"}.get(nm, nm))
    for c_, (ax, g, T) in enumerate(zip(axs[r], fields, TIMES)):
        im = ax.imshow(g, origin="lower", extent=ext, cmap=cm, vmin=0, vmax=vm); frame(ax)
        if r == 0: ax.set_title("t = %.1f h" % (snaps[T][0] / 60))
        if c_ == 0: ax.set_ylabel(lab)
    pos = axs[r, -1].get_position(); cax = fig.add_axes([0.915, pos.y0 + 0.1 * pos.height, 0.014, 0.8 * pos.height]); cbx = fig.colorbar(im, cax=cax); cbx.formatter.set_powerlimits((-2, 3)); cbx.update_ticks()
fig.savefig(out + "_fields.png", dpi=300, bbox_inches="tight", pad_inches=0.03, facecolor="white"); plt.close(fig)
print({T: round(snaps[T][0] / 60, 2) for T in TIMES})
