"""Colony figure (manuscript Fig. 4) from a quarter-colony run: A total biomass and cells against time (x4 = full colony) with the end of exponential growth; B fields of O2, glucose, acetate and CO2 over the quadrant at TIME with the boundaries of the zones of the colony; C-F single-cell maps (growth rate, glucose uptake, acetate exchange, oxygen uptake) in the 80-180 um window, with the same boundaries.
usage: make_figure.py RUNDIR OUT_PREFIX [TIME_MIN=540]      (run COL-QTRB_qtr_s50; the quadrant x >= 0, y >= 0 is simulated, so the zones are arcs centred on the origin)
Conventions of the original notebook (bacterial_colony_figures.ipynb): dry weight = volume x 0.75 x 1.04 pg/um3; radius from the volume as a sphere; brown = dead, grey = living cell with a necrosis rate > 0; uptake and exchange in amol/s = flux (mmol/gDW/h) x volume x 0.75 x 1.04 / 3.6.
Zones (10 um annuli of the distance to the centre, living cells for the fractions of growth and acetate): necrotic core = out to the edge of the inner region where at least 50% of the cells are dead; stalled zone = beyond it, until at least 50% of the living cells grow (rate > 0.05 1/h);
acetate-using ring = beyond it, until fewer than 50% of the living cells take up acetate; proliferating rim = beyond it, out to the edge of the colony (99th percentile of the distance)."""
import matplotlib.ticker
import sys, glob, os, numpy as np, scipy.io as sio, xml.etree.ElementTree as ET, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.collections import PatchCollection
from matplotlib.patches import Circle, Arc, Rectangle, Patch
from matplotlib.lines import Line2D
from matplotlib.colors import LinearSegmentedColormap, Normalize, TwoSlopeNorm
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..", "common")); from pcxml import labels; from plot_style import SC, substrate_cmap, apply_style, FS
run, out = sys.argv[1].rstrip("/"), sys.argv[2]; T = float(sys.argv[3]) if len(sys.argv) > 3 else 540.0
XL = YL = (80, 180); FULL = 4; DRY = 0.75 * 1.04; VIEW = 210.0; AN = 10.0
def load(s):
    root = ET.parse(s).getroot(); lab = labels(root); m = sio.loadmat(s.replace(".xml", "_cells.mat"))["cells"]
    return float(root.find(".//current_time").text), {k: m[lab.index(k)] for k in lab if "[" not in k or k.startswith(("position", "death_rates"))}
ts, bio, cnt, snap = [], [], [], None
for s in sorted(glob.glob(f"{run}/output*[0-9].xml")):
    t, d = load(s); ts.append(t / 60); v = d["total_volume"]; bio.append(FULL * v.sum() * DRY / 1e6); cnt.append(FULL * len(v))
    if abs(t - T) < 1e-6: snap = (s, d)
ts, bio, cnt = map(np.array, (ts, bio, cnt)); s, d = snap
y = np.log(bio[bio > 0]); t = ts[bio > 0]; best = None          # segmented regression of ln(biomass), breakpoint with the smallest residual, >= 3 points per segment
for k in range(3, len(t) - 3):
    r = 0
    for a in (slice(0, k + 1), slice(k, None)):
        c = np.polyfit(t[a], y[a], 1); r += ((y[a] - np.polyval(c, t[a])) ** 2).sum()
    if best is None or r < best[0]: best = (r, k)
Tb = t[best[1]]
x, yv, vol, dead, nec = d["position[0]"], d["position[1]"], d["total_volume"], d["dead"] > 0, d["death_rates[1]"] > 0
rr = np.hypot(x, yv); alive = ~dead & ~nec; gr, ac = d["growth_rate"], d["R_EX_ac_e"]
edges = np.arange(0, rr.max() + AN, AN); ia = np.digitize(rr, edges) - 1; nb = len(edges) - 1
def frac(cond, base):
    return np.array([(cond & base & (ia == i)).sum() / max((base & (ia == i)).sum(), 1) for i in range(nb)])
fdead = frac(dead, np.ones_like(dead)); live = ~dead; fgrow = frac(gr > 0.05, live); fac = frac(ac < -1e-6, live)
i_core = next(i for i in range(nb) if fdead[i] < 0.5); R_core = edges[i_core]                     # end of the inner region with >= 50% dead cells
i_stall = next(i for i in range(i_core, nb) if fgrow[i] >= 0.5); R_stall = edges[i_stall]         # growth resumes
i_ring = next(i for i in range(i_stall, nb) if fac[i] < 0.5); R_ring = edges[i_ring]              # end of the acetate-using ring
R_edge = float(np.percentile(rr, 99)); ZONES = [("Necrotic core", R_core), ("Stalled zone", R_stall), ("Acetate-using ring", R_ring)]
print("transition at %.1f h; biomass %.4f ug, cells %d at the end; zones: core %g, stalled %g, ring %g, edge %.0f um" % (Tb, bio[-1], cnt[-1], R_core, R_stall, R_ring, R_edge))
apply_style(6.3); glc = substrate_cmap("glucose")
fig = plt.figure(figsize=(6.3, 4.9), dpi=200)
gs = fig.add_gridspec(5, 6, height_ratios=[0.30, 1.25, 0.30, 1.5, 0.26], width_ratios=[1.5, 0.85, 1.14, 1.14, 1.14, 1.14], left=0.125, right=0.985, top=0.965, bottom=0.01, wspace=0.28, hspace=0.22)
LS = [("-", "k"), ("--", "k"), (":", "k")]
def zone_arcs(ax, lw=0.67):
    for (nm, R), (ls, c) in zip(ZONES, LS): ax.add_patch(Arc((0, 0), 2 * R, 2 * R, theta1=0, theta2=90, ls=ls, color=c, lw=lw, zorder=5))
ax = fig.add_subplot(gs[1, 0]); ax.text(-0.72, 1.04, "A", transform=ax.transAxes, fontsize=FS["letter"], fontweight="bold")
hb = ax.scatter(ts, bio, c="g", edgecolors="k", s=5, label="Biomass"); hv = ax.axvline(Tb, color="k", ls="--", lw=0.48, label="End of exponential growth (%.1f h)" % Tb)
ax.set_xlabel("Time (h)"); ax.set_ylabel(r"Total biomass ($\mu$g)"); ax2 = ax.twinx(); hc, = ax2.plot(ts, cnt / 1e5, "^-", color="grey", mec="k", ms=2.39, lw=0.72, label="Cells"); ax2.set_ylabel(r"Total cells ($\times 10^5$)"); ax2.spines["top"].set_visible(False)
ax.legend([hb, hc], ["Biomass", "Cells"], frameon=False, loc="upper left", fontsize=FS["tick"]); ax.text(Tb, 1.02, "End of exponential\nphase", transform=ax.get_xaxis_transform(), ha="center", va="bottom", fontsize=FS["label"])
me = sio.loadmat(s.replace(".xml", "_microenvironment0.mat"))["multiscale_microenvironment"]
names = [("oxygen", 4, "Blues", r"O$_2$ (mM)"), ("glucose", 5, glc, "D-glucose (mM)"), ("acetate", 6, "Reds", "Acetate (mM)"), ("CO2", 8, "Purples", r"CO$_2$ (mM)")]
maps = []; xs, ys = np.unique(me[0]), np.unique(me[1]); vis = (me[0] <= VIEW) & (me[1] <= VIEW);
for i, (nm, row, cm, lab) in enumerate(names):
    a = fig.add_subplot(gs[1, i + 2]); g = np.full((len(ys), len(xs)), np.nan); g[np.searchsorted(ys, me[1]), np.searchsorted(xs, me[0])] = me[row]
    im = a.imshow(g, origin="lower", extent=[xs[0] - 1, xs[-1] + 1, ys[0] - 1, ys[-1] + 1], cmap=cm, vmin=0, vmax=me[row][vis].max(), interpolation="bicubic", aspect="equal")
    a.set_xlim(0, VIEW); a.set_ylim(0, VIEW); a.set_anchor("S"); zone_arcs(a); a.add_patch(Rectangle((XL[0], YL[0]), XL[1] - XL[0], YL[1] - YL[0], fill=False, ec="#444444", lw=0.48, ls=(0, (1, 2)), zorder=6))
    a.set_xlabel("x (µm)"); a.set_xticks([0, 100, 200]); a.set_yticks([0, 100, 200]); 
    if i == 0: a.set_ylabel("y (µm)"); a.text(-0.42, 1.04, "B", transform=a.transAxes, fontsize=FS["letter"], fontweight="bold")
    else: a.set_yticklabels([])
    maps.append((a, im, lab))
fig.canvas.draw()
for a, im, lab in maps:      # colorbar directly above its map, same width
    a.apply_aspect(); pos = a.get_position(); cax = fig.add_axes([pos.x0, pos.y1 + 0.008, pos.width, 0.014]); cb = fig.colorbar(im, cax=cax, orientation="horizontal"); cb.set_label(lab, labelpad=7, fontsize=FS["label"])
    cb.ax.xaxis.set_ticks_position("top"); cb.ax.xaxis.set_label_position("top"); cb.ax.tick_params(labelsize=FS["tick"], length=2); cb.formatter.set_powerlimits((-9, 9)); cb.locator = matplotlib.ticker.MaxNLocator(2); cb.update_ticks()
pm = maps[0][0].get_position(); pa = ax.get_position(); ax.set_position([pa.x0, pm.y0, pa.width, pm.height]); ax2.set_position(ax.get_position())     # A as high as the maps of B
sel = (x >= XL[0]) & (x <= XL[1]) & (yv >= YL[0]) & (yv <= YL[1]); rad = np.cbrt(3 * vol / (4 * np.pi))
vals = [("Growth rate (1/h)", gr, plt.cm.Greens), ("Glucose uptake (amol/s)", -d["R_EX_glc__D_e"] * vol * DRY / 3.6, glc), ("Acetate exchange (amol/s)", d["R_EX_ac_e"] * vol * DRY / 3.6, plt.cm.coolwarm_r), ("Oxygen uptake (amol/s)", -d["R_EX_o2_e"] * vol * DRY / 3.6, plt.cm.Blues)]
legs = fig.add_subplot(gs[2, :]); legs.axis("off")
h = [Line2D([], [], color="k", ls=ls, lw=0.76) for ls, c in LS] + [Line2D([], [], color="#444444", ls=(0, (1, 2)), lw=0.57), Patch(fc="saddlebrown", ec="k", lw=0.19), Patch(fc="grey", ec="k", lw=0.19)]
lab_ = ["Core edge, %g µm" % R_core, "Stalled edge, %g µm" % R_stall, "Ring edge, %g µm" % R_ring, "Window of C–F", "Dead cell", "Starving cell"]
legs.legend(h, lab_, loc="upper center", bbox_to_anchor=(0.5, 0.40), ncol=6, frameon=True, fancybox=False, edgecolor="k", framealpha=1, borderpad=0.4, fontsize=FS["legend"], handlelength=1.5, columnspacing=0.9, labelspacing=0.3, handletextpad=0.4).get_frame().set_linewidth(0.6)
bb = gs[3, :].get_position(fig); posn = fig.add_gridspec(1, 4, left=0.03, right=0.985, bottom=bb.y0, top=bb.y1, wspace=0.16)
for j, (lab, v, cm) in enumerate(vals):
    a = fig.add_subplot(posn[0, j]); a.text(-0.1, 1.03, "CDEF"[j], transform=a.transAxes, fontsize=FS["letter"], fontweight="bold"); vv = v[sel & alive]; lo, hi = vv.min(), vv.max()
    norm = TwoSlopeNorm(vmin=lo, vcenter=0.0, vmax=hi) if lo < 0 < hi else Normalize(lo, hi)
    idx = np.where(sel)[0]; patches = [Circle((x[i], yv[i]), rad[i]) for i in idx]; col = ["saddlebrown" if dead[i] else ("grey" if nec[i] else cm(norm(v[i]))) for i in idx]
    a.add_collection(PatchCollection(patches, facecolors=col, edgecolors="k", linewidths=0.14)); zone_arcs(a, lw=0.86)
    a.set_xlim(*XL); a.set_ylim(*YL); a.set_aspect("equal"); a.set_xticks([]); a.set_yticks([]); [sp.set_visible(True) for sp in a.spines.values()]
    a.plot([XL[1] - 55, XL[1] - 5], [YL[1] - 8] * 2, "k-", lw=1.19); a.text(XL[1] - 30, YL[1] - 11, "50 µm", ha="center", va="top", fontsize=FS["label"])
    cax = a.inset_axes([0.1, -0.10, 0.8, 0.04]); cb = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cm), cax=cax, orientation="horizontal"); cb.set_label(lab, labelpad=2, fontsize=FS["label"]); cb.ax.tick_params(labelsize=FS["tick"], length=2)
fig.savefig(out + ".png", bbox_inches="tight", facecolor="white"); fig.savefig(out + ".pdf", bbox_inches="tight", facecolor="white"); print("wrote", out)
