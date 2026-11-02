"""Colony figure (manuscript Fig. 4 layout) from a quarter-colony run: A total biomass and cells vs time (x4 = full colony), B fields at TIME, C-E single-cell maps in the 80-180 um window.
usage: make_figure.py RUNDIR OUT_PREFIX [TIME_MIN=540]      (the crop window lies inside the quadrant x >= 0, y >= 0, so no mirroring is needed)
Conventions of the original notebook (bacterial_colony_figures.ipynb): dry weight = volume x 0.75 x 1.04 pg/um3; radius from the volume as a sphere; brown = dead, grey = living cell with a necrosis rate > 0;
glucose uptake and acetate exchange in fmol/min = flux (mmol/gDW/h) x volume x 0.75 x 1.04 / 60."""
import sys, glob, os, numpy as np, scipy.io as sio, xml.etree.ElementTree as ET, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.collections import PatchCollection
from matplotlib.patches import Circle
from matplotlib.colors import LinearSegmentedColormap, Normalize, TwoSlopeNorm
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..", "common")); from pcxml import labels
run, out = sys.argv[1].rstrip("/"), sys.argv[2]; T = float(sys.argv[3]) if len(sys.argv) > 3 else 540.0
XL = YL = (80, 180); FULL = 4; DRY = 0.75 * 1.04
def load(s):
    root = ET.parse(s).getroot(); lab = labels(root); m = sio.loadmat(s.replace(".xml", "_cells.mat"))["cells"]
    return float(root.find(".//current_time").text), {k: m[lab.index(k)] for k in lab if "[" not in k or k.startswith(("position", "death_rates"))}
ts, bio, cnt, snap = [], [], [], None
for s in sorted(glob.glob(f"{run}/output*[0-9].xml")):
    t, d = load(s); ts.append(t / 60); v = d["total_volume"]; bio.append(FULL * v.sum() * DRY / 1e6); cnt.append(FULL * len(v))
    if abs(t - T) < 1e-6: snap = (s, d)
ts, bio, cnt = map(np.array, (ts, bio, cnt)); s, d = snap
# segmented regression of ln(biomass) (as in the notebook): breakpoint with the smallest residual, >= 3 points per segment
y = np.log(bio[bio > 0]); t = ts[bio > 0]; best = None
for k in range(3, len(t) - 3):
    r = 0
    for a, b in ((slice(0, k + 1), None), (slice(k, None), None)):
        c = np.polyfit(t[a], y[a], 1); r += ((y[a] - np.polyval(c, t[a])) ** 2).sum()
    if best is None or r < best[0]: best = (r, k)
Tb = t[best[1]]; print("transition at %.1f h; biomass %.4f ug, cells %d at the end" % (Tb, bio[-1], cnt[-1]))
plt.rcParams.update({"font.size": 11}); glc = LinearSegmentedColormap.from_list("y", ["#ffffff", "#f1c232"])
fig = plt.figure(figsize=(13, 9.6), dpi=200)
top = fig.add_gridspec(1, 2, width_ratios=[1, 1.9], left=0.07, right=0.97, top=0.95, bottom=0.55, wspace=0.38)
bot = fig.add_gridspec(1, 3, left=0.05, right=0.97, top=0.45, bottom=0.07, wspace=0.42)
ax = fig.add_subplot(top[0, 0]); ax.text(-0.3, 1.04, "A", transform=ax.transAxes, fontsize=22, fontweight="bold"); ax.scatter(ts, bio, c="g", edgecolors="k", s=28, label="Biomass"); ax.axvline(Tb, color="k", ls="--", lw=1)
ax.set_xlabel("Time (h)"); ax.set_ylabel(r"Total biomass ($\mu$g)"); ax2 = ax.twinx(); ax2.plot(ts, cnt, "^-", color="grey", mec="k", ms=5, lw=1.5, label="Cells")
ax2.set_ylabel("Total cells"); ax2.ticklabel_format(style="sci", axis="y", scilimits=(0, 0)); h = ax.get_legend_handles_labels()[0] + ax2.get_legend_handles_labels()[0]
ax.legend(h, ["Biomass", "Cells"], frameon=False, loc="upper left")
for a in (ax, ax2): a.spines["top"].set_visible(False)
me = sio.loadmat(s.replace(".xml", "_microenvironment0.mat"))["multiscale_microenvironment"]
names = [("oxygen", 4, "Blues", r"$O_2$"), ("glucose", 5, glc, "D-Glucose"), ("acetate", 6, "Reds", "Acetate"), ("CO2", 8, "Purples", r"$CO_2$")]
sub = top[0, 1].subgridspec(2, 2, hspace=0.4, wspace=0.55)
xs, ys = np.unique(me[0]), np.unique(me[1])
for i, (nm, row, cm, lab) in enumerate(names):
    a = fig.add_subplot(sub[i // 2, i % 2]); g = np.full((len(ys), len(xs)), np.nan); g[np.searchsorted(ys, me[1]), np.searchsorted(xs, me[0])] = me[row]
    im = a.imshow(g, origin="lower", extent=[xs[0] - 1, xs[-1] + 1, ys[0] - 1, ys[-1] + 1], cmap=cm, vmin=0, vmax=me[row].max(), interpolation="bicubic", aspect="auto")
    a.set_xlim(*XL); a.set_ylim(*YL); a.set_title(lab, fontsize=12, pad=3); a.tick_params(labelsize=9); plt.colorbar(im, ax=a, label="mM", pad=0.02)
    if i == 0: a.text(-0.42, 1.12, "B", transform=a.transAxes, fontsize=22, fontweight="bold")
    if i // 2 == 1: a.set_xlabel("x (µm)")
    if i % 2 == 0: a.set_ylabel("y (µm)")
x, yv, vol, dead, nec = d["position[0]"], d["position[1]"], d["total_volume"], d["dead"] > 0, d["death_rates[1]"] > 0
sel = (x >= XL[0]) & (x <= XL[1]) & (yv >= YL[0]) & (yv <= YL[1]); rad = np.cbrt(3 * vol / (4 * np.pi)); alive = ~dead & ~nec
vals = {"gr": (d["growth_rate"], plt.cm.Greens, "Growth rate (1/h)"),
        "glc": (-d["R_EX_glc__D_e"] * vol * DRY / 60, glc, "Glucose uptake (fmol/min)"),
        "ac": (d["R_EX_ac_e"] * vol * DRY / 60, plt.cm.coolwarm_r, "Acetate exchange (fmol/min)")}
for j, (k, (v, cm, lab)) in enumerate(vals.items()):
    a = fig.add_subplot(bot[0, j]); a.text(-0.08, 1.04, "CDE"[j], transform=a.transAxes, fontsize=22, fontweight="bold"); vv = v[sel & alive]; lo, hi = vv.min(), vv.max()
    norm = TwoSlopeNorm(vmin=lo, vcenter=0.0, vmax=hi) if lo < 0 < hi else Normalize(lo, hi)
    idx = np.where(sel)[0]; patches = [Circle((x[i], yv[i]), rad[i]) for i in idx]
    col = ["saddlebrown" if dead[i] else ("grey" if nec[i] else cm(norm(v[i]))) for i in idx]
    a.add_collection(PatchCollection(patches, facecolors=col, edgecolors="k", linewidths=0.3))
    a.set_xlim(*XL); a.set_ylim(*YL); a.set_aspect("equal"); a.set_xticks([]); a.set_yticks([])
    a.plot([XL[1] - 55, XL[1] - 5], [YL[0] - 5.5] * 2, "k-", lw=2.5, clip_on=False); a.text(XL[1] - 30, YL[0] - 9, "50 µm", ha="center", va="top", fontsize=10)
    sm = plt.cm.ScalarMappable(norm=norm, cmap=cm); plt.colorbar(sm, ax=a, label=lab, fraction=0.045, pad=0.03)
fig.savefig(out + ".png", bbox_inches="tight", facecolor="white"); fig.savefig(out + ".pdf", bbox_inches="tight", facecolor="white"); print("wrote", out)
