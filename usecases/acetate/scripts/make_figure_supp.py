"""Supplementary figure of the acetate switch: glucose (top) and acetate (bottom) fields in the z = 0 slice of the closed box at three times, and the living cells (all z, projected); the min-max of the slice is printed in each panel (the box is well mixed).
usage: make_figure_supp.py RUNDIR OUT_PREFIX [T1=30 T2=120 T3=270 (min)]
Colour scales are shared by the row."""
import sys, os, numpy as np, scipy.io as sio, xml.etree.ElementTree as ET, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..", "common")); from pcxml import labels
run, out = sys.argv[1].rstrip("/"), sys.argv[2]; TS = [float(v) for v in sys.argv[3:6]] or [30.0, 120.0, 270.0]
if len(TS) != 3: TS = [30.0, 120.0, 270.0]
def snap(t):
    base = f"{run}/output{int(round(t / 5)) :08d}"; root = ET.parse(base + ".xml").getroot(); assert abs(float(root.find(".//current_time").text) - t) < 1e-6, base
    me = sio.loadmat(base + "_microenvironment0.mat")["multiscale_microenvironment"]; c = sio.loadmat(base + "_cells.mat")["cells"]; lab = labels(root)
    xs, ys, zs = (np.unique(me[i]) for i in range(3)); k = int(np.argmin(np.abs(zs))); F = {}
    for nm, row in (("glucose", 5), ("acetate", 6)):
        g = np.zeros((len(ys), len(xs))); sel = np.isclose(me[2], zs[k]); g[np.searchsorted(ys, me[1][sel]), np.searchsorted(xs, me[0][sel])] = me[row][sel]; F[nm] = g
    dead = c[lab.index("dead")] > 0
    return xs, ys, F, c[lab.index("position[0]")][~dead], c[lab.index("position[1]")][~dead]
S = [snap(t) for t in TS]; dx = S[0][0][1] - S[0][0][0]; h = S[0][0][-1] + dx / 2
cm = {"glucose": LinearSegmentedColormap.from_list("g", ["#faf8f3", "#d4a31a"]), "acetate": LinearSegmentedColormap.from_list("a", ["#faf6f6", "#cc0000"])}
plt.rcParams.update({"font.size": 11})
fig, ax = plt.subplots(2, 3, figsize=(10.5, 6.4), dpi=200, sharex=True, sharey=True, gridspec_kw=dict(wspace=0.08, hspace=0.12, right=0.88))
for r, nm in enumerate(("glucose", "acetate")):
    vmax = max(s[2][nm].max() for s in S); gm = 0
    for j, s in enumerate(S):
        a = ax[r, j]; im = a.imshow(s[2][nm], origin="lower", extent=[-h, h, -h, h], cmap=cm[nm], vmin=0, vmax=vmax * 1.0001, interpolation="nearest")
        f = s[2][nm]; a.text(0.03, 0.04, "min-max %.3f-%.3f mM" % (max(f.min(), 0.0), max(f.max(), 0.0)), transform=a.transAxes, fontsize=8, color="k", bbox=dict(fc="white", ec="none", alpha=0.7, pad=1.5))
        a.scatter(s[3], s[4], s=9, c="#1f4e79", alpha=0.7, lw=0)
        a.set_xlim(-h, h); a.set_ylim(-h, h)
        if r == 0: a.set_title("T%d  (t = %g h)" % (j + 1, TS[j] / 60))
        if r == 1: a.set_xlabel("X (µm)")
        if j == 0: a.set_ylabel("Y (µm)")
    cb = fig.colorbar(im, ax=ax[r, :], fraction=0.025, pad=0.02); cb.set_label(nm.capitalize() + " (mM)")
fig.text(0.5, 0.005, "Dots: living cells (all z, projected); field = slice z = 0, voxels 4 µm; colour scale shared along each row.", ha="center", fontsize=8)
fig.savefig(out + ".png", bbox_inches="tight", facecolor="white"); fig.savefig(out + ".pdf", bbox_inches="tight", facecolor="white")
for nm in ("glucose", "acetate"): print(nm, [("%.3f-%.3f" % (s[2][nm].min(), s[2][nm].max())) for s in S])
print("cells", [len(s[3]) for s in S])
