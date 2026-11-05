"""Supplementary Figure S6: total biomass (log scale) of a quarter-colony run (x4 = full colony) and the two-line fit of ln(biomass) with the breakpoint of the smallest combined residual (as in make_figure.py).
usage: make_logfigure.py RUNDIR OUT_PREFIX      prints the breakpoint, the doubling times of the two lines and the specific growth rate d ln(biomass)/dt (centred differences, 1 h snapshots)
Biomass conventions as in make_figure.py: dry weight = volume x 0.75 x 1.04 pg/um3, x4 for the full colony."""
import sys, glob, os, numpy as np, scipy.io as sio, xml.etree.ElementTree as ET, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..", "common")); from pcxml import labels
run, out = sys.argv[1].rstrip("/"), sys.argv[2]; FULL = 4; DRY = 0.75 * 1.04
ts, bio = [], []
for s in sorted(glob.glob(f"{run}/output*[0-9].xml")):
    r = ET.parse(s).getroot(); lab = labels(r); m = sio.loadmat(s.replace(".xml", "_cells.mat"))["cells"]
    ts.append(float(r.find(".//current_time").text) / 60); bio.append(FULL * m[lab.index("total_volume")].sum() * DRY / 1e6)
ts, bio = np.array(ts), np.array(bio); ok = bio > 0; t, y = ts[ok], np.log(bio[ok]); best = None
for k in range(3, len(t) - 3):
    c1, c2 = np.polyfit(t[:k + 1], y[:k + 1], 1), np.polyfit(t[k:], y[k:], 1)
    r = ((y[:k + 1] - np.polyval(c1, t[:k + 1])) ** 2).sum() + ((y[k:] - np.polyval(c2, t[k:])) ** 2).sum()
    if best is None or r < best[0]: best = (r, k, c1, c2)
_, k, c1, c2 = best; Tb = t[k]
print("breakpoint %.1f h; doubling time %.2f h before, %.2f h after; biomass %.4f ug at %.1f h" % (Tb, np.log(2) / c1[0], np.log(2) / c2[0], bio[-1], ts[-1]))
mu = np.gradient(y, t); sel = lambda a, b: mu[(t >= a) & (t <= b)]
print("specific growth rate (1/h): 1-3 h %.3f, 3-5 h %.3f, 5-6 h %.3f, 6-8 h %.3f, 8-10 h %.3f" % (sel(1, 3).mean(), sel(3, 5).mean(), sel(5, 6).mean(), sel(6, 8).mean(), sel(8, 10).mean()))
print("exponential-phase growth rate from the first line %.3f 1/h; at the end %.3f 1/h" % (c1[0], mu[-1]))
plt.rcParams.update({"font.size": 12}); fig, ax = plt.subplots(figsize=(6.2, 4.6), dpi=200)
ax.semilogy(ts[ok], bio[ok], "o", c="g", mec="k", ms=5, label="Biomass")
ax.semilogy(t[:k + 1], np.exp(np.polyval(c1, t[:k + 1])), "-", c="#1f77b4", lw=2, label=f"Early fit (doubling {np.log(2) / c1[0]:.2f} h)")
ax.semilogy(t[k:], np.exp(np.polyval(c2, t[k:])), "-", c="#d62728", lw=2, label=f"Late fit (doubling {np.log(2) / c2[0]:.2f} h)")
ax.axvline(Tb, c="k", ls="--", lw=1.2); ax.text(Tb + 0.12, ax.get_ylim()[0] * 1.6, r"$t^*\approx$%.1f h" % Tb, fontsize=11, ha="left")
ax.set_xlabel("Time (h)"); ax.set_ylabel(r"Total biomass ($\mu$g)"); ax.legend(frameon=False, loc="upper left", fontsize=10)
for a in ("top", "right"): ax.spines[a].set_visible(False)
fig.savefig(out + ".png", bbox_inches="tight", facecolor="white"); fig.savefig(out + ".pdf", bbox_inches="tight", facecolor="white"); print("wrote", out)
