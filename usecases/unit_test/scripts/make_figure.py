"""Mass-conservation figure (manuscript Fig. 2 layout) from the closed-batch unit test: A glucose and CO2 (box means), B total dry biomass (solid), analytical endpoint (dotted), cell count (dashed).
usage: make_figure.py RUNDIR OUT_PREFIX [TMAX_H=4]
Analytical endpoint as in dfba_analysis.ipynb: initial biomass + glucose in the box (10 mM x 4096 um3) x 180 g/mol x 0.5197 gDW/g glucose (the e_coli_core aerobic yield, taken from the notebook, not recomputed here).
Dry biomass = sum(volume x (1 - fluid fraction)) x 1.04 pg/um3."""
import sys, glob, os, numpy as np, scipy.io as sio, xml.etree.ElementTree as ET, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..", "common")); from pcxml import labels
run, out = sys.argv[1].rstrip("/"), sys.argv[2]; TMAX = float(sys.argv[3]) if len(sys.argv) > 3 else 4.0
t, glc, co2, bio, n = [], [], [], [], []
for s in sorted(glob.glob(f"{run}/output*[0-9].xml")):
    root = ET.parse(s).getroot(); lab = labels(root); t.append(float(root.find(".//current_time").text) / 60)
    m = sio.loadmat(s.replace(".xml", "_cells.mat"))["cells"]; v, ff = m[lab.index("total_volume")], m[lab.index("fluid_fraction")]
    bio.append((v * (1 - ff)).sum() * 1.04); n.append(len(v))
    me = sio.loadmat(s.replace(".xml", "_microenvironment0.mat"))["multiscale_microenvironment"]; glc.append(me[5].mean()); co2.append(me[7].mean())
t, glc, co2, bio, n = map(np.array, (t, glc, co2, bio, n)); k = t <= TMAX
ana = bio[0] + 10e-3 * 4096e-15 * 180 * 0.5197 * 1e12
print("final biomass %.3f pg, analytical %.3f pg (%.2f%%), cells %d, glucose at the end %.2e mM" % (bio[-1], ana, 100 * (bio[-1] / ana - 1), n[-1], glc[-1]))
plt.rcParams.update({"font.size": 11, "axes.spines.top": False})
fig, (a, b) = plt.subplots(2, 1, figsize=(4.4, 6.4), sharex=True, dpi=200, gridspec_kw=dict(hspace=0.12))
a.plot(t[k], glc[k], color="#f1c232", lw=2.2); a.set_ylabel("Glucose (mM)", color="#f1c232"); a.set_ylim(0, 1.12 * glc.max())
a2 = a.twinx(); a2.plot(t[k], co2[k], color="#8e44ad", lw=2.2); a2.set_ylabel("CO$_2$ (mM)", color="#8e44ad"); a2.set_ylim(0, 1.12 * co2.max()); a2.spines["top"].set_visible(False)
l1, = b.plot(t[k], bio[k], "k-", lw=2.2, label="Simulated total biomass"); l2 = b.axhline(ana, color="k", ls=":", lw=1.6, label="Analytical total biomass")
b.set_ylabel("Total biomass (pg DW)"); b.set_ylim(0, 1.08 * max(ana, bio.max())); b.set_xlabel("Time (h)"); b.set_xlim(0, TMAX)
b2 = b.twinx(); l3, = b2.step(t[k], n[k], where="post", color="#1f78b4", ls="--", lw=2.2, label="Simulated cell count"); b2.set_ylabel("Cell count", color="#1f78b4"); b2.set_ylim(0, 1.15 * n.max()); b2.spines["top"].set_visible(False)
b.legend([l1, l2, l3], [l.get_label() for l in (l1, l2, l3)], loc="upper center", bbox_to_anchor=(0.5, -0.22), frameon=False)
a.text(-0.2, 1.02, "A", transform=a.transAxes, fontsize=20, fontweight="bold"); b.text(-0.2, 1.02, "B", transform=b.transAxes, fontsize=20, fontweight="bold")
fig.savefig(out + ".png", bbox_inches="tight", facecolor="white"); fig.savefig(out + ".pdf", bbox_inches="tight", facecolor="white")
