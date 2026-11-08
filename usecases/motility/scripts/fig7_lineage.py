"""Fig. 7 panels C/E (high supply) as a lineage plot: one trace per cell (ID), migration speed, growth rate and total ATP flux against time, with a fork at every division (the daughter's trace starts at the division time and is joined to its parent, the nearest cell of the previous snapshot). usage: fig7_lineage.py RUN_DIR [OUT_PREFIX=results/fig7/fig7_lineage]"""
import sys, os, glob, numpy as np, scipy.io as sio, xml.etree.ElementTree as ET, matplotlib as mpl
mpl.use("Agg"); import matplotlib.pyplot as plt
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "../../../common")); from pcxml import labels
d = sys.argv[1].rstrip("/"); out = sys.argv[2] if len(sys.argv) > 2 else "results/fig7/fig7_lineage"
snaps = sorted(glob.glob(f"{d}/output[0-9]*[0-9].xml")); T = {}; prev = None; parent = {}; birth = {}
for f in snaps:
    root = ET.parse(f).getroot(); lab = labels(root); t = float(root.find(".//current_time").text); m = sio.loadmat(f.replace(".xml", "_cells.mat"))["cells"]
    col = lambda k: m[lab.index(k)]; ids = col("ID").astype(int); pos = np.c_[col("position[0]"), col("position[1]")]
    for j, i in enumerate(ids):
        if i not in T:
            T[i] = []; birth[i] = t
            if prev is not None:
                pid, ppos = prev; parent[i] = int(pid[np.argmin(((ppos - pos[j]) ** 2).sum(1))])
        T[i].append((t, col("migration_speed")[j], col("growth_rate")[j], col("atp_flux")[j]))
    prev = (ids, pos)
order = sorted(T, key=lambda i: birth[i]); gen = {}; 
for i in order: gen[i] = 0 if i not in parent else gen[parent[i]] + 1
cm = plt.get_cmap("viridis"); fig, axs = plt.subplots(3, 1, figsize=(4.6, 5.6), sharex=True)
for ax, k, lab in zip(axs, (1, 2, 3), ("Migration speed (µm/min)", "Growth rate (h$^{-1}$)", "Total ATP flux\n(mmol/gDW/h)")):
    for i in order:
        a = np.array(T[i]); c = cm(min(gen[i], 3) / 4); t0 = a[0, 0] / 60
        ax.plot(a[:, 0] / 60, a[:, k], color=c, lw=1.6 if gen[i] == 0 else 1.1)
        if i in parent:   # fork: connect the first value of the daughter to the parent's value at the same time
            pa = np.array(T[parent[i]]); v = np.interp(a[0, 0], pa[:, 0], pa[:, k]); ax.plot([t0, t0], [v, a[0, k]], color=c, lw=0.8, ls=":")
    ax.set_ylabel(lab); ax.spines[["top", "right"]].set_visible(False)
axs[-1].set_xlabel("Time (h)"); fig.tight_layout(); fig.savefig(out + ".png", dpi=300); fig.savefig(out + ".svg")
print(len(order), "cells; generations", sorted(set(gen.values())), "divisions at (min)", sorted(round(birth[i]) for i in parent)); print("wrote", out)
