"""Acetate-switch figure from ACE-REF/ACE-SEED (5 mM pulse) and ACE-PULSE (10 mM), dx 4, dt 0.01, 3 seeds: box-mean concentrations and total biomass, after the layout of the
original main figure. A: 5 mM O2 pulse at 210 min (the original protocol); B: 10 mM pulse, 900 min.  /opt/anaconda3/bin/python"""
import sys, os, glob, numpy as np, scipy.io as sio, xml.etree.ElementTree as ET, matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..", "common")); from pcxml import labels
HERE = os.path.dirname(os.path.abspath(__file__)); UC = os.path.join(HERE, ".."); R = UC + "/runs"; OUT = UC + "/results"
EDGE_UM = 36.0; V_L = EDGE_UM ** 3 * 1e-15; DRY = 0.25 * 1.04e-12   # g per um^3 of cell volume (fluid fraction 0.75, 1.04 g/mL)
def series(run):
    T, C, B, J = [], [], [], []
    for s in sorted(glob.glob(f"{R}/{run}/output0*.xml")):
        if "_" in os.path.basename(s)[6:]: continue
        root = ET.parse(s).getroot(); L = labels(root); m = sio.loadmat(s.replace(".xml", "_cells.mat"))["cells"]
        a = m[L.index("dead")] == 0
        T.append(float(root.find(".//current_time").text) / 60)
        B.append(m[L.index("total_volume")][a].sum() * DRY / V_L)
        J.append(np.median(m[L.index("R_EX_ac_e")][a]) if a.any() else 0.0)
        me = sio.loadmat(s.replace(".xml", "_microenvironment0.mat"))["multiscale_microenvironment"]
        C.append([me[4 + i].mean() for i in range(4)])
    return np.array(T), np.array(C), np.array(B), np.array(J)
COL = dict(Glucose="#f0b800", Acetate="#cc0000", Oxygen="#3b82f6", CO2="#8e44ad")
def panel(ax, runs, tmax_h, mark):
    S = [series(r) for r in runs]; T = S[0][0]
    C = np.stack([s[1] for s in S]); B = np.stack([s[2] for s in S]); J = np.stack([s[3] for s in S])
    for i, (nm, k) in enumerate((("Oxygen", 0), ("Glucose", 1), ("Acetate", 2), ("CO2", 3))):
        ax.plot(T, C[:, :, k].mean(0), color=COL[nm], lw=2, label=nm if nm != "CO2" else "CO$_2$")
        ax.fill_between(T, C[:, :, k].min(0), C[:, :, k].max(0), color=COL[nm], alpha=0.25, lw=0)
    ax.set_xlim(0, tmax_h); ax.set_xlabel("Time (h)"); ax.set_ylabel("Mean concentration (mM)")
    ax2 = ax.twinx(); ax2.plot(T, B.mean(0), color="k", lw=2.5, label="Total biomass"); ax2.set_ylabel("Total biomass (g/L)"); ax2.set_ylim(0, B.max() * 1.15)
    for sp in ("top",): ax.spines[sp].set_visible(False); ax2.spines[sp].set_visible(False)
    t = T; j = J.mean(0); ac = C[:, :, 2].mean(0); glc = C[:, :, 1].mean(0)
    t1 = t[np.argmax(j > 0.5)]; t2 = t[np.argmax(glc < 0.05)]; t3 = 210 / 60; t4 = t[(t > t3 + 0.1) & (np.abs(j) < 1e-6)][0]
    for tt, lab in ((t1, "$t_1$"), (t3, "$t_3$"), (t4, "$t_4$")):
        ax.axvline(tt, color="grey", ls=":", lw=1); ax.text(tt + 0.05 * tmax_h / 6, ax.get_ylim()[1] * 0.96, lab, color="#444", va="top")
    return ax2, dict(t1=t1, t2=t2, t3=t3, t4=t4, x0=B[:, 0].mean(), xmax=B.max(), ac_end=ac[-1])
fig, axs = plt.subplots(2, 1, figsize=(9, 8.5))
a2, i1 = panel(axs[0], ["ACE-REF_s50", "ACE-SEED_s51", "ACE-SEED_s52"], 10, "")
b2, i2 = panel(axs[1], [f"ACE-PULSE_o2_10_s{s}" for s in (50, 51, 52)], 15, "")
h, l = axs[0].get_legend_handles_labels(); h2, l2 = a2.get_legend_handles_labels()
fig.legend(h + h2, l + l2, ncol=5, loc="upper center", frameon=False)
axs[0].set_title("A   5 mM O$_2$ pulse at 3.5 h", loc="left", fontsize=11); axs[1].set_title("B   10 mM O$_2$ pulse at 3.5 h", loc="left", fontsize=11)
fig.tight_layout(rect=(0, 0, 1, 0.95)); fig.savefig(OUT + "/fig_acetate_switch.png", dpi=200); fig.savefig(OUT + "/fig_acetate_switch.pdf")
print({k: round(float(v), 3) for k, v in i1.items()}); print({k: round(float(v), 3) for k, v in i2.items()})
