"""Mass-conservation figure (manuscript Fig. 2 layout) from the closed-batch unit test: A glucose and CO2 (box means), B total dry biomass (solid), analytical endpoint (dotted), cell count (dashed).
usage: make_figure.py RUNDIR OUT_PREFIX [TMAX_H=4]
Analytical endpoint (Supplementary Methods S4): initial biomass + (mu / q_glc) x glucose in the box (10 mM x 4096 um3), with mu and q_glc from FBA of the run's own model at the glucose and oxygen Vmax of its configuration (yield 0.0917 gDW per mmol glucose, 0.509 gDW/g, for the shipped e_coli_core; the 0.5197 hard-coded in dfba_analysis.ipynb is not what this model gives and overstates the endpoint by 1.85%).
Dry biomass = sum(volume x (1 - fluid fraction)) x 1.04 pg/um3."""
import sys, glob, os, numpy as np, scipy.io as sio, xml.etree.ElementTree as ET, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..", "common")); from pcxml import labels; from plot_style import SC, apply_style, FS
run, out = sys.argv[1].rstrip("/"), sys.argv[2]; TMAX = float(sys.argv[3]) if len(sys.argv) > 3 else 4.0
t, glc, co2, bio, n = [], [], [], [], []
for s in sorted(glob.glob(f"{run}/output*[0-9].xml")):
    root = ET.parse(s).getroot(); lab = labels(root); t.append(float(root.find(".//current_time").text) / 60)
    m = sio.loadmat(s.replace(".xml", "_cells.mat"))["cells"]; v, ff = m[lab.index("total_volume")], m[lab.index("fluid_fraction")]
    bio.append((v * (1 - ff)).sum() * 1.04); n.append(len(v))
    me = sio.loadmat(s.replace(".xml", "_microenvironment0.mat"))["multiscale_microenvironment"]; glc.append(me[5].mean()); co2.append(me[7].mean())
t, glc, co2, bio, n = map(np.array, (t, glc, co2, bio, n)); k = t <= TMAX
import cobra, warnings; warnings.filterwarnings("ignore")
cfg = ET.parse(f"{run}/PhysiCell_settings.xml").getroot(); uc = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
ex = {e.get("substrate"): (e.find("fba_flux").text.strip().replace("R_", "", 1), float(e.find("Vmax").text)) for e in cfg.findall(".//transport_model/exchange")}
model = cobra.io.read_sbml_model(os.path.join(uc, cfg.find(".//sbml_filename").text.strip()))
for sub in ("glucose", "oxygen"): model.reactions.get_by_id(ex[sub][0]).lower_bound = -ex[sub][1]
sol = model.optimize(); mu, qg = sol.objective_value, -sol.fluxes[ex["glucose"][0]]
print("FBA at the configured Vmax: mu %.4f 1/h, q_glc %.3f mmol/gDW/h, yield %.4f gDW/mmol (%.4f gDW/g)" % (mu, qg, mu / qg, mu / qg * 1000 / 180.0))
ana = bio[0] + (mu / qg) * (10e-3 * 4096e-15 * 1e3) * 1e12
print("final biomass %.3f pg, analytical %.3f pg (%.2f%%), cells %d, glucose at the end %.2e mM" % (bio[-1], ana, 100 * (bio[-1] / ana - 1), n[-1], glc[-1]))
FW, FH = 6.3, 2.5; apply_style(FW)
fig, (a, b) = plt.subplots(1, 2, figsize=(FW, FH), dpi=200, gridspec_kw=dict(wspace=0.75, left=0.10, right=0.90, bottom=0.30, top=0.95))   # panels side by side
a.plot(t[k], glc[k], color=SC["glucose"]); a.set_ylabel("Glucose (mM)", color=SC["glucose"]); a.set_ylim(0, 1.12 * glc.max())
a2 = a.twinx(); a2.plot(t[k], co2[k], color=SC["co2"]); a2.set_ylabel("CO$_2$ (mM)", color=SC["co2"]); a2.set_ylim(0, 1.12 * co2.max()); a2.spines["right"].set_visible(True)
l1, = b.plot(t[k], bio[k], "k-", label="Simulated total biomass"); l2 = b.axhline(ana, color="k", ls=":", lw=1.0, label="Analytical total biomass")
a.set_xlabel("Time (h)"); a.set_xlim(0, TMAX); b.set_ylabel("Total biomass (pg DW)"); b.set_ylim(0, 1.08 * max(ana, bio.max())); b.set_xlabel("Time (h)"); b.set_xlim(0, TMAX)
b2 = b.twinx(); l3, = b2.step(t[k], n[k], where="post", color="#1f78b4", ls="--", label="Simulated cell count"); b2.set_ylabel("Cell count", color="#1f78b4"); b2.set_ylim(0, 1.15 * n.max()); b2.spines["right"].set_visible(True)
fig.legend([l1, l2, l3], [l.get_label() for l in (l1, l2, l3)], loc="lower center", bbox_to_anchor=(0.5, -0.02), handlelength=1.8, columnspacing=1.6, ncol=3, frameon=False)
a.text(-0.30, 1.02, "A", transform=a.transAxes, fontsize=FS["letter"], fontweight="bold"); b.text(-0.30, 1.02, "B", transform=b.transAxes, fontsize=FS["letter"], fontweight="bold")
fig.savefig(out + ".png", bbox_inches="tight", pad_inches=0.03, facecolor="white"); fig.savefig(out + ".pdf", bbox_inches="tight", pad_inches=0.03, facecolor="white")
