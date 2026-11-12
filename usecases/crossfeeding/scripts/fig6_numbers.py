"""Numbers quoted in the Fig. 6 text of the cross-feeding use case, from one run, at the requested times (h).
usage: fig6_numbers.py RUNDIR T_H [T_H ...]
Definitions as in make_figure.py: living cells; total net export (fmol/s) = sum(flux x dry weight) x 1e12 / 3600; per biomass = sum(flux x dw) / sum(dw) (mmol/gDW/h);
growth = mean over the living cells of a species; box means over all voxels of the mesh; positive = export."""
import sys, glob, os, numpy as np, scipy.io as sio, xml.etree.ElementTree as ET
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..", "common")); from pcxml import labels
run = sys.argv[1].rstrip("/"); TS = [float(x) for x in sys.argv[2:]]
EXP = {"glucose": "R_EX_glc_D_e", "H2": "R_EX_h2_e", "CO2": "R_EX_co2_e", "acetate": "R_EX_ac_e", "CH4": "R_EX_ch4_e"}; VAR = {"glucose": 4, "methane": 5, "H2": 7, "acetate": 11}
files = sorted(glob.glob(f"{run}/output*[0-9].xml")); tof = lambda f: float(ET.parse(f).getroot().find(".//current_time").text) / 60
for T in TS:
    f = min(files, key=lambda f: abs(tof(f) - T)); root = ET.parse(f).getroot(); lab = labels(root); c = sio.loadmat(f.replace(".xml", "_cells.mat"))["cells"]; g = lambda k: c[lab.index(k)]
    ty, v, ff, dead = g("cell_type"), g("total_volume"), g("fluid_fraction"), g("dead") > 0; dw = v * (1 - ff) * 1.04e-12; al = ~dead
    me = sio.loadmat(f.replace(".xml", "_microenvironment0.mat"))["multiscale_microenvironment"]; xs, ys = np.unique(me[0]), np.unique(me[1])
    print("== t = %.2f h (%s)" % (tof(f), os.path.basename(f)))
    print("   living cells CB %d, MB %d; dead %d; all %d" % (((ty == 0) & al).sum(), ((ty == 1) & al).sum(), dead.sum(), len(ty)))
    print("   mean growth rate CB %.4f, MB %.4f /h" % tuple(g("growth_rate")[(ty == k) & al].mean() for k in (0, 1)))
    print("   box means (mM): " + ", ".join("%s %.3f" % (k, me[r].mean()) for k, r in VAR.items()))
    for k, nm in ((0, "CB"), (1, "MB")):
        s = (ty == k) & al
        print("   %s total export (fmol/s): " % nm + ", ".join("%s %.2f" % (n, (g(col) * dw * 1e12 / 3600)[s].sum()) for n, col in EXP.items()) + "   | per biomass (mmol/gDW/h): " + ", ".join("%s %.2f" % (n, (g(col) * dw)[s].sum() / dw[s].sum()) for n, col in EXP.items()))
    G = {k: np.zeros((len(ys), len(xs))) for k in ("H2", "acetate", "methane")}
    for k in G: G[k][np.searchsorted(ys, me[1]), np.searchsorted(xs, me[0])] = me[VAR[k]]
    ix = np.clip(np.searchsorted(xs, g("position[0]")), 0, len(xs) - 1); iy = np.clip(np.searchsorted(ys, g("position[1]")), 0, len(ys) - 1); mb = (ty == 1) & al; cb = (ty == 0) & al
    inner = lambda A: A[1:-1, 1:-1]
    print("   maxima H2 %.3f, acetate %.3f, CH4 %.3f (all voxels); without the outer layer: %.3f, %.3f, %.3f" % (G["H2"].max(), G["acetate"].max(), G["methane"].max(), inner(G["H2"]).max(), inner(G["acetate"]).max(), inner(G["methane"]).max()))
    cx, cy = len(xs) // 2, len(ys) // 2
    print("   H2 at the centre voxel %.3f, mean over 5x5 centre %.3f; H2 in the voxels of MB cells: mean %.3f; of CB cells: mean %.3f; CH4 in the voxels of MB cells: mean %.3f, max %.3f; CH4 at centre %.3f" % (G["H2"][cy, cx], G["H2"][cy-2:cy+3, cx-2:cx+3].mean(), G["H2"][iy[mb], ix[mb]].mean(), G["H2"][iy[cb], ix[cb]].mean(), G["methane"][iy[mb], ix[mb]].mean(), G["methane"][iy[mb], ix[mb]].max(), G["methane"][cy, cx]))
    print("   cells outside +-100 um: %d; max |position| %.1f" % ((np.maximum(abs(g("position[0]")), abs(g("position[1]"))) > 100).sum(), np.maximum(abs(g("position[0]")), abs(g("position[1]"))).max()))
