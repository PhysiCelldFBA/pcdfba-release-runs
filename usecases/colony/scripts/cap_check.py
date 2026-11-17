"""Test of the cap rule on a colony run: how often is the recorded uptake bound below the Michaelis-Menten bound, and does the cap number predict it?
usage: cap_check.py RUNDIR [LABEL]     (RUNDIR holds output*.xml, *_cells.mat, *_microenvironment0.mat, PhysiCell_settings.xml)
Per snapshot and substrate with Vmax > 0, for each living cell:
  c      concentration of the substrate in the cell's voxel
  lb_MM  = -Vmax*c/(Km+c)   (mmol/gDW/h, the bound before the voxel cap)       lb = recorded R_EX_*_lb (the bound the LP used)
  ratio  = lb/lb_MM         (1 = Michaelis-Menten limited, < 1 = the cap acted)
  R      = q(c) * X * dt / c with q = Vmax/60 * c/(Km+c) mmol/gDW/min and X the dry biomass of the living cells in that voxel per voxel volume (g/L)
Reports the fraction of cells with ratio < 0.99, the median ratio, the median and the largest R, and the median ratio in bins of R.  Writes results/cap_check_<LABEL>.csv"""
import sys, glob, os, numpy as np, pandas as pd, scipy.io as sio, xml.etree.ElementTree as ET
d = sys.argv[1].rstrip("/"); label = sys.argv[2] if len(sys.argv) > 2 else os.path.basename(d)
cfg = ET.parse(f"{d}/PhysiCell_settings.xml").getroot()
dx = float(cfg.find(".//domain/dx").text); dt = float(cfg.find(".//overall/dt_diffusion").text)
ex = [(e.get("substrate"), e.find("fba_flux").text.strip(), float(e.find("Vmax").text), float(e.find("Km").text)) for e in cfg.findall(".//transport_model/exchange") if float(e.find("Vmax").text) > 0]
rho = 1.04e-12                                                                                     # g/um3
rows, bins = [], []
for xf in sorted(glob.glob(f"{d}/output*[0-9].xml")):
    r = ET.parse(xf).getroot(); t = float(r.find(".//current_time").text)
    lab = {l.text: (int(l.get("index")), int(l.get("size"))) for l in r.findall(".//cell_population/custom/simplified_data/labels/label")}
    m = sio.loadmat(xf.replace(".xml", "_cells.mat"))["cells"]; me = sio.loadmat(xf.replace(".xml", "_microenvironment0.mat"))["multiscale_microenvironment"]
    names = [v.get("name") for v in r.findall(".//microenvironment/domain/variables/variable")]
    col = lambda k: m[lab[k][0]]
    x, y = col("position[0]") if "position[0]" in lab else m[lab["position"][0]], None
    pos = m[lab["position"][0]:lab["position"][0] + 3] if "position" in lab else None
    if pos is None: pos = np.vstack([m[lab[f"position[{i}]"][0]] for i in range(3)])
    alive = col("dead") == 0
    ux, uy = np.unique(np.round(me[0], 6)), np.unique(np.round(me[1], 6))
    ix = np.clip(np.searchsorted(ux, pos[0] - 1e-9) - 0, 0, len(ux) - 1); iy = np.clip(np.searchsorted(uy, pos[1] - 1e-9), 0, len(uy) - 1)
    ix = np.clip(np.rint((pos[0] - ux[0]) / dx).astype(int), 0, len(ux) - 1); iy = np.clip(np.rint((pos[1] - uy[0]) / dx).astype(int), 0, len(uy) - 1)
    vox_of = {(round(a, 6), round(b, 6)): k for k, (a, b) in enumerate(zip(me[0], me[1]))}
    vidx = np.array([vox_of.get((round(ux[i], 6), round(uy[j], 6)), -1) for i, j in zip(ix, iy)])
    ok = alive & (vidx >= 0); vol = me[3][vidx[ok]]
    mass = col("total_volume")[ok] * (1 - col("fluid_fraction")[ok]) * rho                          # g per cell
    tot = np.zeros(me.shape[1]); np.add.at(tot, vidx[ok], mass); X = tot[vidx[ok]] / (vol * 1e-15)  # g/L in the cell's voxel
    for s, fid, V, K in ex:
        if s not in names: continue
        c = me[4 + names.index(s)][vidx[ok]]; lb = col(fid + "_lb")[ok]
        sel = c > 1e-9
        if sel.sum() == 0: continue
        lbmm = -V * c[sel] / (K + c[sel]); ratio = lb[sel] / lbmm; R = (V / 60 * c[sel] / (K + c[sel])) * X[sel] * dt / c[sel]
        rows.append(dict(t=t, substrate=s, n=int(sel.sum()), frac_capped_pct=100 * float((ratio < 0.99).mean()), median_ratio=float(np.median(ratio)), median_R=float(np.median(R)), max_R=float(R.max())))
        edges = [0, 0.01, 0.05, 0.2, 1, 5, np.inf]
        for a, b in zip(edges[:-1], edges[1:]):
            q = (R >= a) & (R < b)
            if q.sum() > 20: bins.append(dict(t=t, substrate=s, R_bin=f"{a}-{b}", n=int(q.sum()), median_ratio=float(np.median(ratio[q])), frac_capped_pct=100 * float((ratio[q] < 0.99).mean())))
os.makedirs("results", exist_ok=True); A, B = pd.DataFrame(rows), pd.DataFrame(bins)
A.to_csv(f"results/cap_check_{label}.csv", index=False); B.to_csv(f"results/cap_check_{label}_bins.csv", index=False)
pd.set_option("display.width", 220); print(f"== {label}: dx {dx}, dt {dt}, exchanges {[e[0] for e in ex]}"); print(A.round(3).to_string(index=False))
