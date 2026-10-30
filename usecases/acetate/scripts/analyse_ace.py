"""Acetate runs on the pin: regime times, concentration observables and uptake-cap fraction per run.
usage: analyse_ace.py [RUN_ID ...]  (default: every runs/ACE-*)  -> results/ace_summary.csv, results/ace_series.csv
Per snapshot (all alive cells): median exchange flux of acetate (R_EX_ac_e, + = secretion), box-mean concentrations, dry biomass per box volume (g/L),
and for each exchange with Vmax > 0 the fraction of cells whose recorded bound (R_EX_*_lb) is below the Michaelis-Menten bound (ratio < 0.99, the cap acted)
and the cap number R = q X dt / c (q at the voxel concentration, X dry biomass of the living cells in the voxel per voxel volume, as in colony/scripts/cap_check.py).
Per run (cap columns: snapshots with t > 0 and box-mean concentration > 0.1 mM only): ferm_start / ferm_end = first / last time the median acetate flux > 0.5; reuse_start = first time < -0.5 (5 min resolution);
ac_peak = maximum box-mean acetate; ac_floor = minimum box-mean acetate after the peak; biomass_end in g/L."""
import sys, os, glob, numpy as np, pandas as pd, scipy.io as sio, xml.etree.ElementTree as ET
HERE = os.path.dirname(os.path.abspath(__file__)); R = os.path.join(HERE, "..", "runs"); OUT = os.path.join(HERE, "..", "results")
names = sys.argv[1:] or sorted(n for n in os.listdir(R) if n.startswith("ACE-"))
RHO = 1.04e-12
series, summ = [], []
for run in names:
    d = f"{R}/{run}"; cfg = ET.parse(f"{d}/PhysiCell_settings.xml").getroot()
    dx = float(cfg.find(".//domain/dx").text); dt = float(cfg.find(".//overall/dt_diffusion").text)
    ex = [(e.get("substrate"), e.find("fba_flux").text.strip(), float(e.find("Vmax").text), float(e.find("Km").text)) for e in cfg.findall(".//transport_model/exchange") if float(e.find("Vmax").text) > 0]
    for xf in sorted(f for f in glob.glob(f"{d}/output*[0-9].xml")):
        r = ET.parse(xf).getroot(); t = float(r.find(".//current_time").text)
        lab = {l.text: (int(l.get("index")), int(l.get("size"))) for l in r.findall(".//cell_population/custom/simplified_data/labels/label")}
        m = sio.loadmat(xf.replace(".xml", "_cells.mat"))["cells"]; me = sio.loadmat(xf.replace(".xml", "_microenvironment0.mat"))["multiscale_microenvironment"]
        names_s = [v.get("name") for v in r.findall(".//microenvironment/domain/variables/variable")]
        col = lambda k: m[lab[k][0]]; pos = np.vstack([col(f"position[{i}]") if f"position[{i}]" in lab else m[lab["position"][0] + i] for i in range(3)])
        alive = col("dead") == 0; ok = alive
        axes = [np.unique(np.round(me[i], 6)) for i in range(3)]
        idx = [np.clip(np.rint((pos[i] - axes[i][0]) / dx).astype(int), 0, len(axes[i]) - 1) for i in range(3)]
        vox = {(round(a, 6), round(b, 6), round(c, 6)): k for k, (a, b, c) in enumerate(zip(me[0], me[1], me[2]))}
        vidx = np.array([vox.get((round(axes[0][i], 6), round(axes[1][j], 6), round(axes[2][k], 6)), -1) for i, j, k in zip(*idx)])
        ok = alive & (vidx >= 0); v = vidx[ok]; volv = me[3]
        mass = col("total_volume")[ok] * (1 - col("fluid_fraction")[ok]) * RHO; tot = np.zeros(me.shape[1]); np.add.at(tot, v, mass)
        X = tot[v] / (volv[v] * 1e-15)
        row = dict(run=run, t=t, n_alive=int(alive.sum()), biomass_gL=float(mass.sum() / (volv.sum() * 1e-15)), ac_flux_med=float(np.median(col("R_EX_ac_e")[alive])) if alive.any() else np.nan)
        for s in ("oxygen", "glucose", "acetate"):
            if s in names_s: row[f"{s}_mean"] = float(me[4 + names_s.index(s)].mean())
        for s, fid, V, K in ex:
            if s not in names_s or not len(v): continue
            c = me[4 + names_s.index(s)][v]; lb = col(fid + "_lb")[ok]; sel = c > 1e-9
            if sel.sum() == 0: continue
            lbmm = -V * c[sel] / (K + c[sel]); ratio = lb[sel] / lbmm; Rn = (V / 60 * c[sel] / (K + c[sel])) * X[sel] * dt / c[sel]
            row[f"{s}_capped_pct"] = 100 * float((ratio < 0.99).mean()); row[f"{s}_R_med"] = float(np.median(Rn)); row[f"{s}_R_max"] = float(Rn.max())
        series.append(row)
S = pd.DataFrame(series); S.to_csv(f"{OUT}/ace_series.csv", index=False)
for run, g in S.groupby("run", sort=False):
    f = g.ac_flux_med.values; t = g.t.values; sec = t[f > 0.5]; reu = t[f < -0.5]; ip = int(np.nanargmax(g.acetate_mean.values))
    row = dict(run=run, snaps=len(g), t_end=t[-1], n_end=int(g.n_alive.iloc[-1]), biomass_end_gL=g.biomass_gL.iloc[-1],
               ferm_start=sec.min() if len(sec) else np.nan, ferm_end=sec.max() if len(sec) else np.nan, reuse_start=reu.min() if len(reu) else np.nan,
               ac_peak=g.acetate_mean.max(), t_ac_peak=t[ip], ac_floor=g.acetate_mean.values[ip:].min(), ac_end=g.acetate_mean.iloc[-1])
    for s in ("oxygen", "glucose", "acetate"):
        if f"{s}_capped_pct" in g:   # only while the substrate is present (box mean > 0.1 mM, t > 0): at exhaustion the voxel limit is the physical limit, and the t = 0 snapshot has no bound yet
            h = g[(g.t > 0) & (g[f"{s}_mean"] > 0.1)]
            row[f"{s}_capped_max_pct"] = h[f"{s}_capped_pct"].max() if len(h) else np.nan; row[f"{s}_R_med_max"] = h[f"{s}_R_med"].max() if len(h) else np.nan; row[f"{s}_R_max"] = h[f"{s}_R_max"].max() if len(h) else np.nan
    summ.append(row)
D = pd.DataFrame(summ); D.to_csv(f"{OUT}/ace_summary.csv", index=False)
pd.set_option("display.width", 250); pd.set_option("display.max_columns", 40); print(D.round(3).to_string(index=False))
