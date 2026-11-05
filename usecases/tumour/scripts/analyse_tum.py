"""Tumour (vessel-to-core column) analysis. usage: analyse_tum.py OUTPREFIX RUNDIR [RUNDIR ...]
Writes OUTPREFIX_series.csv (per run and snapshot), OUTPREFIX_bins.csv (cells per 10 um bin of x at 12, 24, 48 and 72 h), OUTPREFIX_cap.csv (per run, snapshot, substrate) and OUTPREFIX_runs.csv (per run, values at 12/24/48/72 h, ledger).
Geometry: the vessel is at x = 0 and x is the depth. Observables per snapshot:
  living edge        = 95th percentile of x of the living cells (the largest x, living_max, is set by one cell and varies by 40 um between seeds)
  proliferative edge = 95th percentile of x of the living cells with growth_rate > GR_MIN (half of the maximum, 0.0137 1/h, as in experiments.csv)
  dead fraction      = dead / all cells (dead cells are not removed)
Cap: voxel of a cell = (ix, iy, iz) of its position (the box has two layers in z); ratio = recorded lower bound / Michaelis-Menten bound, R = q X dt / c with X the living dry biomass in the voxel (g/L)."""
import sys, glob, os, numpy as np, pandas as pd, scipy.io as sio, xml.etree.ElementTree as ET
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..", "common")); from pcxml import labels
out = sys.argv[1]; GR_MIN = 0.5 * 0.0137; RHO = 1.04e-12; BIN = 10.0; BT = (12, 24, 48, 72)
ser, bins, cap, runs = [], [], [], []
for d in (x.rstrip("/") for x in sys.argv[2:]):
    rid = os.path.basename(d); cfg = ET.parse(f"{d}/PhysiCell_settings.xml").getroot()
    dx = float(cfg.find(".//domain/dx").text); dt = float(cfg.find(".//overall/dt_diffusion").text)
    ex = [(e.get("substrate"), e.find("fba_flux").text.strip(), float(e.find("Vmax").text), float(e.find("Km").text)) for e in cfg.findall(".//transport_model/exchange") if float(e.find("Vmax").text) > 0]
    for xf in sorted(glob.glob(f"{d}/output*[0-9].xml")):
        r = ET.parse(xf).getroot(); lab = labels(r); t = float(r.find(".//current_time").text); g = lambda k: m[lab.index(k)]
        m = sio.loadmat(xf.replace(".xml", "_cells.mat"))["cells"]; me = sio.loadmat(xf.replace(".xml", "_microenvironment0.mat"))["multiscale_microenvironment"]
        names = [v.get("name") for v in r.findall(".//microenvironment/domain/variables/variable")]
        alive = g("dead") == 0; x = g("position[0]"); gr = g("growth_rate"); grow = alive & (gr > GR_MIN)
        row = dict(run=rid, t_h=t / 60, n=int(len(x)), alive=int(alive.sum()), dead=int((~alive).sum()), dead_frac=float((~alive).mean()), living_edge=float(np.percentile(x[alive], 95)) if alive.any() else np.nan, living_max=float(x[alive].max()) if alive.any() else np.nan,
                   prolif_edge=float(np.percentile(x[grow], 95)) if grow.any() else np.nan, n_prolif=int(grow.sum()), growth_mean=float(gr[alive].mean()) if alive.any() else np.nan)
        ser.append(row)
        if round(t / 60, 3) in BT:
            e = np.arange(0, x.max() + BIN, BIN)
            for a, b in zip(e[:-1], e[1:]):
                s = (x >= a) & (x < b); bins.append(dict(run=rid, t_h=t / 60, x0=a, alive=int((s & alive).sum()), dead=int((s & ~alive).sum()), prolif=int((s & grow).sum())))
        ux, uy, uz = (np.unique(np.round(me[k], 6)) for k in range(3))
        idx = {(round(a, 6), round(b, 6), round(c, 6)): k for k, (a, b, c) in enumerate(zip(me[0], me[1], me[2]))}
        ix, iy, iz = (np.clip(np.floor((g(f"position[{k}]") - u[0]) / dx + 0.5).astype(int), 0, len(u) - 1) for k, u in enumerate((ux, uy, uz)))
        vidx = np.array([idx.get((round(ux[i], 6), round(uy[j], 6), round(uz[k], 6)), -1) for i, j, k in zip(ix, iy, iz)])
        ok = alive & (vidx >= 0); mass = g("total_volume")[ok] * (1 - g("fluid_fraction")[ok]) * RHO
        tot = np.zeros(me.shape[1]); np.add.at(tot, vidx[ok], mass); X = tot[vidx[ok]] / (me[3][vidx[ok]] * 1e-15)
        for sub, fid, V, K in ex:
            if sub not in names: continue
            c = me[4 + names.index(sub)][vidx[ok]]; lb = g(fid + "_lb")[ok]; s = c > 1e-9
            if s.sum() == 0: continue
            ratio = lb[s] / (-V * c[s] / (K + c[s])); R = (V / 60 * c[s] / (K + c[s])) * X[s] * dt / c[s]
            cap.append(dict(run=rid, t_h=t / 60, substrate=sub, n=int(s.sum()), frac_capped_pct=100 * float((ratio < 0.99).mean()), median_ratio=float(np.median(ratio)), median_R=float(np.median(R)), max_R=float(R.max())))
S, B, C = pd.DataFrame(ser), pd.DataFrame(bins), pd.DataFrame(cap); S.to_csv(out + "_series.csv", index=False); B.to_csv(out + "_bins.csv", index=False); C.to_csv(out + "_cap.csv", index=False)
for rid, s in S.groupby("run", sort=False):
    row = dict(run=rid); d = [x for x in sys.argv[2:] if os.path.basename(x.rstrip("/")) == rid][0].rstrip("/"); q = s.set_index("t_h")
    for T in BT:
        if T in q.index:
            for k in ("alive", "dead", "dead_frac", "n", "living_edge", "living_max", "prolif_edge", "n_prolif", "growth_mean"): row[f"{k}_{T}h"] = q.loc[T, k]
    row["t_end_h"] = float(s.t_h.max()); c = C[(C.run == rid) & (C.t_h >= 12)]
    if len(c): row["cap_max_R"] = float(c.max_R.max()); row["cap_max_frac_capped_pct"] = float(c.frac_capped_pct.max()); row["cap_max_R_substrate"] = c.loc[c.max_R.idxmax(), "substrate"]
    lg = pd.read_csv(f"{d}/conservation_ledger.csv").iloc[-1]
    res = {col[:-12]: abs(lg[col]) / max(lg[col[:-12] + "_cells_gross"], 1e-300) for col in lg.index if col.endswith("_cells_resid")}
    row["ledger_cells_max_rel"] = max(res.values()); row["ledger_worst"] = max(res, key=res.get)
    runs.append(row)
R = pd.DataFrame(runs); R.to_csv(out + "_runs.csv", index=False); pd.set_option("display.width", 250, "display.max_columns", 60, "display.max_rows", 200)
print(R.T.to_string(float_format=lambda v: "%.4g" % v))
