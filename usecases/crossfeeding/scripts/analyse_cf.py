"""Cross-feeding series analysis. usage: analyse_cf.py OUTPREFIX RUNDIR [RUNDIR ...]
Writes OUTPREFIX_series.csv (per run, snapshot and species), OUTPREFIX_cap.csv (per run, snapshot, species, substrate) and OUTPREFIX_runs.csv (per run: values at 600/1200/1800 min, the
time at which the living M. barkeri number first reaches twice its initial number, exchange signs, cap fractions at 1800 min, ledger residual / gross).
Species: 0 = c_beijerinckii (CB), 1 = m_barkeri (MB). Exchange sign = sign of the median flux of the living cells of the species (mmol/gDW/h, negative = uptake) among cells with a nonzero flux.
Cap: ratio = recorded lower bound / Michaelis-Menten bound -Vmax c/(Km+c) of the species' transport model (1 = Michaelis-Menten limited, < 0.99 = capped by the voxel), cells with c > 1e-9 mM,
as in usecases/colony/scripts/cap_check.py; R = q X dt / c with X the living dry biomass in the voxel (g/L)."""
import sys, glob, os, numpy as np, pandas as pd, scipy.io as sio, xml.etree.ElementTree as ET
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..", "common")); from pcxml import labels
out = sys.argv[1]; SP = {0: "CB", 1: "MB"}; FLUX = {"glucose": "R_EX_glc_D_e", "H2": "R_EX_h2_e", "CO2": "R_EX_co2_e", "acetate": "R_EX_ac_e", "methane": "R_EX_ch4_e"}
RHO = 1.04e-12; ser, cap, runs = [], [], []
for d in (x.rstrip("/") for x in sys.argv[2:]):
    rid = os.path.basename(d); cfg = ET.parse(f"{d}/PhysiCell_settings.xml").getroot()
    dx = float(cfg.find(".//domain/dx").text); dt = float(cfg.find(".//overall/dt_diffusion").text)
    tm = {}
    for cd in cfg.findall(".//cell_definition"):
        k = int(cd.get("ID")); tm[k] = [(e.get("substrate"), e.find("fba_flux").text.strip() if e.find("fba_flux") is not None else None, float(e.find("Vmax").text), float(e.find("Km").text))
                                       for e in cd.findall(".//transport_model/exchange") if float(e.find("Vmax").text) > 0]
    for xf in sorted(glob.glob(f"{d}/output*[0-9].xml")):
        r = ET.parse(xf).getroot(); lab = labels(r); t = float(r.find(".//current_time").text); g = lambda k: m[lab.index(k)]
        m = sio.loadmat(xf.replace(".xml", "_cells.mat"))["cells"]; me = sio.loadmat(xf.replace(".xml", "_microenvironment0.mat"))["multiscale_microenvironment"]
        names = [v.get("name") for v in r.findall(".//microenvironment/domain/variables/variable")]
        ty, alive = g("cell_type").astype(int), g("dead") == 0
        px, py = g("position[0]"), g("position[1]"); ux, uy = np.unique(np.round(me[0], 6)), np.unique(np.round(me[1], 6))
        ix = np.clip(np.rint((px - ux[0]) / dx).astype(int), 0, len(ux) - 1); iy = np.clip(np.rint((py - uy[0]) / dx).astype(int), 0, len(uy) - 1)
        vox = {(round(a, 6), round(b, 6)): k for k, (a, b) in enumerate(zip(me[0], me[1]))}; vidx = np.array([vox.get((round(ux[i], 6), round(uy[j], 6)), -1) for i, j in zip(ix, iy)])
        mass = g("total_volume") * (1 - g("fluid_fraction")) * RHO; ok = alive & (vidx >= 0); tot = np.zeros(me.shape[1]); np.add.at(tot, vidx[ok], mass[ok]); X = np.zeros(len(px)); X[ok] = tot[vidx[ok]] / (me[3][vidx[ok]] * 1e-15)
        for k, s in SP.items():
            sel = alive & (ty == k); row = dict(run=rid, t=t, species=s, alive=int(sel.sum()), dead=int(((ty == k) & ~alive).sum()), growth_mean=float(g("growth_rate")[sel].mean()) if sel.any() else np.nan)
            for nm, col in FLUX.items():
                f = g(col)[sel]; f = f[np.abs(f) > 1e-9]; row["flux_" + nm] = float(np.median(f)) if len(f) else 0.0; row["nflux_" + nm] = len(f)
            for nm in ("H2", "acetate"):
                v = vidx[sel]; v = v[v >= 0]; row["conc_" + nm] = float(np.median(me[4 + names.index(nm)][v])) if len(v) else np.nan
            ser.append(row)
            for sub, fid, V, K in tm[k]:
                if sub not in names or fid is None: continue
                q = ok & (ty == k); c = me[4 + names.index(sub)][vidx[q]]; lb = g(fid + "_lb")[q]; s_ = c > 1e-9
                if s_.sum() == 0: continue
                lbmm = -V * c[s_] / (K + c[s_]); ratio = lb[s_] / lbmm; R = (V / 60 * c[s_] / (K + c[s_])) * X[q][s_] * dt / c[s_]
                cap.append(dict(run=rid, t=t, species=s, substrate=sub, n=int(s_.sum()), frac_capped_pct=100 * float((ratio < 0.99).mean()), frac_capped5_pct=100 * float((ratio < 0.95).mean()), median_ratio=float(np.median(ratio)), median_R=float(np.median(R)), max_R=float(R.max())))
S, C = pd.DataFrame(ser), pd.DataFrame(cap); S.to_csv(out + "_series.csv", index=False); C.to_csv(out + "_cap.csv", index=False)
for rid, s in S.groupby("run", sort=False):
    row = dict(run=rid); d = [x for x in sys.argv[2:] if os.path.basename(x.rstrip("/")) == rid][0].rstrip("/")
    for sp in ("CB", "MB"):
        q = s[s.species == sp].set_index("t")
        for T in (600, 1200, 1800):
            if T in q.index: row[f"{sp}_alive_{T}"] = int(q.loc[T, "alive"]); row[f"{sp}_growth_{T}"] = q.loc[T, "growth_mean"]
        n0 = q.alive.iloc[0]; hit = q.index[q.alive >= 2 * n0]; row[f"{sp}_t_double_min"] = float(hit[0]) if len(hit) else np.nan
        up = q.index[q.alive > n0]; row[f"{sp}_t_first_new_cell_min"] = float(up[0]) if len(up) else np.nan
        last = q.iloc[-1]
        for nm in FLUX: row[f"{sp}_sign_{nm}"] = int(np.sign(last["flux_" + nm])) if last["nflux_" + nm] else 0
        row[f"{sp}_H2_voxel_1800"] = q["conc_H2"].iloc[-1]; row[f"{sp}_ac_voxel_1800"] = q["conc_acetate"].iloc[-1]
    c = C[(C.run == rid) & (C.t == C[C.run == rid].t.max())]
    for _, x in c.iterrows(): row[f"cap_{x.species}_{x.substrate}_pct"] = round(x.frac_capped_pct, 2); row[f"cap5_{x.species}_{x.substrate}_pct"] = round(x.frac_capped5_pct, 2); row[f"R_{x.species}_{x.substrate}_median"] = round(x.median_R, 4)
    for T in (600, 1200, 1800):
        for sp, sub in (("MB", "H2"), ("MB", "acetate"), ("CB", "glucose")):
            z = C[(C.run == rid) & (C.t == T) & (C.species == sp) & (C.substrate == sub)]
            if len(z): row[f"ratio_{sp}_{sub}_{T}"] = round(float(z.median_ratio.iloc[0]), 3); row[f"capped5_{sp}_{sub}_{T}_pct"] = round(float(z.frac_capped5_pct.iloc[0]), 1)
    lg = pd.read_csv(f"{d}/conservation_ledger.csv").iloc[-1]
    for col in lg.index:
        if col.endswith("_cells_resid"): b = col[:-12]; row[f"ledger_{b}"] = abs(lg[col]) / max(lg[b + "_cells_gross"], 1e-300)
    row["t_end_min"] = float(s.t.max()); runs.append(row)
R = pd.DataFrame(runs); R.to_csv(out + "_runs.csv", index=False)
pd.set_option("display.width", 250, "display.max_columns", 60); print(R.T.to_string(float_format=lambda v: '%.4g' % v))
