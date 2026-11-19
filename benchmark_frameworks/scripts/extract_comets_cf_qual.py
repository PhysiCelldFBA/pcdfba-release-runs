"""Criteria inputs from a COMETS CF-REF 10 um run (comets_cfref_10um.py). usage: python extract_comets_cf_qual.py <run dir> <tag>
Reads `biomass` (cycle x y model g), `media_out` (met cycle x y mmol), `flux_out` (cycle x y model_index fluxes in model order,
mmol/gDW/h) and the two .cmd model files; interior = the 120 um domain without the static ring. Writes into <run dir>:
  comets_cfref_<tag>_media_means.csv      time_min, glucose/H2/acetate/methane mean over the interior boxes (mM)
  comets_cfref_<tag>_exchange_signs.csv   time_h, species, metabolite, flux: sum over boxes of biomass x specific flux (mmol/h)
  comets_cfref_<tag>_mb_local.csv         time_min, species, biomass_g, growth_rate (1/h; biomass-weighted mean over boxes with
                                          >= 1e-15 g of ln(B(t+dt_log)/B(t)) / dt_log), H2_mM and acetate_mM (biomass-weighted
                                          mean over the species' boxes)
  comets_cfref_<tag>_checks.json          min box biomass, biomass in the ring, max box biomass / packed mass per box
"""
import sys, json, glob, numpy as np, pandas as pd

d, tag = sys.argv[1].rstrip("/"), sys.argv[2]
timing = json.load(open(glob.glob(f"{d}/comets_cfref_*_timing.json")[0]))
DT, N, dx = timing["timestep_min"], timing["grid"], timing["grid_um"] * 1e-4
V_L = timing["V_box_L"]
SP = ["c_beijerinckii", "m_barkeri"]
MET = {"glucose": "glc_D_e", "H2": "h2_e", "acetate": "ac_e", "methane": "ch4_e"}
EX = {"glucose": "EX_glc_D_e", "acetate": "EX_ac_e", "H2": "EX_h2_e", "CO2": "EX_co2_e", "methane": "EX_ch4_e"}
interior = lambda x, y: (x >= 1) & (x <= N - 2) & (y >= 1) & (y <= N - 2)

b = pd.read_csv(f"{d}/biomass", sep=r"\s+", header=None, names=["cycle", "x", "y", "model", "g"])
b["species"] = b.model.str.replace(".cmd", "", regex=False)
b["t"] = b.cycle * DT
m = pd.read_csv(f"{d}/media_out", sep=r"\s+", header=None, names=["met", "cycle", "x", "y", "mmol"])
m = m[m.met.isin(MET.values())].copy(); m["mM"] = m.mmol / V_L; m["t"] = m.cycle * DT

# media means over the interior (boxes absent from the log hold 0)
mi = m[interior(m.x, m.y)]
mm = mi.groupby(["t", "met"]).mM.sum().unstack(fill_value=0.0) / (N - 2) ** 2
pd.DataFrame({"time_min": mm.index, **{k + "_mM": mm.get(v, pd.Series(0.0, index=mm.index)).values for k, v in MET.items()}}).to_csv(
    f"{d}/comets_cfref_{tag}_media_means.csv", index=False)

# exchange fluxes: reaction order from the .cmd files
def rxn_names(cmd):
    lines = open(cmd).read().splitlines(); i = lines.index("REACTION_NAMES") + 1; out = []
    while lines[i].strip() != "//":
        out.append(lines[i].strip()); i += 1
    return out
names = [rxn_names(f"{d}/{s}.cmd") for s in SP]
rows = []
for line in open(f"{d}/flux_out"):
    v = line.split(); cyc, x, y, k = int(v[0]), int(v[1]), int(v[2]), int(v[3]) - 1
    fl = np.array(v[4:], dtype=float); assert len(fl) == len(names[k]), (len(fl), len(names[k]))
    rows.append((cyc, x, y, SP[k], {nm: fl[names[k].index(r)] for nm, r in EX.items() if r in names[k]}))
bb = b.set_index(["cycle", "x", "y", "species"]).g
ex = {}
for cyc, x, y, s, f in rows:
    w = bb.get((cyc, x, y, s), 0.0)
    for nm, val in f.items():
        ex[(cyc * DT / 60, s, nm)] = ex.get((cyc * DT / 60, s, nm), 0.0) + w * val
pd.DataFrame([dict(time_h=k[0], species=k[1], metabolite=k[2], flux=v) for k, v in sorted(ex.items())]).to_csv(
    f"{d}/comets_cfref_{tag}_exchange_signs.csv", index=False)

# local growth vs local H2 / acetate
conc = m.pivot_table(index=["t", "x", "y"], columns="met", values="mM", aggfunc="sum").fillna(0.0)
ts = sorted(b.t.unique()); out = []
for s in SP:
    q = b[b.species == s].set_index(["t", "x", "y"]).g
    for t0, t1 in zip(ts[:-1], ts[1:]):
        a = q.loc[t0]; a = a[a >= 1e-15]
        nxt = q.loc[t1].reindex(a.index).fillna(0.0)
        ok = nxt > 0
        r = np.log(nxt[ok] / a[ok]) / ((t1 - t0) / 60)
        c = conc.loc[t0].reindex(a.index).fillna(0.0)
        w = a / a.sum()
        out.append(dict(time_min=t0, species=s, biomass_g=float(q.loc[t0].sum()),
                        growth_rate=float((r * a[ok]).sum() / a[ok].sum()) if ok.any() else np.nan,
                        H2_mM=float((w * c.get(MET["H2"], 0.0)).sum()), acetate_mM=float((w * c.get(MET["acetate"], 0.0)).sum())))
pd.DataFrame(out).to_csv(f"{d}/comets_cfref_{tag}_mb_local.csv", index=False)

ring = b[~interior(b.x, b.y)]
packed = {"c_beijerinckii": 2.13e-5 * dx ** 2, "m_barkeri": 2.37e-5 * dx ** 2}
tot = b.groupby(["t", "x", "y"]).g.sum()
json.dump({"min_box_biomass_g": float(b.g.min()), "ring_biomass_g_max": float(ring.groupby("t").g.sum().max()) if len(ring) else 0.0,
           "max_box_total_over_packed_cb": float(tot.max() / packed["c_beijerinckii"]),
           "occupied_boxes_end": int((b[b.t == ts[-1]].groupby(["x", "y"]).g.sum() >= 1e-15).sum())},
          open(f"{d}/comets_cfref_{tag}_checks.json", "w"), indent=1)
print(open(f"{d}/comets_cfref_{tag}_checks.json").read())
