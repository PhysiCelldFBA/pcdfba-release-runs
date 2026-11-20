"""Snapshot fields of BM-CF-COMETS for the microenvironment figures of S12 (make_microenv.py).
usage: python3 extract_comets_fields.py <run dir> <out.npz>
Reads `media_out` (met cycle x y mmol) and `biomass` (cycle x y model g) of the run and the grid of its *_timing.json, and writes
  t_h                       (T,)            times of the media log (every MediaLogRate cycles)
  glucose, H2, acetate, methane (T, N, N)  concentration per box (mM = mmol / V_box_L); boxes absent from the log hold 0
  c_beijerinckii, m_barkeri (T, N, N)      biomass per box (g)
  grid_um, N                               box width (um) and boxes per side, including the static ring (index 0 and N-1)
Arrays are indexed [t, x, y] as in COMETS."""
import sys, json, glob, numpy as np, pandas as pd

d, out = sys.argv[1].rstrip("/"), sys.argv[2]
tm = json.load(open(glob.glob(f"{d}/comets_cfref_*_timing.json")[0]))
DT, N, V_L = tm["timestep_min"], tm["grid"], tm["V_box_L"]
MET = {"glucose": "glc_D_e", "H2": "h2_e", "acetate": "ac_e", "methane": "ch4_e"}

m = pd.read_csv(f"{d}/media_out", sep=r"\s+", header=None, names=["met", "cycle", "x", "y", "mmol"])
m = m[m.met.isin(MET.values())].copy()
# media_out is 1-based (static ring at 1 and N; checked 2026-11-20: H2 at cycle 0 in 2..N-1, glucose 1 mM at 1 and N), biomass 0-based
assert m.x.min() == 1 and m.y.min() == 1, "media_out expected 1-based"
m["x"] -= 1; m["y"] -= 1
b = pd.read_csv(f"{d}/biomass", sep=r"\s+", header=None, names=["cycle", "x", "y", "model", "g"])
b["species"] = b.model.str.replace(".cmd", "", regex=False)
assert m.x.max() < N and m.y.max() < N and b.x.max() < N and b.y.max() < N, "coordinates outside the grid"
cyc = np.array(sorted(m.cycle.unique())); ci = {c: i for i, c in enumerate(cyc)}
res = {"t_h": cyc * DT / 60.0, "grid_um": tm["grid_um"], "N": N}
for k, v in MET.items():
    a = np.zeros((len(cyc), N, N)); q = m[m.met == v]
    a[q.cycle.map(ci).values, q.x.values, q.y.values] = q.mmol.values / V_L; res[k] = a
for s in ("c_beijerinckii", "m_barkeri"):
    a = np.zeros((len(cyc), N, N)); q = b[(b.species == s) & b.cycle.isin(ci)]
    a[q.cycle.map(ci).values, q.x.values, q.y.values] = q.g.values; res[s] = a
np.savez_compressed(out, **res)
print(out, len(cyc), "snapshots", "x range", m.x.min(), m.x.max(), "biomass cycles", b.cycle.nunique())
