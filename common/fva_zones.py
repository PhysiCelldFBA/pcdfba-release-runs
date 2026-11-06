"""Share of cells with an undetermined exchange (FVA range > TOL) by distance bin, from a fva_replay *_cells.csv. usage: fva_zones.py CELLS.csv [AXIS=x] [EDGES=0,60,120,160,200,400] [FRACTION=0.99999] [TOL=1e-2] [OUT.csv]
Distance is the position along AXIS (x for the tumour: distance from the vessel at x = 0). Only feasible cells enter; the number of infeasible cells is printed. Writes the table if OUT.csv is given."""
import sys, pandas as pd
f = sys.argv[1]; ax = sys.argv[2] if len(sys.argv) > 2 else "x"; ed = [float(v) for v in (sys.argv[3] if len(sys.argv) > 3 else "0,60,120,160,200,400").split(",")]
fr = float(sys.argv[4]) if len(sys.argv) > 4 else 0.99999; tol = float(sys.argv[5]) if len(sys.argv) > 5 else 1e-2
d = pd.read_csv(f); d = d[d.fraction == fr]; inf = d[d.obj_fba.isna()].cell.nunique(); d = d[d.obj_fba.notna()].copy()
d["zone"] = pd.cut(d[ax], [ed[0] - 1e-9] + ed[1:]); d["und_pct"] = (d.width > tol) * 100
t = d.pivot_table(index="zone", columns="exchange", values="und_pct", aggfunc="mean", observed=True).round(0); t["cells"] = d.groupby("zone", observed=True).cell.nunique()
print(f"{f}: t = {d.t.iloc[0]:g} min, feasible {d.cell.nunique()}, infeasible {inf}; undetermined = range > {tol:g} at fraction {fr}"); pd.set_option("display.width", 250); print(t.to_string())
if len(sys.argv) > 6: t.to_csv(sys.argv[6])
