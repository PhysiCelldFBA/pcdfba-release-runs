"""Do the internal bounds of the shipped trimmed M. barkeri iMG746 bind in the cross-feeding runs? usage: mb_caps.py RUNDIR [SNAPSHOTS=10,30,60] [N=200 cells]
For living M. barkeri cells of the snapshots the problem of the cell is rebuilt as in common/fva_replay.py (recorded lower bounds of the recorded exchanges, other bounds from the shipped trimmed file) and solved (gurobi) with three sets of internal bounds, everything else identical: the shipped trimmed model; 'relaxed' (each changed internal bound widened to the wider of the shipped and the full iMG746 value, so never tighter than shipped; this removes the caps); 'full_bounds' (the internal bounds of the full iMG746, unit_test/models_full, which also fixes ATPM at 2 and is tighter than shipped in places).
Reports: growth with and without the caps (and against the recorded growth_rate); the number of cells in which the optimum changes by more than 1e-6; the caps that are at their bound in the shipped optimum (flux = bound, bound differs from the full model). Writes results/mb_caps_<run>.csv (one row per cell)."""
import sys, os, numpy as np, pandas as pd, scipy.io as sio, xml.etree.ElementTree as ET, warnings, cobra
warnings.filterwarnings("ignore")
UC = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); sys.path.insert(0, os.path.join(UC, "..", "..", "common")); from pcxml import labels; from fva_replay import cell_models, rid
d = sys.argv[1].rstrip("/"); snaps = [int(x) for x in (sys.argv[2] if len(sys.argv) > 2 else "10,30,60").split(",")]; n = int(sys.argv[3]) if len(sys.argv) > 3 else 200
models = cell_models(f"{d}/PhysiCell_settings.xml", d); ty = [t for t, p in models.items() if "barkeri" in p][0]
Mt = cobra.io.read_sbml_model(models[ty]); Mf = cobra.io.read_sbml_model(os.path.join(UC, "..", "unit_test", "models_full", "M_barkeri_model_iMG746.xml")); Mt.solver = "gurobi"
Mb = Mt.copy(); Mb.solver = "gurobi"; chg = {}
for r in Mb.reactions:
    if not r.boundary and r.id in Mf.reactions:
        f = Mf.reactions.get_by_id(r.id)
        if r.bounds != f.bounds: chg[r.id] = (f.bounds, r.bounds); r.bounds = f.bounds
Mr = Mt.copy(); Mr.solver = "gurobi"                                   # caps removed: every changed internal bound is widened to the wider of the shipped and the full value (never tighter than the shipped model)
for i in chg:
    r, f = Mr.reactions.get_by_id(i), Mf.reactions.get_by_id(i); r.bounds = (min(r.lower_bound, f.lower_bound), max(r.upper_bound, f.upper_bound))
rows = []; rng = np.random.default_rng(1)
for s in snaps:
    xf = f"{d}/output{s:08d}.xml"; root = ET.parse(xf).getroot(); lab = labels(root); t = float(root.find(".//current_time").text); m = sio.loadmat(xf.replace(".xml", "_cells.mat"))["cells"]; col = lambda k: m[lab.index(k)]
    ex = [k for k in lab if k.startswith("R_EX") and not k.endswith("_lb")]; ok = np.flatnonzero((np.asarray(col("cell_type")).astype(int) == ty) & (np.asarray(col("dead")) == 0)); ok = np.sort(rng.choice(ok, min(n, len(ok)), replace=False))
    for j in ok:
        res = {}
        for tag, M in (("shipped", Mt), ("relaxed", Mr), ("full_bounds", Mb)):
            with M as MM:
                for k in ex:
                    r = rid(MM, k)
                    if r is not None: r.lower_bound = float(col(k + "_lb")[j])
                v = MM.slim_optimize(error_value=np.nan); res[tag] = v
                if tag == "shipped" and not np.isnan(v):
                    fl = MM.optimize().fluxes; atb = [i for i in chg if abs(fl[i]) > 1e-9 and (abs(fl[i] - MM.reactions.get_by_id(i).upper_bound) < 1e-7 or abs(fl[i] - MM.reactions.get_by_id(i).lower_bound) < 1e-7)]
                    res["n_at_cap"] = len(atb); res["caps"] = ";".join(sorted(atb))
        rows.append(dict(snap=s, t_min=t, cell=int(col("ID")[j]), growth_rec=float(col("growth_rate")[j]), growth_shipped=res["shipped"], growth_relaxed=res["relaxed"], growth_full_bounds=res["full_bounds"], n_at_cap=res.get("n_at_cap"), caps=res.get("caps")))
D = pd.DataFrame(rows); D["ratio"] = D.growth_relaxed / D.growth_shipped; D["changed"] = (D.growth_relaxed - D.growth_shipped).abs() > 1e-6; D["ratio_full"] = D.growth_full_bounds / D.growth_shipped
out = f"{UC}/results/mb_caps_{os.path.basename(d)}.csv"; D.to_csv(out, index=False)
print(f"capped reactions (bounds differ from the full model): {len(chg)}; MB cells {len(D)}; replay growth against recorded: max rel diff {((D.growth_shipped - D.growth_rec).abs() / D.growth_rec.abs().replace(0, np.nan)).max():.2g}")
for s, g in D.groupby("snap"): print(f"t = {g.t_min.iloc[0]:g} min: {len(g)} cells; growth changes by > 1e-6 when the caps are relaxed in {int(g.changed.sum())}; ratio relaxed/shipped median {g.ratio.median():.4g} (min {g.ratio.min():.4g}, max {g.ratio.max():.4g}); full-model bounds (also ATPM 2 fixed) median {g.ratio_full.median():.4g}; cells with >= 1 cap at its bound {int((g.n_at_cap > 0).sum())}, median number {g.n_at_cap.median():g}")
c = pd.Series(";".join(D.caps.dropna()).split(";")).value_counts(); print("caps at their bound, by number of cell-times:"); print(c.head(12).to_string())
