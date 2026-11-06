"""Free uptake of substrates that BioFVM does not track. usage: free_uptake_check.py RUNDIR [SNAPSHOT=last] [N=40 cells per cell type]
The addon sets the lower bound of the exchanges listed in <transport_model> (the BioFVM densities) from the Michaelis-Menten bound; every other exchange (or demand) reaction keeps the bound of the SBML file, so a negative lower bound there lets the cell take up a substrate from nowhere, which the voxel ledger cannot see.
For living cells of the snapshot the problem is rebuilt as in fva_replay.py and solved with Gurobi twice: as it was (open) and with the lower bound of every unmodelled exchange with a negative lower bound set to 0 (closed). Reported per cell type: the unmodelled exchanges with a negative lower bound; growth closed / open (1 = the free uptake is not used for growth; the model is infeasible when closed if nan);
and for each of those exchanges the median uptake in the open optimum and the largest uptake that the optimum allows (FVA minimum at fraction 0.99999). Writes RUNDIR/../results/free_uptake_<run>.csv relative to the use case (results/)."""
import sys, os, numpy as np, pandas as pd, scipy.io as sio, xml.etree.ElementTree as ET, warnings, cobra
warnings.filterwarnings("ignore")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); from pcxml import labels; from fva_replay import cell_models, rid, snapshots
from cobra.flux_analysis import flux_variability_analysis as fva
d = sys.argv[1].rstrip("/"); snap = sys.argv[2] if len(sys.argv) > 2 else "last"; n = int(sys.argv[3]) if len(sys.argv) > 3 else 40
uc = os.path.dirname(os.path.dirname(os.path.abspath(d))); xf = snapshots(d)[-1] if snap == "last" else f"{d}/output{int(snap):08d}.xml"
root = ET.parse(xf).getroot(); lab = labels(root); m = sio.loadmat(xf.replace(".xml", "_cells.mat"))["cells"]; col = lambda k: m[lab.index(k)]
set_ = ET.parse(f"{d}/PhysiCell_settings.xml").getroot(); models = cell_models(f"{d}/PhysiCell_settings.xml", d)
ex = [k for k in lab if k.startswith("R_EX") and not k.endswith("_lb")]; types = np.asarray(col("cell_type")).astype(int); alive = np.asarray(col("dead")) == 0; rows = []
for i, cd in enumerate(set_.findall(".//cell_definitions/cell_definition")):
    ty = int(cd.get("ID", i))
    if ty not in models: continue
    modeled = {e.find("fba_flux").text.strip() for e in cd.findall(".//transport_model/exchange")}; modeled |= {x[2:] for x in modeled if x.startswith("R_")}
    M = cobra.io.read_sbml_model(models[ty]); M.solver = "gurobi"; free = [r for r in M.reactions if r.boundary and r.lower_bound < 0 and r.id not in modeled and r.id.lower().replace("ex_", "ex_") not in {k[2:] for k in ex}]
    ok = np.flatnonzero(alive & (types == ty))[:n]
    if not len(ok) or not free: print(cd.get("name"), "cells", len(ok), "unmodelled uptake reactions", len(free)); continue
    res = {r.id: [] for r in free}; g = []; gh = []; fmin = {r.id: [] for r in free}
    for j in ok:
        with M as MM:
            for k in ex:
                r = rid(MM, k)
                if r is not None: r.lower_bound = float(col(k + "_lb")[j])
            go = MM.slim_optimize(error_value=np.nan); sol = MM.optimize() if not np.isnan(go) else None
            if sol is not None:
                for r in free: res[r.id].append(sol.fluxes[r.id])
                F = fva(MM, reaction_list=[MM.reactions.get_by_id(r.id) for r in free], fraction_of_optimum=0.99999, processes=1)
                for r in free: fmin[r.id].append(F.loc[r.id, "minimum"])
            lbs = {r.id: r.lower_bound for r in free}
            for r in free: MM.reactions.get_by_id(r.id).lower_bound = 0.0
            gc = MM.slim_optimize(error_value=np.nan)
            for r in free: MM.reactions.get_by_id(r.id).lower_bound = lbs[r.id] if r.id in ("EX_h_e", "EX_h2o_e") else 0.0
            g_B = MM.slim_optimize(error_value=np.nan)            # all closed except H+ and H2O
            for r in free: MM.reactions.get_by_id(r.id).lower_bound = 0.0 if r.id in ("EX_h_e", "EX_h2o_e") else lbs[r.id]
            g_A = MM.slim_optimize(error_value=np.nan)            # only H+ and H2O closed
        g.append((go, gc)); gh.append((g_A, g_B))
    g = np.array(g); ratio = g[:, 1] / g[:, 0]
    print(f"{cd.get('name')} ({os.path.basename(models[ty])}), t = {float(root.find('.//current_time').text):g} min, {len(ok)} cells: unmodelled exchanges with lb < 0: {len(free)}; growth closed / open: median {np.nanmedian(ratio):.4g} (min {np.nanmin(ratio):.4g}, max {np.nanmax(ratio):.4g}), infeasible when closed: {int(np.isnan(g[:, 1]).sum())}")
    gh = np.array(gh); print(f"    only H+ and H2O closed: growth ratio median {np.nanmedian(gh[:, 0] / g[:, 0]):.4g} (infeasible {int(np.isnan(gh[:, 0]).sum())}); everything except H+ and H2O closed: median {np.nanmedian(gh[:, 1] / g[:, 0]):.4g} (infeasible {int(np.isnan(gh[:, 1]).sum())})")
    for r in free:
        if res[r.id]: rows.append(dict(run=os.path.basename(d), t=float(root.find(".//current_time").text), cell_type=cd.get("name"), reaction=r.id, lb=r.lower_bound, median_flux_open=float(np.median(res[r.id])), median_fva_min=float(np.median(fmin[r.id])),
                                       uptake_used_pct=100 * float(np.mean(np.array(res[r.id]) < -1e-6)), growth_ratio_closed_open=float(np.nanmedian(ratio))))
D = pd.DataFrame(rows); os.makedirs(f"{uc}/results", exist_ok=True); out = f"{uc}/results/free_uptake_{os.path.basename(d)}_{os.path.basename(xf)[6:14]}.csv"; D.to_csv(out, index=False)
pd.set_option("display.width", 220); print(D[["cell_type", "reaction", "lb", "median_flux_open", "median_fva_min", "uptake_used_pct"]].round(4).to_string(index=False)); print("wrote", out)
