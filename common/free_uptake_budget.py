"""Budget of the uptake of substrates that BioFVM does not track. usage: free_uptake_budget.py RUNDIR [N=200 cells per type and snapshot] [EVERY=1 snapshot step] [SOLVENT=EX_h_e,EX_h2o_e]
For every exchange or demand reaction with a negative lower bound that is not in <transport_model> (the BioFVM densities), over the snapshots of a run: replay (as fva_replay.py, Gurobi) a random sample of living cells of each type; per reaction
 - the uptake of the open optimum (FBA) and the largest uptake the optimum allows (FVA minimum, fraction 0.99999);
 - the amount drawn by the population: sum over the snapshots (trapezoid in time) of N_living x mean(uptake x dry mass of the cell), mmol, with the dry mass of the addon (volume x (1 - fluid fraction) x cell_density), and the same divided by the volume of the domain (mM): the mean concentration that the domain would have to hold to supply it;
 - 'unused' when no sampled cell could take it up at any snapshot (largest allowed uptake < 1e-6 mmol/gDW/h); closing every unused reaction (except SOLVENT) is tested: growth closed / open in each sampled cell; mmol_per_gDW_formed = uptake / growth rate, the amount of the substrate taken up per gram dry weight of biomass formed (independent of the domain volume).
Writes results/free_uptake_budget_<run>.csv and prints the table. The sample gives the population mean; the statistical error is the standard error of the mean flux x mass over the sampled cells (column se_pct)."""
import sys, os, numpy as np, pandas as pd, scipy.io as sio, xml.etree.ElementTree as ET, warnings, cobra
warnings.filterwarnings("ignore")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); from pcxml import labels; from fva_replay import cell_models, rid, snapshots
from cobra.flux_analysis import flux_variability_analysis as fva
d = sys.argv[1].rstrip("/"); n = int(sys.argv[2]) if len(sys.argv) > 2 else 200; every = int(sys.argv[3]) if len(sys.argv) > 3 else 1
solvent = set((sys.argv[4] if len(sys.argv) > 4 else "EX_h_e,EX_h2o_e").split(",")); uc = os.path.dirname(os.path.dirname(os.path.abspath(d)))
S = ET.parse(f"{d}/PhysiCell_settings.xml").getroot(); models = cell_models(f"{d}/PhysiCell_settings.xml", d); g = lambda a: float(S.find(f".//domain/{a}").text)
Vdom = (g("x_max") - g("x_min")) * (g("y_max") - g("y_min")) * (g("z_max") - g("z_min")) * 1e-15      # litres
if g("z_max") - g("z_min") < 1e-9: Vdom = (g("x_max") - g("x_min")) * (g("y_max") - g("y_min")) * g("dz") * 1e-15
CD = {}
for i, cd in enumerate(S.findall(".//cell_definitions/cell_definition")):
    ty = int(cd.get("ID", i))
    if ty not in models: continue
    M = cobra.io.read_sbml_model(models[ty]); M.solver = "gurobi"; modeled = {e.find("fba_flux").text.strip() for e in cd.findall(".//transport_model/exchange")}; modeled |= {x[2:] for x in modeled if x.startswith("R_")}
    CD[ty] = dict(name=cd.get("name"), M=M, dens=float(cd.find(".//growth_model/cell_density").text), free=[r.id for r in M.reactions if r.boundary and r.lower_bound < 0 and r.id not in modeled])
xs = snapshots(d)[::every]; rng = np.random.default_rng(1); rec = []; T = []
for xf in xs:
    root = ET.parse(xf).getroot(); lab = labels(root); m = sio.loadmat(xf.replace(".xml", "_cells.mat"))["cells"]; col = lambda k: m[lab.index(k)]; t = float(root.find(".//current_time").text); T.append(t)
    ex = [k for k in lab if k.startswith("R_EX") and not k.endswith("_lb")]; types = np.asarray(col("cell_type")).astype(int); alive = np.asarray(col("dead")) == 0
    for ty, c in CD.items():
        ix = np.flatnonzero(alive & (types == ty))
        if not len(ix) or not c["free"]: continue
        ok = np.sort(rng.choice(ix, min(n, len(ix)), replace=False)); M = c["M"]
        for j in ok:
            mass = float(col("total_volume")[j]) * (1 - float(col("fluid_fraction")[j])) * c["dens"] * 1e-12            # g
            with M as MM:
                for k in ex:
                    r = rid(MM, k)
                    if r is not None: r.lower_bound = float(col(k + "_lb")[j])
                go = MM.slim_optimize(error_value=np.nan)
                if np.isnan(go): continue
                fl = MM.optimize().fluxes; F = fva(MM, reaction_list=[MM.reactions.get_by_id(i) for i in c["free"]], fraction_of_optimum=0.99999, processes=1)
                unused = [i for i in c["free"] if i not in solvent and F.loc[i, "minimum"] > -1e-6]
                for i in unused: MM.reactions.get_by_id(i).lower_bound = 0.0
                gc = MM.slim_optimize(error_value=np.nan)
            for i in c["free"]: rec.append(dict(t=t, type=c["name"], ty=ty, n_alive=len(ix), n_samp=len(ok), cell=int(col("ID")[j]), mass=mass, reaction=i, flux=float(fl[i]), fva_min=float(F.loc[i, "minimum"]), growth=go, growth_closed_unused=gc))
D = pd.DataFrame(rec); out = []
for (nm, rx), q in D.groupby(["type", "reaction"]):
    per_t = q.assign(up=np.maximum(-q.flux, 0) * q.mass, net=-q.flux * q.mass).groupby("t").agg(up=("up", "mean"), net=("net", "mean"), n=("n_alive", "first"), s=("up", "sem")).reset_index()
    tot = np.trapz(per_t.up * per_t.n, per_t.t / 60.0) if len(per_t) > 1 else np.nan; net = np.trapz(per_t.net * per_t.n, per_t.t / 60.0) if len(per_t) > 1 else np.nan
    se = np.sqrt(np.trapz((per_t.s * per_t.n) ** 2, per_t.t / 60.0)) if len(per_t) > 1 else np.nan
    out.append(dict(run=os.path.basename(d), type=nm, reaction=rx, lb=CD[q.ty.iloc[0]]["M"].reactions.get_by_id(rx).lower_bound, cells_sampled=len(q), max_uptake_allowed=-q.fva_min.min(), cells_using_pct=100 * (q.flux < -1e-6).mean(),
                    mean_uptake_when_used=-q.flux[q.flux < -1e-6].mean() if (q.flux < -1e-6).any() else 0.0, drawn_mmol=tot, net_mmol=net, se_pct=100 * se / tot if tot and tot > 0 else np.nan, domain_L=Vdom, mean_conc_mM=tot / Vdom if Vdom else np.nan,
                    mmol_per_gDW_formed=float((np.maximum(-q.flux, 0) * q.mass)[q.growth > 1e-6].sum() / (q.growth * q.mass)[q.growth > 1e-6].sum()) if (q.growth > 1e-6).any() else np.nan,
                    unused=bool(q.fva_min.min() > -1e-6), solvent=rx in solvent))
O = pd.DataFrame(out); gr = D.drop_duplicates(["t", "cell", "type"]); gr = gr[gr.growth > 1e-4]; ratio = (gr.growth_closed_unused / gr.growth).replace([np.inf, -np.inf], np.nan)
os.makedirs(f"{uc}/results", exist_ok=True); fn = f"{uc}/results/free_uptake_budget_{os.path.basename(d)}.csv"; O.to_csv(fn, index=False); D.to_csv(fn.replace(".csv", "_cells.csv"), index=False)
pd.set_option("display.width", 250); pd.set_option("display.float_format", lambda x: f"{x:.4g}")
print(f"{os.path.basename(d)}: {len(xs)} snapshots, {len(gr)} cell-times replayed, domain {Vdom:.3g} L; growth with every unused exchange closed / open (cell-times with growth > 1e-4 /h): min {np.nanmin(ratio):.9f}, max {np.nanmax(ratio):.9f}")
print(O.drop(columns=["run", "domain_L", "cells_sampled"]).to_string(index=False)); print("wrote", fn)
