"""UT-TRIM: does trimming change the behaviour of a model? usage: trim_check.py -> results/trim_check.csv, results/trim_removed.csv
For each pair (full model from dfba_sbml_models, models_full/; trimmed model as shipped to the runs) the full model gets the medium of the trimmed one: every exchange (boundary) reaction of the trimmed model copies its bounds, every other exchange of the full model is closed (both bounds 0); the other reactions keep the bounds of the full file. Then
 - FBA optimum of both (objective of the file);
 - FVA at fraction 0.99999 (gurobi) on the exchanges of the trimmed model, compared limit by limit;
 - FVA (same fraction) of the reactions that the trimming removed, in the full model: how many can carry flux (range > 1e-6) in that medium.
The trimmed M. barkeri that the runs read differs from the one in dfba_sbml_models (sha256 in the table)."""
import os, hashlib, warnings, numpy as np, pandas as pd, cobra
from cobra.flux_analysis import flux_variability_analysis as fva
warnings.filterwarnings("ignore")
UC = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); R = os.path.join(UC, ".."); F = 0.99999
P = [("C. beijerinckii iCB925", "models_full/iCB925.xml", "crossfeeding/config/shipped/sbml_models/C_beijerinckii_model_iCB925_trimmed.xml"),
     ("M. barkeri iMG746", "models_full/M_barkeri_model_iMG746.xml", "crossfeeding/config/shipped/sbml_models/M_barkeri_model_iMG746_trimmed.xml"),
     ("M. barkeri iMG746 (trimmed file of dfba_sbml_models)", "models_full/M_barkeri_model_iMG746.xml", "unit_test/models_full/M_barkeri_model_iMG746_trimmed_repo.xml"),
     ("MCF7 Zielinski 2017 core", "models_full/MCF7_Zielinski_2017_core.sbml", "tumour/config/shipped/MCF7_Zielinski_2017_core_trimmed.sbml")]
sha = lambda p: hashlib.sha256(open(p, "rb").read()).hexdigest()[:12]
rows, rem = [], []
for name, full, trim in P:
    pf, pt = os.path.join(UC, full), os.path.join(R, trim); Mf, Mt = cobra.io.read_sbml_model(pf), cobra.io.read_sbml_model(pt)
    for M in (Mf, Mt): M.solver = "gurobi"
    low = {r.id.upper(): r for r in Mf.reactions}                        # exchange ids differ in case (iCB925: ex_ against EX_)
    mp = {r.id: low[r.id.upper()] for r in Mt.reactions if r.id.upper() in low}   # trimmed reaction -> full reaction
    ex_t = [r.id for r in Mt.reactions if r.boundary and r.id in mp]
    for r in Mf.reactions:                                              # medium of the trimmed model, other exchanges closed
        if r.boundary: r.bounds = (0, 0)
    for i in ex_t: mp[i].bounds = Mt.reactions.get_by_id(i).bounds
    # variant B: the trimmed topology with the internal bounds of the full model (separates the removal of reactions from the tightening of bounds)
    Mb = Mt.copy(); Mb.solver = "gurobi"
    for r in Mb.reactions:
        if not r.boundary and r.id in mp: r.bounds = mp[r.id].bounds
    out = {}
    for tag, M in (("trimmed", Mt), ("trimmed_fullbounds", Mb), ("full", Mf)): out[tag] = M.slim_optimize()
    def fv(M, ids): return fva(M, reaction_list=ids, fraction_of_optimum=F, processes=1)
    ids_f = [mp[i].id for i in ex_t]; fT, fB, fF = fv(Mt, ex_t), fv(Mb, ex_t), fv(Mf, ids_f); fF.index = ex_t
    dd = lambda a, b: pd.concat([(a.minimum - b.minimum).abs(), (a.maximum - b.maximum).abs()], axis=1).max(axis=1)
    dTF, dBF = dd(fT, fF), dd(fB, fF)
    removed = [r.id for r in Mf.reactions if r.id.upper() not in {x.upper() for x in Mt.reactions.list_attr("id")}]
    fr = fva(Mf, reaction_list=removed, fraction_of_optimum=F, processes=1) if removed and out["full"] > 1e-9 else pd.DataFrame(columns=["minimum", "maximum"])
    act = fr[(fr.maximum - fr.minimum).abs() > 1e-6]
    for i in act.index: rem.append(dict(model=name, reaction=i, minimum=fr.loc[i, "minimum"], maximum=fr.loc[i, "maximum"]))
    rows.append(dict(model=name, full_sha=sha(pf), trimmed_sha=sha(pt), objective_full=",".join(r.id for r in Mf.reactions if r.objective_coefficient), objective_trimmed=",".join(r.id for r in Mt.reactions if r.objective_coefficient),
                     reactions_full=len(Mf.reactions), reactions_trimmed=len(Mt.reactions), removed=len(removed), removed_with_range=len(act), only_in_trimmed=len([r for r in Mt.reactions if r.id not in mp]),
                     internal_bounds_changed=sum(1 for r in Mt.reactions if not r.boundary and r.id in mp and r.bounds != mp[r.id].bounds), opt_full=out["full"], opt_trimmed_fullbounds=out["trimmed_fullbounds"], opt_trimmed=out["trimmed"],
                     fva_maxdiff_topology_only=dBF.max(), fva_maxdiff_shipped=dTF.max(), exchanges=len(ex_t)))
D = pd.DataFrame(rows); D.to_csv(f"{UC}/results/trim_check.csv", index=False); pd.DataFrame(rem).to_csv(f"{UC}/results/trim_removed.csv", index=False); pd.set_option("display.width", 250); print(D.T.to_string())
