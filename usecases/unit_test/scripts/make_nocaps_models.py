"""Build the trimmed models without the FVA caps (models_full/*_trimmed_nocaps.xml). usage: make_nocaps_models.py
Starting from the trimmed file the runs read (same reactions and metabolites, exchange bounds of the medium), every internal reaction gets the bounds of the full model (dfba_sbml_models, models_full/), so the trimming only removes blocked reactions and gap metabolites; ATPM lower bound 0 (no maintenance requirement; metabolic death is disabled in the cross-feeding configuration); the upper bound of ATPM is the full model's. Exchange reactions are not touched (they carry the medium), biomass keeps the bounds of the full model (the addon sets the objective upper bound to max_growth_rate). Log: results/nocaps_models.csv."""
import os, hashlib, warnings, pandas as pd, cobra
warnings.filterwarnings("ignore")
UC = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); R = os.path.join(UC, "..")
P = [("models_full/iCB925.xml", "crossfeeding/config/shipped/sbml_models/C_beijerinckii_model_iCB925_trimmed.xml", "models_full/C_beijerinckii_model_iCB925_trimmed_nocaps.xml"),
     ("models_full/M_barkeri_model_iMG746.xml", "crossfeeding/config/shipped/sbml_models/M_barkeri_model_iMG746_trimmed.xml", "models_full/M_barkeri_model_iMG746_trimmed_nocaps.xml")]
sha = lambda p: hashlib.sha256(open(p, "rb").read()).hexdigest()[:12]; rows = []
for full, trim, out in P:
    Mf, Mt = cobra.io.read_sbml_model(os.path.join(UC, full)), cobra.io.read_sbml_model(os.path.join(R, trim)); low = {r.id.upper(): r for r in Mf.reactions}; n = 0
    for r in Mt.reactions:
        f = low.get(r.id.upper())
        if f is None or r.boundary: continue
        if r.bounds != f.bounds: n += 1
        r.bounds = f.bounds
    if "ATPM" in Mt.reactions: Mt.reactions.ATPM.lower_bound = 0.0
    cobra.io.write_sbml_model(Mt, os.path.join(UC, out)); M2 = cobra.io.read_sbml_model(os.path.join(UC, out))
    assert [r.id for r in M2.reactions] == [r.id for r in Mt.reactions] and len(M2.metabolites) == len(Mt.metabolites)
    rows.append(dict(file=out, from_trimmed=trim, sha256=sha(os.path.join(UC, out)), reactions=len(M2.reactions), metabolites=len(M2.metabolites), internal_bounds_reset=n, atpm_bounds=str(M2.reactions.ATPM.bounds), objective=",".join(r.id for r in M2.reactions if r.objective_coefficient)))
D = pd.DataFrame(rows); D.to_csv(f"{UC}/results/nocaps_models.csv", index=False); print(D.T.to_string())
