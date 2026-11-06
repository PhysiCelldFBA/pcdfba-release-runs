"""UT-ELEM: elemental balance of the metabolic models used (cobra check_mass_balance, charge included). usage: elem_balance.py   -> results/elem_balance.csv (one row per model), results/elem_unbalanced.csv (the unbalanced internal reactions).
Internal = not an exchange/demand/sink reaction (boundary). A reaction is checkable when every metabolite has a chemical formula; the biomass reaction is reported separately. Imbalance is in elements (and charge); H and charge are given as in the model."""
import os, hashlib, warnings, pandas as pd, cobra
warnings.filterwarnings("ignore")
UC = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); R = os.path.join(UC, "..")
M = [("e_coli_core (unit test, acetate, colony, motility)", "unit_test/config/shipped/Ecoli_core.xml"), ("e_coli_core ATPM 8.39", "unit_test/models/Ecoli_core_ATPM8.39.xml"), ("iML1515", "unit_test/models/iML1515.xml"),
     ("C. beijerinckii iCB925 (trimmed)", "crossfeeding/config/shipped/sbml_models/C_beijerinckii_model_iCB925_trimmed.xml"), ("M. barkeri iMG746 (trimmed)", "crossfeeding/config/shipped/sbml_models/M_barkeri_model_iMG746_trimmed.xml"),
     ("MCF7 Zielinski 2017 core (trimmed)", "tumour/config/shipped/MCF7_Zielinski_2017_core_trimmed.sbml")]
rows, bad = [], []
for name, rel in M:
    p = os.path.join(R, rel); m = cobra.io.read_sbml_model(p); bio = [r.id for r in m.reactions if r.objective_coefficient != 0]
    bio = sorted(set(bio) | {r.id for r in m.reactions if "biomass" in r.id.lower()}); internal = [r for r in m.reactions if not r.boundary and r.id not in bio]; chk = unb = nof = 0
    for r in internal:
        if any(not x.formula for x in r.metabolites): nof += 1; continue
        chk += 1; imb = r.check_mass_balance()
        if imb: unb += 1; bad.append(dict(model=name, reaction=r.id, imbalance=str(imb), n_metabolites=len(r.metabolites)))
    rows.append(dict(model=name, file=rel, sha256=hashlib.sha256(open(p, "rb").read()).hexdigest()[:12], reactions=len(m.reactions), metabolites=len(m.metabolites), exchange=sum(r.boundary for r in m.reactions), internal=len(internal), checkable=chk, balanced=chk - unb,
                     unbalanced=unb, no_formula=nof, biomass_reactions=",".join(bio), metabolites_without_formula=sum(not x.formula for x in m.metabolites)))
D = pd.DataFrame(rows); os.makedirs(f"{UC}/results", exist_ok=True); D.to_csv(f"{UC}/results/elem_balance.csv", index=False); pd.DataFrame(bad).to_csv(f"{UC}/results/elem_unbalanced.csv", index=False)
pd.set_option("display.width", 250); print(D.drop(columns=["file"]).to_string(index=False))
