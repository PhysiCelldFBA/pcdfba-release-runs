
import os, sys, time, json, resource
import cobra
import cometspy as cm
import pandas as pd
import numpy as np

t0 = time.time()

COMETS_HOME = os.environ["COMETS_HOME"]

# ---- load model, identify ids robustly ----
model = cobra.io.read_sbml_model("Ecoli_core.xml")
rxn_ids = [r.id for r in model.reactions]
def find_id(*cands):
    for c_ in cands:
        if c_ in rxn_ids:
            return c_
    raise KeyError(cands)

biomass_id = find_id("R_Biomass_Ecoli_core", "Biomass_Ecoli_core")
ex_glc = find_id("R_EX_glc__D_e", "EX_glc__D_e")
ex_o2  = find_id("R_EX_o2_e", "EX_o2_e")
ex_ac  = find_id("R_EX_ac_e", "EX_ac_e")
ex_co2 = find_id("R_EX_co2_e", "EX_co2_e")
print("ids:", biomass_id, ex_glc, ex_o2, ex_ac, ex_co2)

model.objective = biomass_id

cmodel = cm.model(model)
cmodel.id = "ecoli"
cmodel.optimizer = "GLOP"

# open uptake (pDFBA re-opens every step; COMETS applies continuous MM cap)
cmodel.change_bounds(ex_glc, -1000., 1000.)
cmodel.change_bounds(ex_o2,  -1000., 1000.)
cmodel.change_bounds(ex_ac,  -1000., 1000.)
cmodel.change_bounds(ex_co2, -1000., 1000.)

# Michaelis-Menten kinetics, exact pDFBA config values.
# CORRECTED: cometspy's change_km() expects Km in mol/L (M), not mmol/L (mM) --
# unlike change_vmax(), which is correctly in mmol/gDW/h. pDFBA's config values
# (0.02, 0.002, 1.0) are in mM; divide by 1000 to pass the intended M value.
cmodel.change_vmax(ex_glc, 8.0);  cmodel.change_km(ex_glc, 0.02/1000.0)
cmodel.change_vmax(ex_o2, 18.0);  cmodel.change_km(ex_o2, 0.002/1000.0)
cmodel.change_vmax(ex_ac, 5.0);   cmodel.change_km(ex_ac, 1.0/1000.0)

# growth-rate cap (pDFBA max_growth_rate = 0.86 /h)
cmodel.change_bounds(biomass_id, 0.0, 0.86)

# initial biomass (v2, 2026-11-18): the dry biomass of the 20 PhysiCelldFBA founders, 0.1449 g/L x 4.6656e-11 L = 6.76e-12 g
# (0.338 pg per cell = volume x (1 - fluid fraction) x density). v1 used 1.352 pg per cell (fluid fraction left out), 4x too much.
init_biomass_g = 0.1448902606310014 * 4.6656e-11
cmodel.initial_pop = [[0, 0, init_biomass_g]]

# ---- layout: single well-mixed box = pDFBA's 36x36x36 um closed box ----
layout = cm.layout([cmodel])
layout.grid = [1, 1]
V_L = 4.6656e-11  # 36x36x36 um^3 box volume in liters

def mmol(conc_mM):
    return conc_mM * V_L

# derive metabolite ids from the exchange reactions themselves, not endswith()
# matching on bare metabolite ids (endswith("o2_e") wrongly matches "co2_e" too)
met_glc = [met.id for met in model.reactions.get_by_id(ex_glc).metabolites][0]
met_o2  = [met.id for met in model.reactions.get_by_id(ex_o2).metabolites][0]
met_ac  = [met.id for met in model.reactions.get_by_id(ex_ac).metabolites][0]
met_co2 = [met.id for met in model.reactions.get_by_id(ex_co2).metabolites][0]
print("met ids:", met_glc, met_o2, met_ac, met_co2)

layout.set_specific_metabolite(met_glc, mmol(8.0))
layout.set_specific_metabolite(met_o2,  mmol(5.0))
layout.set_specific_metabolite(met_ac,  mmol(0.0))
layout.set_specific_metabolite(met_co2, mmol(0.0))

# also open all other exchanged metabolites to a large, non-limiting static value
# (non-curated-exchange starvation fix, as established for this scenario)
all_ex_mets = set()
for r in model.reactions:
    if r.id.startswith("R_EX_") or r.id.startswith("EX_"):
        for met_obj in r.metabolites:
            all_ex_mets.add(met_obj.id)
curated = {met_glc, met_o2, met_ac, met_co2}
for met_id in all_ex_mets - curated:
    try:
        layout.set_specific_metabolite(met_id, mmol(1000.0))
    except Exception as e:
        print("skip", met_id, e)

params = cm.params()
params.set_param("timeStep", 1.0/60.0)   # 1 min, matches pDFBA save interval granularity
params.set_param("spaceWidth", 36e-4)    # 36 um box edge (cm) -- diffusion irrelevant (1 box)
params.set_param("minSpaceBiomass", 1e-15)   # v2: COMETS default (~2.5e-11 g) would delete the 6.8e-12 g founders
params.set_param("maxCycles", 210)       # leg 1: 0 -> 210 min
params.set_param("writeTotalBiomassLog", True)
params.set_param("writeMediaLog", True)
params.set_param("MediaLogRate", 1)
params.set_param("totalBiomassLogRate", 1)
params.set_param("maxSpaceBiomass", 10)  # large cap, not meant to bind


def _fix_commas(path):
    with open(path) as f:
        txt = f.read()
    if ',' in txt:
        txt2 = txt.replace(',', '.')
        with open(path, 'w') as f:
            f.write(txt2)

def run_leg(layout_obj, params_obj, tag, wd="./"):
    suffix = "_" + tag
    layout_obj.write_necessary_files(wd, suffix)
    c_global = wd + ".current_global" + suffix
    c_package = wd + ".current_package" + suffix
    c_script = wd + ".current_script" + suffix
    params_obj.write_params(c_global, c_package)
    with open(c_script, "w") as f:
        f.write("load_comets_parameters .current_global" + suffix + "\n")
        f.write("load_package_parameters .current_package" + suffix + "\n")
        f.write("load_layout .current_layout" + suffix)
    cmd = f'bash -c "cd {wd} && {COMETS_HOME}/bin/comets_scr_fixed {c_script}"'
    import subprocess as sp
    p = sp.run(cmd, shell=True, capture_output=True, text=True)
    print("RUN", tag, "rc=", p.returncode)
    if p.returncode != 0:
        print("STDOUT TAIL:", p.stdout[-3000:])
        print("STDERR TAIL:", p.stderr[-3000:])
        raise RuntimeError(f"COMETS leg {tag} failed")
    tbm_path = wd + params_obj.all_params["TotalBiomassLogName"]
    media_path = wd + params_obj.all_params["MediaLogName"]
    _fix_commas(tbm_path)
    _fix_commas(media_path)
    tbm = pd.read_csv(tbm_path, sep=r"\s+", header=None, names=["cycle"] + layout_obj.get_model_ids())
    media = pd.read_csv(media_path, sep=r"\s+", names=["metabolite","cycle","x","y","conc_mmol"])
    return tbm, media

tbm1, media1 = run_leg(layout, params, "leg1")
print("leg1 done, cycles:", tbm1.cycle.max())

# ---- extract final state at end of leg1, build leg2 ----
last_cycle = int(tbm1.cycle.max())
final_biomass_g = float(tbm1.loc[tbm1.cycle==last_cycle, "ecoli"].values[0])
m1_last = media1[media1.cycle==last_cycle]
def get_conc_mmol(met):
    row = m1_last[m1_last.metabolite==met]
    return float(row.conc_mmol.values[0]) if len(row) else 0.0

glc_mmol_1 = get_conc_mmol(met_glc)
o2_mmol_1  = get_conc_mmol(met_o2)
ac_mmol_1  = get_conc_mmol(met_ac)
co2_mmol_1 = get_conc_mmol(met_co2)
print("end of leg1:", "biomass_g=",final_biomass_g, "glc=",glc_mmol_1,"o2=",o2_mmol_1,"ac=",ac_mmol_1,"co2=",co2_mmol_1)

# leg 2: oxygen force-set to exactly 5 mM (inject_density semantics), others carried over
cmodel2 = cm.model(model)
cmodel2.id = "ecoli"
cmodel2.optimizer = "GLOP"
cmodel2.change_bounds(ex_glc, -1000., 1000.)
cmodel2.change_bounds(ex_o2,  -1000., 1000.)
cmodel2.change_bounds(ex_ac,  -1000., 1000.)
cmodel2.change_bounds(ex_co2, -1000., 1000.)
cmodel2.change_vmax(ex_glc, 8.0);  cmodel2.change_km(ex_glc, 0.02/1000.0)
cmodel2.change_vmax(ex_o2, 18.0);  cmodel2.change_km(ex_o2, 0.002/1000.0)
cmodel2.change_vmax(ex_ac, 5.0);   cmodel2.change_km(ex_ac, 1.0/1000.0)
cmodel2.change_bounds(biomass_id, 0.0, 0.86)
cmodel2.initial_pop = [[0, 0, final_biomass_g]]

layout2 = cm.layout([cmodel2])
layout2.grid = [1, 1]
layout2.set_specific_metabolite(met_glc, glc_mmol_1)
layout2.set_specific_metabolite(met_o2,  mmol(5.0))   # <-- the pulse: force-set to 5 mM
layout2.set_specific_metabolite(met_ac,  ac_mmol_1)
layout2.set_specific_metabolite(met_co2, co2_mmol_1)
for met_id in all_ex_mets - curated:
    try:
        layout2.set_specific_metabolite(met_id, mmol(1000.0))
    except Exception as e:
        print("skip2", met_id, e)

params2 = cm.params()
params2.set_param("timeStep", 1.0/60.0)
params2.set_param("spaceWidth", 36e-4)
params2.set_param("minSpaceBiomass", 1e-15)   # v2: COMETS default (~2.5e-11 g) would delete the 6.8e-12 g founders
params2.set_param("maxCycles", 390)   # leg 2: 210 -> 600 min
params2.set_param("writeTotalBiomassLog", True)
params2.set_param("writeMediaLog", True)
params2.set_param("MediaLogRate", 1)
params2.set_param("totalBiomassLogRate", 1)
params2.set_param("maxSpaceBiomass", 10)

tbm2, media2 = run_leg(layout2, params2, "leg2")
print("leg2 done, cycles:", tbm2.cycle.max())

# ---- assemble full time series ----
tbm1 = tbm1.rename(columns={"ecoli":"biomass_g"})
tbm2 = tbm2.rename(columns={"ecoli":"biomass_g"})
tbm1["time_min"] = tbm1.cycle * 1.0
tbm2["time_min"] = 210 + tbm2.cycle * 1.0

def pivot_media(media, offset):
    p = media.pivot_table(index="cycle", columns="metabolite", values="conc_mmol", aggfunc="first").reset_index()
    p["time_min"] = offset + p.cycle * 1.0
    return p

p1 = pivot_media(media1, 0)
p2 = pivot_media(media2, 210)

full_bio = pd.concat([tbm1[["time_min","biomass_g"]], tbm2[["time_min","biomass_g"]]], ignore_index=True)
full_media = pd.concat([p1, p2], ignore_index=True)
full = pd.merge(full_bio, full_media, on="time_min", how="outer").sort_values("time_min").reset_index(drop=True)

# convert mmol -> mM using box volume
for met, newcol in [(met_glc,"glucose_mM"),(met_o2,"oxygen_mM"),(met_ac,"acetate_mM"),(met_co2,"CO2_mM")]:
    if met in full.columns:
        full[newcol] = full[met] / V_L
full["time_h"] = full.time_min/60.0

full.to_csv("comets_aceref_v2_timeseries.csv", index=False)
print(full.tail(10).to_string())

elapsed = time.time()-t0
peak_mb = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024.0
timing = {"elapsed_s": elapsed, "peak_mem_mb": peak_mb}
with open("comets_aceref_v2_timing.json","w") as f:
    json.dump(timing, f)
print("TIMING", timing)
