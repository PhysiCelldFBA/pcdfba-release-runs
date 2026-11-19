# CF-REF mapped into COMETS 2.12.5 at the tool's own scale (2026-11-19). Based on comets_cfref_48x48_v2.py.
# Same biology as PhysiCelldFBA CF-REF (models, kinetics, biomass caps, medium, Dirichlet edge with H2 / acetate = 0,
# 120 x 120 um, 30 h, initial dry mass and positions), but with boxes that hold a small population (GRID_UM, default 10 um)
# instead of the literal 2.5 um reconstruction, and with COMETS's own biomass motion:
#   MOTION=none        no biomass motion (flowDiffRate 0)
#   MOTION=convection  "Convection 2D" (COMETS's colony-expansion model). Per COMETS 2.12.5 bytecode
#                      (Utility.pressure2D), p = (E/xi) (1 - packedDensity dx^2 / rho)^1.5 for rho > packedDensity dx^2,
#                      rho = TOTAL biomass of all species in the box (g), dx = spaceWidth (cm); each species moves with its
#                      share b_i / rho of d(b)/dt = laplacian(p) (Utility.advection2D). packedDensity is in g/cm2 and pushing
#                      starts at packedDensity dx^2 g per box. p, and so E/xi, is in g cm2/s: E/xi is NOT scale-free.
#     packedDensity: dry mass per area of a packed monolayer of the PhysiCelldFBA cells (CF-REF_s50.xml): volume V
#       (C.b. 1.3, M.b. 1.8 um3), dry mass V (1 - 0.75) 1.04e-12 g/um3, sphere footprint pi r^2 at hexagonal disc packing
#       0.9069 -> C.b. 2.13e-5, M.b. 2.37e-5 g/cm2 (about 2e-11 g per 10 um box).
#     elasticModulus: the COMETS petri-dish protocol example (cometspy docs, petri_dish.md: E 1e-10, xi 1, packedDensity 0.022,
#       spaceWidth 0.02 cm, about 1e-6 g per box) used E/xi = 1e-10 directly; on 10 um boxes holding about 1e-11 g that
#       pushes 1e5 times harder and the explicit update blew up at the onset of pushing (6-6.5 h, run of 2026-11-19).
#       We keep instead the example's dimensionless pushing rate k = (E/xi) / (dx^2 packedDensity dx^2) = 0.0284 /s
#       (K_PUSH), i.e. E = k dx^2 packedDensity dx^2 (C.b. 6.0e-19, M.b. 6.7e-19 at 10 um); k is about 600 x the C.b.
#       growth rate (quasi-static packing) and k dt = 0.43 at dt 0.25 min (explicit limit about 0.58).
#     frictionConstant 1.0, convDiffConstant 0, as in the example. Not tuned against PhysiCelldFBA.
# Grid: GRID_UM boxes over the 120 um domain plus a 1-box static ring OUTSIDE it (Dirichlet: glucose, CO2, NH4, Pi, SO4 at
# their initial values; H2, acetate, methane at 0). Initial biomass: cf_ref_box_biomass.json (5 um boxes) x BIOSCALE,
# each old box summed into the new box containing it. Box volume = spaceWidth^3 (COMETS 2D): a 10 um box is 8x deeper
# than the 1.25 um PhysiCelldFBA slab, so the same biomass is 8x less concentrated per volume.
# SCENARIO=cfref (default): the reduced CF-REF of the review series (120 um, 30 h, H2 / acetate start at 0, Km(H2) of
#   M. barkeri 0.05 mM, trimmed models, initial biomass from cf_ref_box_biomass.json).
# SCENARIO=fig6: the consortium of manuscript Fig. 6 (run CF3-NC50_s50): 200 um box, 50 h, H2 and acetate start at 1 mM
#   (boundary 0 mM), Km(H2) of M. barkeri 0.01 mM, models without the FVA caps (*_trimmed_nocaps.xml), initial biomass
#   from the PhysiCelldFBA cell file (CELLS, default cells_CF3-REF_s50.csv: the 111 cells of the run), one cell's dry mass
#   V (1 - 0.75) 1.04e-12 g/um3 (C.b. V 1.3, M.b. 1.8 um3) put in the box that contains it.
# Logs: total biomass, per-box biomass and media every LOG_MIN (default 30) min; fluxes every FLUX_H (default 10) h.

import os, time, json, resource, glob, subprocess as sp
import cobra
import cometspy as cm
import pandas as pd

t0 = time.time()
GRID_UM = float(os.environ.get("GRID_UM", "10"))
MOTION = os.environ.get("MOTION", "none")
TIMESTEP_MIN = float(os.environ.get("TIMESTEP_MIN", "0.25"))
SCENARIO = os.environ.get("SCENARIO", "cfref"); assert SCENARIO in ("cfref", "fig6"), SCENARIO
T_END_H = float(os.environ.get("T_END_H", "30" if SCENARIO == "cfref" else "50"))
LOG_MIN = float(os.environ.get("LOG_MIN", "30"))
FLUX_H = float(os.environ.get("FLUX_H", "10"))
NTHREADS = int(os.environ.get("NTHREADS", "1"))
MINSPACE = float(os.environ.get("MINSPACE", "1e-25"))
BIOSCALE = float(os.environ.get("BIOSCALE", "0.25"))
DEFVMAX = float(os.environ.get("DEFVMAX", "1000"))
DEFKM = float(os.environ.get("DEFKM", "1e-9"))
TAG = os.environ.get("OUT_TAG", f"10um_{MOTION}")
PACKED = {"c_beijerinckii": 2.13e-5, "m_barkeri": 2.37e-5}   # g/cm2, see header
K_PUSH = float(os.environ.get("K_PUSH", "0.0284"))           # 1/s, dimensionless pushing rate of the COMETS example
FRICTION, CONVDIFF = 1.0, 0.0
_dx = GRID_UM * 1e-4
ELASTIC = {k: K_PUSH * _dx ** 2 * v * _dx ** 2 * FRICTION for k, v in PACKED.items()}   # g cm2/s
assert MOTION in ("none", "convection"), MOTION
print("GRID_UM", GRID_UM, "MOTION", MOTION, "TIMESTEP_MIN", TIMESTEP_MIN, "T_END_H", T_END_H, "LOG_MIN", LOG_MIN,
      "FLUX_H", FLUX_H, "SCENARIO", SCENARIO, "NTHREADS", NTHREADS, "MINSPACE", MINSPACE, "BIOSCALE", BIOSCALE, "TAG", TAG,
      "K_PUSH", K_PUSH, "ELASTIC", ELASTIC)
COMETS_HOME = os.environ["COMETS_HOME"]

SFX = "_nocaps" if SCENARIO == "fig6" else ""
cb_model = cobra.io.read_sbml_model(f"C_beijerinckii_model_iCB925_trimmed{SFX}.xml")
mb_model = cobra.io.read_sbml_model(f"M_barkeri_model_iMG746_trimmed{SFX}.xml")


def find_id(model, *cands):
    ids = [r.id for r in model.reactions]
    for c_ in cands:
        if c_ in ids:
            return c_
    raise KeyError(cands)


cb_bio = find_id(cb_model, "R_biomass", "biomass")
mb_bio = find_id(mb_model, "R_Mb_biomass_65", "Mb_biomass_65")
cb_model.objective = cb_bio
mb_model.objective = mb_bio
cb_ex = {s: find_id(cb_model, f"R_EX_{s}_e", f"EX_{s}_e") for s in ["glc_D", "h2", "nh4", "pi", "so4", "ac", "co2"]}
mb_ex = {s: find_id(mb_model, f"R_EX_{s}_e", f"EX_{s}_e") for s in ["co2", "h2", "ch4", "nh4", "pi", "ac"]}

# kinetics (Km in M, Vmax in mmol/gDW/h), unchanged from CF-REF
cb_kinetics = {"glc_D": (1e-3, 7.5), "nh4": (1e-6, 15.0), "pi": (1e-6, 15.0), "so4": (1e-6, 15.0)}
mb_kinetics = {"co2": (1e-7, 100.0), "h2": (1e-5 if SCENARIO == "fig6" else 5e-5, 20.0), "nh4": (1e-7, 10.0), "pi": (1e-7, 10.0), "ac": (1e-7, 3.0)}


def make(model, mid, kin, ex, bio, cap):
    m = cm.model(model)
    m.id = mid
    m.optimizer = "GLOP"
    for sub, (km, vmax) in kin.items():
        m.change_bounds(ex[sub], -1000., 1000.)
        m.change_vmax(ex[sub], vmax); m.change_km(ex[sub], km)
    m.change_bounds(bio, 0.0, cap)
    if MOTION == "convection":
        m.add_convection_parameters(packedDensity=PACKED[mid], elasticModulus=ELASTIC[mid],
                                    frictionConstant=FRICTION, convDiffConstant=CONVDIFF)
    return m


cb_cmodel = make(cb_model, "c_beijerinckii", cb_kinetics, cb_ex, cb_bio, 0.4)
mb_cmodel = make(mb_model, "m_barkeri", mb_kinetics, mb_ex, mb_bio, 0.12)

# grid: n interior boxes over 120 um + 1-box static ring outside
DOMAIN_UM = 200.0 if SCENARIO == "fig6" else 120.0
n_in = int(round(DOMAIN_UM / GRID_UM)); assert abs(n_in * GRID_UM - DOMAIN_UM) < 1e-9, "GRID_UM must divide the domain"
N = n_in + 2
dx_cm = GRID_UM * 1e-4
V_box_L = dx_cm ** 3 * 1e-3

pop = {"c_beijerinckii": {}, "m_barkeri": {}}
if SCENARIO == "fig6":
    # one PhysiCelldFBA cell -> its dry mass in the box that contains it (domain from -DOMAIN_UM/2)
    cells = pd.read_csv(os.environ.get("CELLS", "cells_CF3-REF_s50.csv"))
    VOL = {"c_beijerinckii": 1.3, "m_barkeri": 1.8}   # um3, CF-REF_s50.xml / CF3-REF_s50.xml
    for r in cells.itertuples():
        nx = int((r.x + DOMAIN_UM / 2) // GRID_UM) + 1; ny = int((r.y + DOMAIN_UM / 2) // GRID_UM) + 1
        pop[r.type][(nx, ny)] = pop[r.type].get((nx, ny), 0.0) + VOL[r.type] * (1 - 0.75) * 1.04e-12
    for sp_ in pop:
        want = (cells.type == sp_).sum() * VOL[sp_] * 0.25 * 1.04e-12
        assert abs(sum(pop[sp_].values()) - want) <= 1e-9 * want, (sp_, sum(pop[sp_].values()), want)
else:
    # old 5 um boxes (bx, by) of the 24 x 24 grid covering the same 120 um -> new interior box
    OLD_UM = 5.0
    src = json.load(open("cf_ref_box_biomass.json"))
    for key, d in src.items():
        bx, by = [int(v) for v in key.split(",")]
        nx, ny = int((bx + 0.5) * OLD_UM // GRID_UM) + 1, int((by + 0.5) * OLD_UM // GRID_UM) + 1
        for sp_, v in d.items():
            pop[sp_][(nx, ny)] = pop[sp_].get((nx, ny), 0.0) + v * BIOSCALE * 1e-12
    for sp_ in pop:
        want = sum(d.get(sp_, 0.0) for d in src.values()) * BIOSCALE * 1e-12
        assert abs(sum(pop[sp_].values()) - want) <= 1e-9 * want, (sp_, sum(pop[sp_].values()), want)
cb_cmodel.initial_pop = [[x, y, b] for (x, y), b in pop["c_beijerinckii"].items()]
mb_cmodel.initial_pop = [[x, y, b] for (x, y), b in pop["m_barkeri"].items()]
print("grid", N, "x", N, "(interior", n_in, ") boxes; initial boxes cb", sorted(pop["c_beijerinckii"]), "mb", sorted(pop["m_barkeri"]))
print("initial biomass g: cb", sum(pop["c_beijerinckii"].values()), "mb", sum(pop["m_barkeri"].values()))

layout = cm.layout([cb_cmodel, mb_cmodel])
layout.grid = [N, N]
bc_mM = {"glucose": 1, "CO2": 1, "NH4": 5, "Pi": 5, "SO4": 5, "H2": 0, "acetate": 0, "methane": 0}   # boundary (ring)
ic_mM = dict(bc_mM, **({"H2": 1, "acetate": 1} if SCENARIO == "fig6" else {}))                  # initial, inside
met_name = {"glucose": "glc_D_e", "CO2": "co2_e", "NH4": "nh4_e", "Pi": "pi_e",
            "SO4": "so4_e", "H2": "h2_e", "acetate": "ac_e", "methane": "ch4_e"}
diff_um2_min = {"glucose": 40000.0, "methane": 90000, "CO2": 115000, "H2": 270000,
                "NH4": 118000, "Pi": 53000, "SO4": 64000, "acetate": 72000}
for key, mM in ic_mM.items():
    layout.set_specific_metabolite(met_name[key], mM * V_box_L)
    layout.set_specific_metabolite_diffusion(met_name[key], diff_um2_min[key] * 1e-8 / 60.0)

# all other exchanged metabolites at 1000 mM (as PhysiCelldFBA leaves unmapped exchanges at their SBML bounds)
curated = set(met_name.values())
all_ex_mets = set()
for model in (cb_model, mb_model):
    for r in model.reactions:
        if r.id.startswith("R_EX_") or r.id.startswith("EX_"):
            all_ex_mets.update(m.id for m in r.metabolites)
for met_id in sorted(all_ex_mets - curated):
    try:
        layout.set_specific_metabolite(met_id, 1000.0 * V_box_L)
    except Exception as e:
        print("skip", met_id, e)

# Dirichlet: static ring outside the 120 um domain
for key, mM in bc_mM.items():
    for i in range(N):
        for j in (0, N - 1):
            layout.set_specific_static_at_location(met_name[key], (i, j), mM * V_box_L)
            layout.set_specific_static_at_location(met_name[key], (j, i), mM * V_box_L)

# Media diffusion: COMETS uses an ADI scheme (Utility.diffuseFick: alternating diffuseHorizontal / diffuseVertical), which
# is stable but oscillates (and the negative values are clipped to 0) when r = D dt_sub / dx^2 > 1. With the default 10
# substeps, H2 at 10 um / 0.25 min has r = 68 and the field became a checkerboard of ~1.6 mM and 0 mM boxes (run of
# 2026-11-19). numDiffPerStep is therefore set so that r <= R_DIFF (default 1) for the fastest metabolite.
R_DIFF = float(os.environ.get("R_DIFF", "1"))
D_max_cm2_s = max(diff_um2_min.values()) * 1e-8 / 60.0
NDIFF = int(os.environ.get("NDIFF", max(10, int(-(-D_max_cm2_s * TIMESTEP_MIN * 60 / dx_cm ** 2 // R_DIFF)))))
print("numDiffPerStep", NDIFF, "r_max", D_max_cm2_s * TIMESTEP_MIN * 60 / NDIFF / dx_cm ** 2)
cyc = lambda minutes: max(1, int(round(minutes / TIMESTEP_MIN)))
params = cm.params()
params.set_param("timeStep", TIMESTEP_MIN / 60.0)   # hours
params.set_param("spaceWidth", dx_cm)
params.set_param("maxCycles", int(os.environ.get("CF_MAXCYCLES", cyc(T_END_H * 60))))
params.set_param("minSpaceBiomass", MINSPACE)
params.set_param("defaultVmax", DEFVMAX)
params.set_param("defaultKm", DEFKM)
params.set_param("cellSize", 1e-16)
params.set_param("growthDiffRate", 0.0)
params.set_param("numRunThreads", NTHREADS)
params.set_param("numDiffPerStep", NDIFF)
if MOTION == "none":
    params.set_param("flowDiffRate", 0.0)
else:
    params.set_param("biomassMotionStyle", "Convection 2D")
    params.set_param("flowDiffRate", 0.0)
for k, v in [("writeTotalBiomassLog", True), ("totalBiomassLogRate", cyc(LOG_MIN)),
             ("writeBiomassLog", True), ("BiomassLogRate", cyc(LOG_MIN)),
             ("writeMediaLog", True), ("MediaLogRate", cyc(LOG_MIN)),
             ("writeFluxLog", True), ("FluxLogRate", cyc(FLUX_H * 60))]:
    params.set_param(k, v)

wd = "./"
layout.write_necessary_files(wd)
params.write_params(wd + "global_params.txt", wd + "package_params.txt")
layout_fname = glob.glob(wd + ".current_layout*")[0].split("/")[-1]
with open(wd + "comets_script", "w") as fh:
    fh.write("\n".join(["load_comets_parameters global_params.txt", "load_package_parameters package_params.txt",
                        f"load_layout {layout_fname}"]))
if os.environ.get("WRITE_ONLY") == "1":
    print("WRITE_ONLY: inputs written, COMETS not run"); raise SystemExit(0)

p = sp.run(f'bash -c "{COMETS_HOME}/bin/comets_scr_fixed ./comets_script"', shell=True, capture_output=True, text=True)
print("rc=", p.returncode)
open("comets_stdout.txt", "w").write(p.stdout + "\n--- STDERR ---\n" + p.stderr)
if p.returncode != 0:
    print("STDERR TAIL:", p.stderr[-2000:]); raise RuntimeError("COMETS run failed")
for f in glob.glob("*_out") + ["biomass"]:
    if os.path.exists(f):
        txt = open(f).read()
        if "," in txt:
            open(f, "w").write(txt.replace(",", "."))
tbm = pd.read_csv(params.all_params["TotalBiomassLogName"], sep=r"\s+", header=None, names=["cycle"] + layout.get_model_ids())
tbm["time_min"] = tbm.cycle * TIMESTEP_MIN
tbm.to_csv(f"comets_cfref_{TAG}_biomass.csv", index=False)
print(tbm.tail(3).to_string())
elapsed = time.time() - t0
json.dump({"elapsed_s": elapsed, "peak_mem_mb": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.0,
           "threads": NTHREADS, "scenario": SCENARIO, "domain_um": DOMAIN_UM, "grid_um": GRID_UM, "grid": N, "motion": MOTION, "timestep_min": TIMESTEP_MIN,
           "minSpaceBiomass": MINSPACE, "bioscale": BIOSCALE, "V_box_L": V_box_L, "numDiffPerStep": NDIFF,
           "packedDensity": PACKED if MOTION == "convection" else None,
           "elasticModulus": ELASTIC if MOTION == "convection" else None, "k_push": K_PUSH if MOTION == "convection" else None,
           "frictionConstant": FRICTION if MOTION == "convection" else None},
          open(f"comets_cfref_{TAG}_timing.json", "w"))
print("TIMING", elapsed)
