"""Write runs_def.csv and overrides_def.csv: the reduced-size cross-feeding runs (CF-*).
Reduced size: box +-60 um instead of +-100 um (120 is a multiple of 1.25, 0.8 and 1.6, so every dx fits an integer number of voxels) and the shipped cell cluster cut to a half-size cluster of the same density (scripts/make_cells.py, 111 cells instead of 353).
Shipped base: dx 1.25, dz 1.25, dt_diffusion = intracellular_dt 0.001, dt_mechanics 0.1, dt_phenotype 6. In every run dz stays 1.25 um (slab thickness constant); only dx = dy changes.
usage: design_def.py"""
import csv, os
UC = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); P = "/PhysiCell_settings"
runs, ov = [], []
run = lambda rid, eid, seed=50: runs.append((rid, eid, seed)); put = lambda s, x, v: ov.append((s, x, v))
CDS = [P + f"/cell_definitions/cell_definition[@name='{n}']/phenotype/intracellular/settings" for n in ("c_beijerinckii", "m_barkeri")]
SB = ("./config/shipped/sbml_models/C_beijerinckii_model_iCB925_trimmed.xml", "./config/shipped/sbml_models/M_barkeri_model_iMG746_trimmed.xml")
# all CF experiments: reduced box, half-size cluster, 1800 min, snapshots every 30 min, no SVG, ledger on, 112 threads
put("*", "@cells.scale", "half")
put("*", P + "/domain/x_min", -60); put("*", P + "/domain/x_max", 60); put("*", P + "/domain/y_min", -60); put("*", P + "/domain/y_max", 60)
put("*", P + "/overall/max_time", 1800); put("*", P + "/save/SVG/enable", "false"); put("*", P + "/save/full_data/interval", 30)
put("*", P + "/parallel/omp_num_threads", 56)
put("*", P + "/initial_conditions/cell_positions/folder", "./config")
for c, s in zip(CDS, SB):
    put("*", c + "/sbml_filename", s); put("*", "+" + c, "<dfba_conservation_ledger>true</dfba_conservation_ledger>")
def setdt(rid, dt):
    put(rid, P + "/overall/dt_diffusion", dt)
    for c in CDS: put(rid, c + "/intracellular_dt", dt)
def setdx(rid, dx):
    put(rid, P + "/domain/dx", dx); put(rid, P + "/domain/dy", dx)
# CF-REF: dx 1.25, dt 0.001 (shipped), seeds 50 51 52 and the same-seed repeat
for rid, seed in (("CF-REF_s50", 50), ("CF-SEED_s51", 51), ("CF-SEED_s52", 52), ("CF-REP_s50", 50)): run(rid, rid.split("_")[0], seed)
# CF-DX: voxel size at dt 0.001 (dx 1.6 done; dx 0.8, below the cell size, is not run: dropped 2026-11-04, the series continues as CF3)
for dx in (1.6,): rid = f"CF-DX_dx{dx}_dt0.001_s50"; run(rid, "CF-DX"); setdx(rid, dx)
# CF-DT: time step at dx 1.25 (dt 0.002 done; dt 0.0005 not run: dropped 2026-11-04, the finest step is run as CF3-DT at dx 1.6)
for dt in (0.002,): rid = f"CF-DT_dx1.25_dt{dt}_s50"; run(rid, "CF-DT"); setdt(rid, dt)
# CF-DXDT (dx 0.8 and 1.6 at dt 0.0005) was never run and is dropped 2026-11-04; its role is taken by CF3-DT
# CF-FIG: preview of the manuscript figure from one run (dx 1.6, dt 0.001) on a debug node, 112 threads, own run id (the same setting is CF-DX_dx1.6_dt0.001_s50)
rid = "CF-FIG_dx1.6_dt0.001_s50"; run(rid, "CF-FIG"); setdx(rid, 1.6); put(rid, P + "/parallel/omp_num_threads", 112)
# CF2: second design after the first series (RESULTS.md): the +-60 um box and the zero boundaries made H2 / acetate leave through the edges.
# Box +-100 um (the paper box, 200 is a multiple of 1.25, 0.8 and 1.6) with the same 111 cells; H2 and acetate start at 1 mM (CO2 already starts at 1 mM and stays held at 1 mM at the boundary as shipped; methane 0, glucose and mineral nutrients as shipped);
# the H2, acetate and methane boundaries stay 0 as shipped. CF2-KM also lowers the Km of M. barkeri for H2 from 0.05 to 0.01 mM.
for rid, eid, km in (("CF2-BOX_s50", "CF2-BOX", None), ("CF2-KM_s50", "CF2-KM", 0.01)):
    run(rid, eid)
    for ax in ("x", "y"): put(rid, P + f"/domain/{ax}_min", -100); put(rid, P + f"/domain/{ax}_max", 100)
    for sub in ("H2", "acetate"): put(rid, P + f"/microenvironment_setup/variable[@name='{sub}']/initial_condition", 1)
    if km: put(rid, P + "/cell_definitions/cell_definition[@name='m_barkeri']/phenotype/intracellular/transport_model/exchange[@substrate='H2']/Km", km)
# CF3: convergence series on the CF2-KM design (box +-100 um, H2 / acetate initially 1 mM, Km of M. barkeri for H2 0.01 mM), with voxels not smaller than the largest living cell.
# Largest living cell volume in CF2-KM: C. beijerinckii 2.26, M. barkeri 2.98 um3 (dz 1.25 um), so the voxel volume dx*dx*1.25 must be >= 2.98: dx >= 1.55 um; dx 1.25 (1.95 um3) is outside, dx 0.8 is not run.
# Reference dx 1.6 (3.2 um3), dt 0.001; seed 51; dx 2.0 (5.0 um3); dt 0.002 and 0.0005 at dx 1.6. 200 is a multiple of 1.6 and 2.0.
KM = P + "/cell_definitions/cell_definition[@name='m_barkeri']/phenotype/intracellular/transport_model/exchange[@substrate='H2']/Km"
def cf3(rid, eid, seed=50, dx=1.6, dt=0.001):
    run(rid, eid, seed)
    for ax in ("x", "y"): put(rid, P + f"/domain/{ax}_min", -100); put(rid, P + f"/domain/{ax}_max", 100)
    for sub in ("H2", "acetate"): put(rid, P + f"/microenvironment_setup/variable[@name='{sub}']/initial_condition", 1)
    put(rid, KM, 0.01); setdx(rid, dx)
    if dt != 0.001: setdt(rid, dt)
cf3("CF3-REF_s50", "CF3-REF"); cf3("CF3-SEED_s51", "CF3-SEED", 51); cf3("CF3-DX_dx2.0_dt0.001_s50", "CF3-DX", dx=2.0)
for dt in (0.002, 0.0005): cf3(f"CF3-DT_dx1.6_dt{dt}_s50", "CF3-DT", dt=dt)
# CF-NC / CF3-NC (2026-11-06): the same two designs with the trimmed models without the FVA caps (models/*_trimmed_nocaps.xml, built by unit_test/scripts/make_nocaps_models.py:
# trimmed topology, internal bounds of the full models, exchange bounds of the medium unchanged, ATPM lower bound 0). Everything else as CF-REF and CF3-REF, seed 50.
NC = ("./models/C_beijerinckii_model_iCB925_trimmed_nocaps.xml", "./models/M_barkeri_model_iMG746_trimmed_nocaps.xml")
def nocaps(rid):
    for c, f in zip(CDS, NC): put(rid, c + "/sbml_filename", f)
run("CF-NC_s50", "CF-NC"); nocaps("CF-NC_s50")
cf3("CF3-NC_s50", "CF3-NC"); nocaps("CF3-NC_s50")
# CF3-NF (2026-11-06): CF3-NC with no-flux boundaries for H2 and acetate (Dirichlet condition off on every face; glucose, CO2 and methane as shipped) and H2 / acetate initially 0, so the only source of H2 and acetate is C. beijerinckii and nothing leaves the box.
# CF3-NC50: CF3-NC run for 3000 min (50 h, the duration in the manuscript text) to measure the cost and the late growth of M. barkeri.
cf3("CF3-NF_s50", "CF3-NF"); nocaps("CF3-NF_s50")
for sub in ("H2", "acetate"):
    V = P + f"/microenvironment_setup/variable[@name='{sub}']"
    put("CF3-NF_s50", V + "/initial_condition", 0)  # with a closed box the initial 1 mM would be an H2 / acetate pool for M. barkeri: start at 0, so C. beijerinckii is the only source
    put("CF3-NF_s50", V + "/Dirichlet_boundary_condition/@enabled", "False")
    for f in ("xmin", "xmax", "ymin", "ymax", "zmin", "zmax"): put("CF3-NF_s50", V + f"/Dirichlet_options/boundary_value[@ID='{f}']/@enabled", "False")
# CF3-DB (2026-11-07): CF3-NC with H2 and acetate held at 0.001 mM on every face (Dirichlet) and initially 0.001 mM, between the open box at 0 (CF3-NC) and the closed box (CF3-NF).
# 0.001 mM is of the order of the concentration that the open box itself keeps one voxel inside the edge at 0.5 h (H2 0.0017, acetate 0.0008 mM; scripts/border_profile.py) and is a tenth of the Km of M. barkeri for H2.
cf3("CF3-DB_s50", "CF3-DB"); nocaps("CF3-DB_s50")
for sub in ("H2", "acetate"):
    V = P + f"/microenvironment_setup/variable[@name='{sub}']"
    put("CF3-DB_s50", V + "/initial_condition", 0.001); put("CF3-DB_s50", V + "/Dirichlet_boundary_condition", 0.001)
    for f in ("xmin", "xmax", "ymin", "ymax", "zmin", "zmax"): put("CF3-DB_s50", V + f"/Dirichlet_options/boundary_value[@ID='{f}']", 0.001)
cf3("CF3-NC50_s50", "CF3-NC50"); nocaps("CF3-NC50_s50"); put("CF3-NC50_s50", P + "/overall/max_time", 3000)
# CF3-NC-DT (2026-11-07): the dt series of the main scenario (no caps, open at 0 mM): dt 0.002 and 0.0005 min at dx 1.6, seed 50, to measure the dependence on dt where M. barkeri is capped (cap ratio 0.14 at 10 h in CF3-NC and CF3-DB).
for dt in (0.002, 0.0005): rid = f"CF3-NC-DT_dx1.6_dt{dt}_s50"; cf3(rid, "CF3-NC-DT", dt=dt); nocaps(rid)
# CF3-BR (2026-11-07): H2 and acetate held on every face and initially at the value one voxel inside the edge of the open box (CF3-NC) when M. barkeri makes its first new cell (1050 min): H2 0.0238 -> 0.024 mM, acetate 0.0129 -> 0.013 mM
# (CF3-NC, scripts/border_profile.py, results/cf3_border.csv). Above the Km of M. barkeri for H2 (0.01 mM): the surroundings are pre-loaded, a bracketing case and not a test of cross-feeding.
cf3("CF3-BR_s50", "CF3-BR"); nocaps("CF3-BR_s50")
for sub, v in (("H2", 0.024), ("acetate", 0.013)):
    V = P + f"/microenvironment_setup/variable[@name='{sub}']"
    put("CF3-BR_s50", V + "/initial_condition", v); put("CF3-BR_s50", V + "/Dirichlet_boundary_condition", v)
    for f in ("xmin", "xmax", "ymin", "ymax", "zmin", "zmax"): put("CF3-BR_s50", V + f"/Dirichlet_options/boundary_value[@ID='{f}']", v)
# CF3-NC-AC8 (2026-11-13): CF3-NC (same seed, 1800 min) with the acetate Vmax of M. barkeri 8 instead of 3 mmol/gDW/h (iAF692 maximum acetate uptake, Feist et al. 2006; the replay of CF3-NC50 cells found every living M. barkeri cell at the 3.0 bound).
cf3("CF3-NC-AC8_s50", "CF3-NC-AC8"); nocaps("CF3-NC-AC8_s50")
put("CF3-NC-AC8_s50", P + "/cell_definitions/cell_definition[@name='m_barkeri']/phenotype/intracellular/transport_model/exchange[@substrate='acetate']/Vmax", 8.0)
# CF3-NC-ACVMAX, CF3-NC-ACKM, CF3-NC-DIFF (2026-11-18): sensitivity of the main scenario (CF3-NC: no caps, open at 0 mM, dx 1.6, dt 0.001, seed 50, 1800 min) for R2.M3 / R3.M3.
# CF3-NC-ACVMAX: acetate Vmax of M. barkeri 3 -> 0.75 (x0.25; CF3-NC-AC8 is the upward arm). CF3-NC-ACKM: acetate Km of M. barkeri 0.0001 -> 0.01 and 0.3 mM.
# CF3-NC-DIFF: diffusion coefficients of H2 (270000 um2/min) and acetate (72000) together x0.5 and x2; glucose, CO2, methane and the minerals as shipped.
AC = P + "/cell_definitions/cell_definition[@name='m_barkeri']/phenotype/intracellular/transport_model/exchange[@substrate='acetate']/"
# CF3-NC-CTL: same-config rerun of CF3-NC in the same job (the binary of CF3-NC_s50, 1f60748dc4ee2e45, is no longer on MN5; this job uses f9c90c60aec019bd), own run id as CF-REP.
cf3("CF3-NC-CTL_s50", "CF3-NC-CTL"); nocaps("CF3-NC-CTL_s50")
rid = "CF3-NC-ACVMAX_0.75_s50"; cf3(rid, "CF3-NC-ACVMAX"); nocaps(rid); put(rid, AC + "Vmax", 0.75)
for km in (0.01, 0.3): rid = f"CF3-NC-ACKM_{km}_s50"; cf3(rid, "CF3-NC-ACKM"); nocaps(rid); put(rid, AC + "Km", km)
DIF = lambda s: P + f"/microenvironment_setup/variable[@name='{s}']/physical_parameter_set/diffusion_coefficient"
for tag, f in (("x0.5", 0.5), ("x2", 2)):
    rid = f"CF3-NC-DIFF_{tag}_s50"; cf3(rid, "CF3-NC-DIFF"); nocaps(rid)
    for s, d in (("H2", 270000), ("acetate", 72000)): put(rid, DIF(s), int(d * f))
with open(f"{UC}/runs_def.csv", "w", newline="") as f: w = csv.writer(f, lineterminator="\n"); w.writerow(["run_id", "experiment_id", "seed"]); w.writerows(runs)
with open(f"{UC}/overrides_def.csv", "w", newline="") as f: w = csv.writer(f, lineterminator="\n"); w.writerow(["scope", "xpath", "value"]); w.writerows(ov)
print(len(runs), "runs,", len(ov), "override rows")
