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
# CF-DX: voxel size at dt 0.001
for dx in (0.8, 1.6): rid = f"CF-DX_dx{dx}_dt0.001_s50"; run(rid, "CF-DX"); setdx(rid, dx)
# CF-DT: time step at dx 1.25
for dt in (0.002, 0.0005): rid = f"CF-DT_dx1.25_dt{dt}_s50"; run(rid, "CF-DT"); setdt(rid, dt)
# CF-DXDT: finer step with the other voxel sizes
for dx in (0.8, 1.6): rid = f"CF-DXDT_dx{dx}_dt0.0005_s50"; run(rid, "CF-DXDT"); setdx(rid, dx); setdt(rid, 0.0005)
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
with open(f"{UC}/runs_def.csv", "w", newline="") as f: w = csv.writer(f, lineterminator="\n"); w.writerow(["run_id", "experiment_id", "seed"]); w.writerows(runs)
with open(f"{UC}/overrides_def.csv", "w", newline="") as f: w = csv.writer(f, lineterminator="\n"); w.writerow(["scope", "xpath", "value"]); w.writerows(ov)
print(len(runs), "runs,", len(ov), "override rows")
