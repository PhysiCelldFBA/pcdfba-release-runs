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
with open(f"{UC}/runs_def.csv", "w", newline="") as f: w = csv.writer(f, lineterminator="\n"); w.writerow(["run_id", "experiment_id", "seed"]); w.writerows(runs)
with open(f"{UC}/overrides_def.csv", "w", newline="") as f: w = csv.writer(f, lineterminator="\n"); w.writerow(["scope", "xpath", "value"]); w.writerows(ov)
print(len(runs), "runs,", len(ov), "override rows")
