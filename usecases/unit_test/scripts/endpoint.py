"""Closed-batch endpoint of the mass-conversion test: simulated final dry biomass against the FBA-derived analytical endpoint (Supplementary Methods S4).
usage: endpoint.py OUT.csv RUNDIR [RUNDIR ...]
Analytical endpoint = initial biomass + (mu / q_glc) x glucose in the box, with mu and q_glc from FBA of the run's own model with every exchange listed in its configuration set to its Vmax (lower bound = -Vmax);
initial biomass = simulated biomass at t = 0 (cells x volume x (1 - fluid fraction) x 1.04 pg/um3); glucose in the box = initial concentration x box volume (4096 um3 -> use the domain of the run)."""
import sys, glob, os, csv, numpy as np, scipy.io as sio, xml.etree.ElementTree as ET, cobra, warnings
warnings.filterwarnings("ignore")
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..", "common")); from pcxml import labels
UC = os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."); rows = []
for d in (x.rstrip("/") for x in sys.argv[2:]):
    cfg = ET.parse(f"{d}/PhysiCell_settings.xml").getroot()
    ex = {e.get("substrate"): (e.find("fba_flux").text.strip().replace("R_", "", 1), float(e.find("Vmax").text)) for e in cfg.findall(".//transport_model/exchange")}
    model = cobra.io.read_sbml_model(os.path.join(UC, cfg.find(".//sbml_filename").text.strip()))
    for sub, (rid, vmax) in ex.items(): model.reactions.get_by_id(rid).lower_bound = -vmax
    sol = model.optimize(); mu, qg = sol.objective_value, -sol.fluxes[ex["glucose"][0]]
    g0 = float([v for v in cfg.findall(".//microenvironment_setup/variable") if v.get("name") == "glucose"][0].find("initial_condition").text)
    vol = np.prod([float(cfg.find(f".//domain/{a}_max").text) - float(cfg.find(f".//domain/{a}_min").text) for a in ("x", "y", "z")])
    bio, nc = [], []
    for s in sorted(glob.glob(f"{d}/output*[0-9].xml")):
        lab = labels(ET.parse(s).getroot()); m = sio.loadmat(s.replace(".xml", "_cells.mat"))["cells"]
        bio.append((m[lab.index("total_volume")] * (1 - m[lab.index("fluid_fraction")])).sum() * 1.04); nc.append(m.shape[1])
    ana = bio[0] + (mu / qg) * (g0 * 1e-3 * vol * 1e-15 * 1e3) * 1e12
    rows.append(dict(run=os.path.basename(d), model=os.path.basename(cfg.find(".//sbml_filename").text.strip()), mu=round(mu, 4), q_glc=round(qg, 3), yield_gdw_per_g=round(mu / qg * 1000 / 180.0, 4), initial_pg=round(bio[0], 3), simulated_pg=round(bio[-1], 3), analytical_pg=round(ana, 3), diff_pct=round(100 * (bio[-1] / ana - 1), 3), cells=nc[-1]))
w = csv.DictWriter(open(sys.argv[1], "w", newline=""), list(rows[0])); w.writeheader(); w.writerows(rows)
for r in rows: print(r)
