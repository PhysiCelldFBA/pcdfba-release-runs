"""Build the PhysiCell config of every planned run of a use case: base config + overrides, nothing else.
usage: build_configs.py USECASE_DIR [RUN_ID ...] [--out DIR]     (default out: USECASE_DIR/config)

Tables (all in USECASE_DIR):
  experiments.csv  experiment_id, title, question, feeds, base, rationale     base = path of the base PhysiCell_settings.xml, relative to USECASE_DIR
  runs*.csv        run_id, experiment_id, seed                                 the planned runs (what was intended)
  overrides*.csv   scope, xpath, value                                         scope = * (all runs) | experiment_id | run_id; applied in this order
                   (files runs_*.csv / overrides_*.csv are written by scripts/design_*.py; all matching files are read)
Rules: every xpath must match exactly one element (a typo or a changed base fails loudly, never silently).
  Derived, never listed in overrides.csv:  options/random_seed = the run's seed (and the user parameter of legacy bases);
  save/folder = ./runs/<run_id>;  initial_conditions/cell_positions/filename = cells_<run_id>.csv (only if scripts/make_cells.py exists).
  A key starting with + adds an element (as the first child): xpath = +PARENT_XPATH, value = the XML fragment; it fails if the parent is not unique or the
  element already exists (then the base has caught up and the row must go).
  Keys starting with @ are not XML: they are arguments for USECASE_DIR/scripts/make_cells.py (@cells.n = number of cells).
Writes <out>/<run_id>.xml and <out>/cells_<run_id>.csv."""
import sys, os, csv, json, subprocess, glob
from lxml import etree
args = [a for a in sys.argv[1:] if not a.startswith("--")]; uc = args[0]; only = args[1:]
out = sys.argv[sys.argv.index("--out") + 1] if "--out" in sys.argv else os.path.join(uc, "config")
if "--out" in sys.argv: only = [a for a in only if a != out]
rd = lambda n: list(csv.DictReader(open(os.path.join(uc, n))))
rdg = lambda stem: [r for f in sorted(glob.glob(os.path.join(uc, stem + "*.csv"))) for r in csv.DictReader(open(f))]
exps = {e["experiment_id"]: e for e in rd("experiments.csv")}; runs = rdg("runs"); ov = rdg("overrides")
os.makedirs(out, exist_ok=True)
def setone(root, xp, val, ctx):
    m = root.xpath(xp)
    if len(m) != 1: sys.exit(f"{ctx}: xpath {xp} matches {len(m)} elements (must be exactly 1)")
    m[0].text = str(val)
for r in runs:
    rid, eid = r["run_id"], r["experiment_id"]
    if only and rid not in only: continue
    if eid not in exps: sys.exit(f"{rid}: experiment {eid} not in experiments.csv")
    tree = etree.parse(os.path.join(uc, exps[eid]["base"]), etree.XMLParser(remove_blank_text=False)); root = tree.getroot()
    extra = {}
    for scope in ("*", eid, rid):
        for o in (o for o in ov if o["scope"] == scope):
            if o["xpath"].startswith("@"): extra[o["xpath"][1:]] = o["value"]
            elif o["xpath"].startswith("+"):
                par = root.xpath(o["xpath"][1:]); new = etree.fromstring(o["value"])
                if len(par) != 1: sys.exit(f"{rid} [{scope}]: parent {o['xpath'][1:]} matches {len(par)} elements (must be exactly 1)")
                if par[0].find(new.tag) is not None: sys.exit(f"{rid} [{scope}]: <{new.tag}> already exists under {o['xpath'][1:]}; remove this row")
                par[0].insert(0, new)   # first child, as the earlier generators did
            else: setone(root, o["xpath"], o["value"], f"{rid} [{scope}]")
    setone(root, "/PhysiCell_settings/options/random_seed", r["seed"], rid)
    if root.xpath("/PhysiCell_settings/user_parameters/random_seed"):   # legacy bases only; since #27 the seed is <options><random_seed> alone
        setone(root, "/PhysiCell_settings/user_parameters/random_seed", r["seed"], rid)
    setone(root, "/PhysiCell_settings/save/folder", f"./runs/{rid}", rid)
    if os.path.isfile(os.path.join(uc, "scripts", "make_cells.py")):   # only use cases with a cells hook get a per-run cells file
        setone(root, "/PhysiCell_settings/initial_conditions/cell_positions/filename", f"cells_{rid}.csv", rid)
    cfg = os.path.join(out, rid + ".xml"); tree.write(cfg, xml_declaration=True, encoding="UTF-8")
    hook = os.path.join(uc, "scripts", "make_cells.py")
    if os.path.isfile(hook): subprocess.run([sys.executable, hook, cfg, os.path.join(out, f"cells_{rid}.csv"), json.dumps(extra)], check=True)
    print("built", rid)
