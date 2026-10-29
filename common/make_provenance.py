"""Write usecases/<uc>/provenance.csv: what actually happened to each planned run (runs.csv), read from the run's own output.
usage: make_provenance.py USECASE_DIR
binary -> commit comes from USECASE_DIR/binaries.tsv (binary_sha<TAB>commit<TAB>note); an unlisted binary gets commit '?'.
A planned run with no output is status 'not run'; no end line in run.log is 'partial'; rc != 0 is 'failed'."""
import sys, os, re, csv, glob, hashlib, xml.etree.ElementTree as ET
uc = sys.argv[1]; rd = lambda n: list(csv.DictReader(open(os.path.join(uc, n))))
bt = {l.split("\t")[0]: l.rstrip("\n").split("\t") for l in open(os.path.join(uc, "binaries.tsv")) if l.strip() and not l.startswith("#")}
cols = ["run_id","experiment_id","commit","binary_sha","config_sha","seed","dx","dt_diffusion","intracellular_dt","max_time","threads","machine","start","end","status","output_path"]
rows = []
for r in (x for f in sorted(glob.glob(os.path.join(uc, "runs*.csv"))) for x in csv.DictReader(open(f))):
    rid = r["run_id"]; out = os.path.join(uc, "runs", rid); sf, lf = (os.path.join(out, n) for n in ("PhysiCell_settings.xml", "run.log"))
    row = dict.fromkeys(cols, ""); row.update(run_id=rid, experiment_id=r["experiment_id"], output_path=out)
    if not (os.path.isfile(sf) and os.path.isfile(lf)): row["status"] = "not run"; rows.append(row); continue
    x = ET.parse(sf).getroot(); f = lambda p: (x.find(p).text or "").strip() if x.find(p) is not None else ""
    log = open(lf, errors="replace").read(); m = re.search(r"start=(\S+) binary=(\S+) host=(\S+)", log); e = re.search(r"end=(\S+) rc=(\d+)", log)
    row.update(binary_sha=m.group(2) if m else "", commit=bt.get(m.group(2) if m else "", ["", "?"])[1], config_sha=hashlib.sha256(open(sf, "rb").read()).hexdigest()[:16],
        seed=f(".//options/random_seed"), dx=f(".//domain/dx"), dt_diffusion=f(".//overall/dt_diffusion"), intracellular_dt=f(".//intracellular_dt"),
        max_time=f(".//overall/max_time"), threads=f(".//parallel/omp_num_threads"), machine=m.group(3) if m else "", start=m.group(1) if m else "",
        end=e.group(1) if e else "", status=("ok" if e.group(2) == "0" else "failed") if e else "partial")
    rows.append(row)
with open(os.path.join(uc, "provenance.csv"), "w", newline="") as fh: w = csv.DictWriter(fh, cols); w.writeheader(); w.writerows(rows)
print(len(rows), "rows ->", os.path.join(uc, "provenance.csv"))
