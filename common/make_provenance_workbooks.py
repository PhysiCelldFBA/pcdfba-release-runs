"""Provenance workbooks: one Excel file per use case (provenance/<uc>_provenance.xlsx) and a CSV copy of every sheet (provenance/<uc>/<sheet>.csv), plus provenance/all_runs.csv.
usage: make_provenance_workbooks.py [USECASE ...]     (run from the repository root; default: all use cases with a provenance.csv)
Sheets: README; experiments (experiments.csv); runs (provenance.csv + the planned run's experiment, seed and title); overrides (the rows of every overrides*.csv: scope, xpath, value, file);
binaries (binaries.tsv); parameters_varying (every parameter of the PhysiCell settings that differs between runs of the use case: one row per parameter, one column per run); parameters_full (all parameters of every run).
Parameters are read from runs/<run>/PhysiCell_settings.xml, the configuration the run actually used (leaf elements and their attributes, except units and description), as XPath-like keys, repeated elements keyed by name / ID / substrate / type.
Run folders without a planned run are listed in README. Nothing in runs/ is modified."""
import sys, os, csv, glob, datetime, subprocess, xml.etree.ElementTree as ET
import pandas as pd
SKIP_ATTR = {"units", "description", "name", "ID", "substrate", "source"}
def flatten(path):
    root = ET.parse(path).getroot(); out = {}
    def key(e, siblings):
        for a in ("name", "ID", "substrate", "type"):
            if e.get(a) is not None: return f"{e.tag}[{a}={e.get(a)}]"
        same = [s for s in siblings if s.tag == e.tag]
        return e.tag if len(same) == 1 else f"{e.tag}[{same.index(e) + 1}]"
    def walk(e, p):
        for a, v in e.attrib.items():
            if a not in SKIP_ATTR and not (a == "type" and len(e) == 0 and False): out[f"{p}/@{a}"] = v
        kids = list(e)
        if not kids:
            if (e.text or "").strip(): out[p] = e.text.strip()
        for k in kids: walk(k, f"{p}/{key(k, kids)}")
    walk(root, root.tag); return out
def build(uc):
    base = os.path.join("usecases", uc); rd = lambda n: pd.read_csv(os.path.join(base, n), dtype=str, keep_default_na=False)
    exp = rd("experiments.csv"); prov = rd("provenance.csv")
    planned = pd.concat([rd(os.path.basename(f)) for f in sorted(glob.glob(os.path.join(base, "runs*.csv"))) if os.path.basename(f) != "provenance.csv"], ignore_index=True).drop_duplicates("run_id")
    title = dict(zip(exp.experiment_id, exp.title)); runs = prov.merge(planned[["run_id", "seed"]].rename(columns={"seed": "seed_planned"}), on="run_id", how="left"); runs.insert(2, "experiment_title", runs.experiment_id.map(title).fillna(""))
    ov = []
    for f in sorted(glob.glob(os.path.join(base, "overrides*.csv"))):
        d = rd(os.path.basename(f)); d["file"] = os.path.basename(f); ov.append(d)
    ov = pd.concat(ov, ignore_index=True) if ov else pd.DataFrame(columns=["scope", "xpath", "value", "file"])
    bn = pd.read_csv(os.path.join(base, "binaries.tsv"), sep="\t", header=None, comment="#", names=["binary_sha", "commit", "note"], dtype=str, keep_default_na=False)
    cols, missing = {}, []
    for rid in prov.run_id:
        p = os.path.join(base, "runs", rid, "PhysiCell_settings.xml")
        if not os.path.isfile(p): p2 = os.path.join(base, "config", rid + ".xml"); p = p2 if os.path.isfile(p2) else None
        if p is None: missing.append(rid); continue
        cols[rid] = flatten(p)
    full = pd.DataFrame(cols); full.index.name = "parameter"; full = full.sort_index()
    var = full[full.drop([i for i in full.index if i.endswith("/save/folder")], errors="ignore").nunique(axis=1, dropna=False).reindex(full.index).fillna(0) > 1]
    extra = sorted(set(os.listdir(os.path.join(base, "runs"))) - set(prov.run_id))
    try: code = subprocess.run(["git", "-C", "code", "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()
    except Exception: code = ""
    readme = pd.DataFrame({"item": ["use case", "generated", "pinned code commit (submodule code/)", "runs in provenance.csv", "runs with a configuration read", "parameters (all / varying between runs)", "run folders without a planned run", "sheets"],
        "value": [uc, datetime.date.today().isoformat() + " by common/make_provenance_workbooks.py", code, len(prov), len(cols), f"{len(full)} / {len(var)}", "; ".join(extra) if extra else "none",
                  "experiments: what each experiment asks and its pass criterion; runs: what happened to each run (commit, binary, config hash, machine, times, rc, status); overrides: the changes on top of the base configuration; binaries: binary hash to commit; parameters_varying / parameters_full: the configuration each run used (see the script docstring)"]})
    if missing: readme.loc[len(readme)] = ["runs without a configuration file", "; ".join(missing)]
    return {"README": readme, "experiments": exp, "runs": runs, "overrides": ov, "binaries": bn, "parameters_varying": var.reset_index(), "parameters_full": full.reset_index()}
def main():
    ucs = sys.argv[1:] or sorted(u for u in os.listdir("usecases") if os.path.isfile(os.path.join("usecases", u, "provenance.csv"))); os.makedirs("provenance", exist_ok=True); allruns = []
    for uc in ucs:
        sheets = build(uc); os.makedirs(os.path.join("provenance", uc), exist_ok=True)
        with pd.ExcelWriter(os.path.join("provenance", f"{uc}_provenance.xlsx"), engine="openpyxl") as xw:
            for n, df in sheets.items():
                df.to_excel(xw, sheet_name=n, index=False); ws = xw.sheets[n]; ws.freeze_panes = "B2" if n.startswith("parameters") else "A2"
                for i, c in enumerate(df.columns):
                    w = min(60, max(len(str(c)), *(len(str(v)) for v in df[c].head(200))) + 2) if len(df) else len(str(c)) + 2
                    ws.column_dimensions[ws.cell(1, i + 1).column_letter].width = w if n != "parameters_full" and n != "parameters_varying" or i == 0 else 22
                df.to_csv(os.path.join("provenance", uc, n + ".csv"), index=False)
        r = sheets["runs"].copy(); r.insert(0, "use_case", uc); allruns.append(r)
        print(f"{uc:13s} runs {len(sheets['runs']):4d}  parameters {len(sheets['parameters_full']):5d}  varying {len(sheets['parameters_varying']):4d}")
    pd.concat(allruns, ignore_index=True).to_csv(os.path.join("provenance", "all_runs.csv"), index=False); print("wrote provenance/")
if __name__ == "__main__": main()
