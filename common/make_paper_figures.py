"""Build paper_figures/: the figures of the manuscript under their manuscript names, from paper_figures_manifest.csv.
usage: make_paper_figures.py [--paper-figures DIR_OF_THE_MANUSCRIPT_FIGURES] [--out paper_figures]    (run from the repository root)
status in the manifest: identical = the regenerated output is the file of the manuscript; replace = the regenerated output differs and is to replace it (copied here);
older_runs / external = not an output of the pinned runs: the manuscript's own file is copied (needs --paper-figures) and flagged.
Writes OUT/<manuscript_file> and OUT/MANIFEST.csv (manuscript_file, status, copied_from, sha256, differs_from_manuscript, note)."""
import sys, os, csv, shutil, hashlib, argparse
ap = argparse.ArgumentParser(); ap.add_argument("--paper-figures", default=None); ap.add_argument("--out", default="paper_figures"); a = ap.parse_args()
sha = lambda p: hashlib.sha256(open(p, "rb").read()).hexdigest()
os.makedirs(a.out, exist_ok=True); rows = []; miss = []
for r in csv.DictReader(open("paper_figures_manifest.csv")):
    f = r["manuscript_file"]; ms = os.path.join(a.paper_figures, f) if a.paper_figures else None
    src = r["source_file"] if r["status"] in ("identical", "replace") else (ms or "")
    if not src or not os.path.isfile(src): miss.append((f, src)); rows.append((f, r["status"], src, "", "", "MISSING SOURCE")); continue
    shutil.copy2(src, os.path.join(a.out, f)); h = sha(src); d = ""
    if ms and os.path.isfile(ms): d = "no" if sha(ms) == h else "yes"
    rows.append((f, r["status"], src, h, d, r["note"]))
w = csv.writer(open(os.path.join(a.out, "MANIFEST.csv"), "w", newline=""), lineterminator="\n"); w.writerow(["manuscript_file", "status", "copied_from", "sha256", "differs_from_manuscript", "note"]); w.writerows(rows)
print(len(rows) - len(miss), "of", len(rows), "figures copied to", a.out); [print("  missing:", m) for m in miss]
bad = [r for r in rows if r[1] == "identical" and r[4] == "yes"]; print("identical rows that differ from the manuscript file:", len(bad))
for s in ("identical", "replace", "older_runs", "external"): print(f"  {s:11s}", sum(1 for r in rows if r[1] == s))
