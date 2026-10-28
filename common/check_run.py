"""Check that a run is what the tables say: the PhysiCell_settings.xml saved in the run's output must equal (canonical XML) the config
rebuilt from base + overrides.csv by build_configs.py. Exit 1 on any difference.
usage: check_run.py USECASE_DIR RUN_ID"""
import sys, os, subprocess, tempfile
from lxml import etree
uc, rid = sys.argv[1:3]
def canon(p):
    t = etree.parse(p, etree.XMLParser(remove_blank_text=True, remove_comments=True))
    for e in t.iter():
        e.text = e.text.strip() if e.text and e.text.strip() else None; e.tail = None
    return etree.tostring(t, method="c14n")
with tempfile.TemporaryDirectory() as d:
    subprocess.run([sys.executable, os.path.join(os.path.dirname(__file__), "build_configs.py"), uc, rid, "--out", d], check=True, stdout=subprocess.DEVNULL)
    same = canon(os.path.join(d, rid + ".xml")) == canon(os.path.join(uc, "runs", rid, "PhysiCell_settings.xml"))
print(rid, "OK" if same else "MISMATCH: saved settings differ from the rebuilt config"); sys.exit(0 if same else 1)
