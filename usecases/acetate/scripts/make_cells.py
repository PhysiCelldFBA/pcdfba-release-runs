"""Initial cells of the acetate runs: n cells uniform in the box (1 um from the walls), seeded by the run seed.
usage: make_cells.py CONFIG.xml OUT.csv '{"cells.n": "20"}'   (box, per axis, and seed are read from the generated config; 1 um from every wall)"""
import sys, json, random
from lxml import etree
cfg, out, extra = sys.argv[1], sys.argv[2], json.loads(sys.argv[3])
r = etree.parse(cfg).getroot(); g = lambda p: float(r.xpath(p)[0].text)
box = [(g(f"/PhysiCell_settings/domain/{a}_min") + 1, g(f"/PhysiCell_settings/domain/{a}_max") - 1) for a in "xyz"]; seed = int(r.xpath("/PhysiCell_settings/options/random_seed")[0].text)
rnd = random.Random(seed); rows = ["x,y,z,type"]
for _ in range(int(extra["cells.n"])): rows.append("%.4f,%.4f,%.4f,ecoli" % tuple(rnd.uniform(lo, hi) for lo, hi in box))
open(out, "w").write("\n".join(rows) + "\n")
