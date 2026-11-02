"""Initial cells of the cross-feeding runs, taken from the shipped lattice (config/shipped/200cells.csv, 353 cells: C. beijerinckii x -24..-5, M. barkeri x 6..25 um, y -9..10 um).
usage: make_cells.py CONFIG.xml OUT.csv '{"cells.scale": "half"}'
  cells.scale = full   the shipped table, unchanged
  cells.scale = half   a half-size cluster at the same lattice spacing (same cell density) and the same gap between the two species: x kept within about 11-12 um of the gap
                       on each side, y rows of the lower half, shifted up to centre on y = 0; ~1/4 of the cells, to go with a box of half the side"""
import sys, json, os, csv
cfg, out, extra = sys.argv[1], sys.argv[2], json.loads(sys.argv[3])
src = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "config", "shipped", "200cells.csv")
rows = list(csv.DictReader(open(src)))
scale = extra.get("cells.scale", "full")
if scale == "half":
    ys = sorted({float(r["y"]) for r in rows}); ylo, yhi = ys[0], ys[len(ys) // 2 - 1]; shift = -(ylo + yhi) / 2
    keep = []
    for r in rows:
        x, y = float(r["x"]), float(r["y"])
        if y > yhi + 1e-9: continue
        if r["type"] == "c_beijerinckii" and x < -17.5: continue
        if r["type"] == "m_barkeri" and x > 16.0: continue
        r["y"] = repr(round(y + shift, 6)); keep.append(r)
    rows = keep
elif scale != "full": sys.exit("cells.scale must be full or half")
with open(out, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys()), lineterminator="\n"); w.writeheader(); w.writerows(rows)
