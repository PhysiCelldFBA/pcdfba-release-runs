"""Cells file of the unit-test runs. Only the update-order runs read one (cell_positions enabled in order_base.xml):
the row order of the CSV is the order of all_cells. usage: make_cells.py CONFIG.xml OUT.csv '{"cells.order": "AAABBB..."}'
Runs without @cells.order write nothing."""
import sys, json
cfg, out, extra = sys.argv[1], sys.argv[2], json.loads(sys.argv[3])
if "cells.order" in extra:
    open(out, "w", newline="\n").write("x,y,z,type\n" + "".join(f"0,0,0,{t}\n" for t in extra["cells.order"]))   # LF: CRLF makes PhysiCell reject every type
