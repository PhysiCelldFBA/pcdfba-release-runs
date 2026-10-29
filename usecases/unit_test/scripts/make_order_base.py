"""Derive base/order_base.xml (the update-order test: four cell types A B C N in one voxel, cells from a CSV) from the shipped
voxel-sharing 'mixed' configuration, the way gen_update_order_configs.py does. Rerun after every repin of code/.
types: A 1.3 um3, glucose Vmax 8 | B 2.6 um3, 8 | C 1.3 um3, 24 | N 1.3 um3, 0 (cannot take glucose)."""
import re, os
UC = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
src = open(f"{UC}/config/shipped/voxel_sharing_test/Unit_test_voxel_sharing_mixed.xml").read()
m1 = re.search(r'( *)<cell_definition name="ecoli1" ID="0">.*?</cell_definition>\n', src, re.S)
m2 = re.search(r'( *)<cell_definition name="ecoli2" ID="1">.*?</cell_definition>\n', src, re.S)
assert m1 and m2 and m2.start() == m1.end(), "unexpected layout of the mixed config"
types = {"A": (1.3, "8.0"), "B": (2.6, "8.0"), "C": (1.3, "24.0"), "N": (1.3, "0.0")}; defs = ""
for i, (t, (vol, vmax)) in enumerate(types.items()):
    b = m1.group(0).replace('name="ecoli1" ID="0"', f'name="{t}" ID="{i}"').replace('<cell_adhesion_affinity name="ecoli1">', f'<cell_adhesion_affinity name="{t}">')
    b = re.sub(r'(<volume>\s*<total units="micron\^3">)[\d.]+', lambda m: m.group(1) + str(vol), b)
    b = re.sub(r'(<exchange substrate="glucose">.*?<Vmax units="mmol/gDW cell/h">)[\d.]+', lambda m: m.group(1) + vmax, b, flags=re.S); defs += b
cfg = src[:m1.start()] + defs + src[m2.end():]
assert '<cell_positions type="csv" enabled="false">' in cfg
cfg = cfg.replace('<cell_positions type="csv" enabled="false">', '<cell_positions type="csv" enabled="true">')
os.makedirs(f"{UC}/base", exist_ok=True); open(f"{UC}/base/order_base.xml", "w").write(cfg)
print("base/order_base.xml:", cfg.count("<cell_definition "), "cell definitions")
