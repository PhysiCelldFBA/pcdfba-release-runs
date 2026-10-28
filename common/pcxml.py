"""Small helpers for reading PhysiCell output."""
def labels(root):
    """Names of the columns of the *_cells.mat matrix of one snapshot, from its output*.xml root (vector labels expand to name[0], name[1], ...)."""
    lab = []
    for L in root.find(".//cell_population/custom/simplified_data").find("labels").findall("label"):
        k = int(L.get("size")); n = L.text.strip(); lab += [n] if k == 1 else [f"{n}[{j}]" for j in range(k)]
    return lab
