"""Mass-balance cross-check of a closed-box consortium run (CF3-NF): carbon balance from the conservation ledger and the dry mass of the cells, and the stoichiometric bound of the methane yield per glucose from the two metabolic models.
usage: mass_balance_check.py RUNDIR [OUT_PREFIX]
Ledger units: 1 mM x um^3 = 1 amol (1e-18 mol). Cumulative net exchange of the cells at the last row: glucose (taken by C. beijerinckii), methane (M. barkeri), CO2, H2 and acetate (net of both species); in the closed box (H2, acetate: no flux boundary) the net production equals the content of the box (H2_F, acetate_F).
Carbon in = 6 x glucose taken; carbon in the products = 2 x acetate + CO2 + CH4; the rest is the carbon of the biomass made, compared with the dry mass gained by the living cells (V (1 - fluid fraction) x 1.04e-12 g/um^3) at a carbon fraction of 0.40-0.52 g/gDW (assumed range; the biomass reactions of the models have no complete formulas).
Methane: compared with the complete-conversion limit of 3 CH4 per glucose (Buswell) and with what the H2 and acetate left in the box could still give."""
import sys, os, glob, numpy as np, pandas as pd, scipy.io as sio, xml.etree.ElementTree as ET
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..", "common")); from pcxml import labels
d = sys.argv[1].rstrip("/"); out = sys.argv[2] if len(sys.argv) > 2 else None
L = pd.read_csv(f"{d}/conservation_ledger.csv"); r = L.iloc[-1]; amol = lambda k: float(r[k])
G, CH4, CO2, H2, AC = -amol("glucose_cells_dF"), amol("methane_cells_dF"), amol("CO2_cells_dF"), amol("H2_cells_dF"), amol("acetate_cells_dF")
print(f"t = {r.t_min:.0f} min: glucose taken {G:.4g} amol, CH4 {CH4:.4g}, CO2 net {CO2:.4g}, H2 net {H2:.4g} (box content {amol('H2_F'):.4g}), acetate net {AC:.4g} (box content {amol('acetate_F'):.4g})")
Cin = 6 * G; Cprod = 2 * AC + CO2 + CH4; Cbio = Cin - Cprod
print(f"carbon in {Cin:.4g} amol C; in acetate {2*AC:.4g}, CO2 {CO2:.4g}, CH4 {CH4:.4g}; carbon of the biomass made by difference {Cbio:.4g} amol C = {Cbio*12e-18*1e9:.3f} ng C")
def dry(xf):
    root = ET.parse(xf).getroot(); lab = labels(root); m = sio.loadmat(xf.replace(".xml", "_cells.mat"))["cells"]; g = lambda k: m[lab.index(k)]
    alive = g("dead") == 0; return float((g("total_volume") * (1 - g("fluid_fraction")) * 1.04e-12)[alive].sum()) * 1e9  # ng
fs = sorted(glob.glob(f"{d}/output[0-9]*[0-9].xml")); d0, d1 = dry(fs[0]), dry(fs[-1]); print(f"dry mass of the living cells: {d0:.3f} ng at the start, {d1:.3f} ng at the end, gained {d1-d0:.3f} ng")
for fc in (0.40, 0.45, 0.50, 0.52): print(f"  carbon fraction {fc}: biomass carbon gained {(d1-d0)*fc:.3f} ng C ({(d1-d0)*fc*1e-9/12/1e-18:.4g} amol C), against {Cbio:.4g} by difference: ratio {(d1-d0)*fc*1e-9/12/1e-18/Cbio:.2f}")
print(f"methane per glucose taken: {CH4/G:.3f} mol/mol; H2 accumulated per glucose {H2/G:.3f}; acetate accumulated per glucose {AC/G:.3f}")
# Methane: observed against the stoichiometric limit of the complete conversion of glucose to methane (C6H12O6 -> 3 CH4 + 3 CO2, 24 electrons / 8 per CH4) and the methane that the H2 and acetate left in the box could still give (4 H2 + CO2 -> CH4; acetate -> CH4 + CO2).
buswell = 3.0; potential = CH4 / G + AC / G + H2 / G / 4
print(f"observed methane per glucose {CH4/G:.3f}; methane that the acetate ({AC/G:.3f} per glucose) and the H2 ({H2/G:.3f}) left in the box could still give: {AC/G:.3f} + {H2/G/4:.3f}; observed + potential {potential:.3f} against the limit {buswell:.1f} (complete conversion, no biomass)")
print(f"carbon of the biomass made per glucose: {Cbio/G:.2f} C of 6 ({100*Cbio/Cin:.0f}% of the carbon taken)")
# note: a community LP of the two models (C. beijerinckii glucose -1, M. barkeri takes at most what C. beijerinckii releases) gives 4.0 CH4 per glucose, above the limit of 3, because the models are not constrained by thermodynamics and have H2O and H+ open; it is not used as a bound.
if out:
    pd.DataFrame([dict(run=os.path.basename(d), t_min=float(r.t_min), glucose_amol=G, CH4_amol=CH4, CO2_net_amol=CO2, H2_net_amol=H2, acetate_net_amol=AC, C_in=Cin, C_products=Cprod, C_biomass_by_difference=Cbio, dry_gain_ng=d1-d0, CH4_per_glucose=CH4/G, potential_CH4_per_glucose=potential, buswell_limit=buswell)]).to_csv(out + "_mass_balance.csv", index=False)
