"""Phenotypic phase plane of the MCF7 core model (Supplementary Figure, MCF7 phase plane; S6.2): growth rate as a function of the lactate
secretion flux and of the oxygen uptake flux, both fixed, over 0-0.7 and 0-3.0 mmol/gDW/h on a 150 x 150 grid (infeasible = 0 growth).
Left: heatmap of the grid; right: surface of the same grid.  Model: models/MCF7_Zielinski_2017_core.sbml (PhysiCelldFBA commit b9a337e6, md5 4ca64020...).
Constraints as in the notebook scripts/reduced_cost_analysis.ipynb of PhysiCelldFBA (section Phenotype Phase Space Analysis): DM_gudac_c_ upper bound 0,
biomass 0-0.1, ATP dump closed, ATPM lower bound 1.07, alanine, citrate and proline secretion closed; maximum growth 0.0268 1/h.
The previous figure drew the heatmap of this grid with the axis limits of an unrelated earlier cell (O2 -0.4 to -0.1, lactate 0-0.56); the data are the same.
usage: make_figure_ppp.py MODEL.sbml OUT_PREFIX      (about 3 min; cobrapy, one LP per grid point)"""
import sys, logging, numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from cobra.io import read_sbml_model
logging.disable(logging.CRITICAL)
model_file, out = sys.argv[1], sys.argv[2]
m = read_sbml_model(model_file, f_replace={'F_REACTION': lambda x: x})
m.reactions.R_DM_gudac_c_.upper_bound = 0.0
m.reactions.R_biomass_reaction.bounds = (0.0, 0.1)
m.reactions.R_DM_atp_c_.bounds = (0.0, 0.0)
m.reactions.R_ATPM.lower_bound = 1.07
for r in ("R_EX_ala_L_e", "R_EX_cit_e", "R_EX_pro_L_e"): m.reactions.get_by_id(r).upper_bound = 0.0
n = 150; lac_vals, o2_vals = np.linspace(0.0, 0.7, n), np.linspace(0.0, 3.0, n); LAC, O2 = np.meshgrid(lac_vals, o2_vals)
lac_rxn, o2_rxn = m.reactions.R_EX_lac_L_e, m.reactions.R_EX_o2_e; mu = np.zeros((n, n))
for i in range(n):
    for j in range(n):
        lac_rxn.bounds = (LAC[i, j], LAC[i, j]); o2_rxn.bounds = (-O2[i, j], -O2[i, j]); s = m.optimize(raise_error=False)
        if s.status == "optimal": mu[i, j] = s.fluxes["R_biomass_reaction"]
print("maximum growth rate %.5f 1/h; %d of %d grid points infeasible or without growth" % (mu.max(), (mu == 0).sum(), mu.size))
np.savetxt(out + "_grid.csv", np.column_stack([LAC.ravel(), O2.ravel(), mu.ravel()]), delimiter=",", header="lactate_secretion_flux,o2_uptake_flux,growth_rate_1_per_h", comments="", fmt="%.6g")
fig = plt.figure(figsize=(12, 5), dpi=300)
ax1 = fig.add_subplot(1, 2, 1); ax1.set_box_aspect(0.75)
im = ax1.imshow(mu, extent=[lac_vals.min(), lac_vals.max(), o2_vals.min(), o2_vals.max()], origin="lower", aspect="auto", cmap="Greens")
ax1.set_xlabel("Lactate secretion flux"); ax1.set_ylabel("O2 uptake flux"); ax1.grid()
cb = fig.colorbar(im, ax=ax1, shrink=0.5, pad=0.02); cb.set_label("Growth rate")
ax2 = fig.add_subplot(1, 2, 2, projection="3d")
ax2.plot_surface(LAC, O2, mu, cmap="Greens", linewidth=0.3, edgecolors="k", antialiased=True)
ax2.set_xlabel("Lactate secretion flux"); ax2.set_ylabel("O2 uptake flux"); ax2.set_zlabel("Growth rate"); ax2.view_init(elev=35, azim=-25)
fig.tight_layout(); fig.savefig(out + ".png"); fig.savefig(out + ".pdf"); print("wrote", out)
