# Mahadevan, Edwards & Doyle 2002, Biophys J 83:1331-1340 (the reference of the acetate-switch use case)

Notes taken from the PDF (`Mahadevan_2002.pdf`, a local copy in this folder, ignored by git: publisher copyright). Page numbers are the journal's.

## System (p. 1333, Eq. 5; p. 1335, Eq. 6b)
- Well-mixed batch bioreactor, 10 h, E. coli on glucose; acetate secreted first and re-used later.
- dGlc/dt = A_Glc v X; dAc/dt = A_Ac v X; dO2/dt = A_O2 v X + kLa (0.21 - O2); dX/dt = mu X.
- Gas-phase oxygen constant at 0.21 mM; kLa = 7.5 1/h (12.5 1/h in the sensitivity run).
- Initial: glucose 10.8 mM, acetate 0.4 mM, O2 0.21 mM, X0 = 0.001 g/L (biomass 1000-fold over the batch).
- Glucose uptake: Michaelis-Menten, 10 Glc / (Km + Glc) mmol/gDW/h, Km = 0.015 mM (Wong et al. 1997).
- O2 uptake <= 15 mmol/gDW/h; transport proportional to the concentration difference. Acetate transport rapid.
- Rate-of-change constraints on the fluxes (vdot_max = 0.1, 0.3, 0.3, 0.1 for the four pathways); only the SOA (static optimisation, LP per interval, 10000 intervals) and DOA formulations.
- Network: 54 metabolites, 85 reactions reduced to four pathways (v1 aerobic with O2 35, v2 aerobic, v3 glucose -> acetate + biomass, v4 acetate -> biomass).

## Results to reproduce (qualitatively; Figs 4-5, p. 1336 and text)
- Phase 1, 0 - 4.6 h: glucose and O2 uptake bounds limit growth (aerobic).
- Phase 2, 4.6 - 6.9 h: O2 near zero, transport-limited (kLa); glucose still consumed, acetate produced (peak near 4 mM at about 7 h in the figures).
- Glucose nearly exhausted at 6.9 h; then acetate utilisation, growth linear and limited by O2 transport; acetate back to zero at 9.5 h (8 h for kLa = 12.5).
- Biomass reaches about 1 g/L at 10 h. Data: Varma & Palsson 1994b.
- Sensitivity: glucose Vmax 11 mmol/gDW/h -> acetate re-use not active within the batch; kLa most sensitive parameter.
