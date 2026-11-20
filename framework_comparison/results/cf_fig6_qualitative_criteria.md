# Fig. 6 consortium benchmark: qualitative criteria (fixed 2026-11-19, before the COMETS and BacArena Fig. 6 runs)

Replaces CF-REF as the consortium of the benchmark (decision 2026-11-19: benchmark the configuration of manuscript Fig. 6).
PhysiCelldFBA run: `CF3-NC50_s50` (200 um box, dx 1.6 um, dt 0.001 min, 50 h, H2 and acetate start at 1 mM with 0 mM on
the boundary, glucose and CO2 1 mM, Km(H2) of M. barkeri 0.01 mM, models without the FVA caps; Fig. 6 shows it to 45 h).
COMETS 2.12.5: `comets_cfref_10um.py` with SCENARIO=fig6 (20 x 20 boxes of 10 um plus the fixed-concentration ring,
dt 0.25 min), MOTION=none, convection (k = 0.028 /s) and convection with k/4. BacArena 1.8.2:
`run_bacarena_fig6_80x80.R` (80 x 80 sites of 2.5 um, closed arena, hourly steps).

The criteria are those of `cf_qualitative_criteria.md` (CF-REF), with the times extended to the 50 h of the run:

| # | criterion | evaluated as | pass if |
|---|---|---|---|
| 1 | both species grow | biomass fold change 0-50 h per species (PhysiCelldFBA: biomass from the integrated mean growth rate) | > 1 for both species |
| 2 | *C. beijerinckii* takes up glucose and releases acetate, H2 and CO2 | sign of the species' total (PhysiCelldFBA: median over cells) exchange flux at 10, 20, 30, 40 and 50 h | glucose < 0; acetate, H2, CO2 > 0 at all five times |
| 3 | *M. barkeri* takes up H2 and acetate and releases methane | as 2 | H2, acetate < 0; methane > 0 at all five times |
| 4a | *M. barkeri* growth follows *C. beijerinckii* growth | *M. barkeri* specific growth rate at 50 h against 5 h | rate(50 h) > rate(5 h) |
| 4b | *M. barkeri* growth rises with the H2 at its own locations | Spearman rho between *M. barkeri* specific growth rate and H2 at its boxes / sites / cells, across snapshots 0.5-50 h | rho > 0 |

Definitions as in `cf_qualitative_criteria.md` (specific growth rate per tool, H2 at the species' locations, one
(rate, H2) pair per snapshot). Known difference stated in the text, not scored: H2 and acetate start at 1 mM, so in the
first hours M. barkeri grows on the initial H2 and acetate, not only on what C. beijerinckii releases; 4b covers the whole run.

Amendment 2026-11-19, after the CF-REF COMETS results and before any Fig. 6 COMETS output was read: the specific growth
rate of M. barkeri in COMETS is the population rate, ln(B(t+30 min)/B(t)) / 0.5 h of the species' total biomass, as for
BacArena, not the biomass-weighted per-box rate. Reason: under Convection 2D the per-box ratio also counts biomass pushed
between boxes (CF-REF, k/4: 0.019 1/h per box against 0.042 1/h for the population at 29.5 h; without biomass motion the
two definitions agree to 1e-4 1/h). The H2 for 4b is still the biomass-weighted mean over the species' boxes.
