# Unit test: results synthesis

Status 2026-10-29. All 70 runs of `runs_def.csv` were run on the pinned code (`review-runs-v1`, commit 7f2b1d0a; binary f397aa14e87ace45 built from it) and all finished (rc 0); `provenance.csv` has one row per run. They were run **locally** (8-core machine, single-thread runs 6 at a time), not on MN5; the 16-thread run of UT-THREADS was oversubscribed on 8 cores, so its timing means nothing (its result does). The full printout of the checks is `results/checks_local_pin.txt`, the CONV-BOX comparison is `results/ut_convbox.csv`, BASE and THREADS are in `results/ut_base_threads.csv`.
Not done: UT-DOMAIN (needs the addon hotfix). UT-TRIM done 2026-11-06 (section at the end). UT-ELEM done 2026-11-06 (section at the end).

## Results

- **dx/dt invariance with the cap inactive (check A, 45 comparisons against the finest run printed in `results/checks_local_pin.txt` (dx 0.5, dt 0.001), limit 0.5%; the cap-number analysis and UT-CONV-BOX use the plateau dx 4, dt 0.001, where the cap cannot act, and give the same pass/fail pattern: dx >= 1 within 0.072% in volume, 0.046% in glucose taken):** 40 pass. The five failures are all dx 0.5 um at dt 0.01 and 0.02 (volume 12% and 49% from the plateau, glucose taken 16% and 63%; at dt 0.02 two of four cells remain). dx 1, 2 and 4 pass at every dt from 0.001 to 0.02 (volume within 0.054%, glucose taken within 0.041%). dx 0.5 is kept as the failure case.
- **Cap active (glucose D = 0, 1 mM; check B):** 31 of 32 pass; the voxel of the cell is emptied to within 0.035%. One failure: dx 1, dt 0.001, lowest glucose -3.5e-4 mM against a limit of -2.7e-5 mM (2 mu_max dt c0).
- **Ledger residual (check C):** 48 inactive and 47 capped checks pass; one failure: dx 0.5, dt 0.001, oxygen residual 2.3e-4 of gross against a limit of 2.7e-5 (the next largest is 4.6e-5, within its limit).
- **Carbon balance (check D):** 16 of 16 pass (largest difference 0.05%, limit 0.5%). **Elemental balance (E):** e_coli_core, 0 of 74 internal reactions unbalanced.
- **UT-CONV-BOX (8 cells growing to 32 in a 16 um box, 50 mM glucose, 26 mM left at 120 min):** compared with the plateau at 30, 60, 90 and 120 min, volume agrees within 0.03% for every dx >= 1 and every dt, and for dx 0.5 at dt 0.001 and 0.005 (0.046%); cell counts are identical. The six comparisons outside the 5% criterion (of 48 at 30, 60 and 90 min) are all dx 0.5 with dt 0.01 and 0.02: volume -9% and -20% at 30 min, -23% and -37% at 60 min, -24% and -51% at 90 min, with cells lost at 60 and 90 min.
- **Voxel sharing (1, 4, 8 cells, mixed):** all pass; glucose stays above -8e-6 mM; the 4-cell and mixed cases draw identical glucose (relative difference 0), so cells that cannot take a substrate do not dilute the others.
- **UT-ORDER (12 cells in one voxel, 5 creation orders, 1 and 4 threads):** all agree with the identity order at 1 thread to at most 1.1e-15 (criterion 1e-12); the lowest concentration is -3.4e-6 mM (solver tolerance).
- **UT-THREADS (1, 4, 16 threads):** identical in every printed digit (volume 45.0398, 32 cells, mean growth rate 0.732744).
- **UT-BASE (5 identical single-thread runs, default test):** identical results. The default test uses up its glucose (about 3e-61 at the end), so its end-state growth rate is 0.

## Cap number (UT-CAPRULE)

R = n q m_dw dt / (c V_vox) = q X dt / c is the fraction of a voxel's substrate that the cells in it ask for in one step (n cells in the voxel at the start; q at the run's own glucose; `scripts/cap_rule.py`, `results/ut_vs_plateau.csv`, `results/fig_cap_rule.png`). For these cells R = 4.5 n dt / dx^3 at 10 mM glucose.

| Experiment | R passing | R failing |
|---|---|---|
| UT-DXDT-A (1 cell, 10 mM; limit 0.5% at 120 min) | up to 0.18 (deviation <= 0.12%) | 0.36 (-12.3%), 0.72 (-48.7%) |
| UT-CONV-BOX (8 cells, 50 mM; limit 5% at 30-90 min) | up to 0.29 (<= 0.044%) | 0.58 (24.4%), 1.15 (51.3%) |

Both experiments change between 0.29 (largest R that passes) and 0.36 (smallest that fails), so the initial threshold is **R* of about 0.3**, with a margin of roughly 0.2-0.5 because the grid has few points around it. Below R* the deviation from the plateau is under 0.15%, above it at least 12%. R grows with the biomass in the voxel, so an initial R of 0.3 reaches R = 1, where the cap starts to act on a cell, once the biomass has grown a few times. In the colony the cap acts on a cell at its current R (`../colony/RESULTS.md`).

## Safe zone

Two conditions, both necessary:
1. **Voxel at least as large as the largest cell.** A cell is a point source in the voxel that holds its centre, which is only consistent if it fits in a voxel. Voxel volume (dx dy dz) >= the largest volume a living cell reaches before dividing (about 2.3-2.9 um3 for the cells here; take it from the output, not from the configured starting volume). With several cell types the largest sets the voxel.
2. **Initial R below about 0.3** for the substrates the cells take up (above that, the uptake cap and not the kinetics sets the dynamics; see the colony).

Outside the first condition: dx 0.5 um (voxel 0.06 times the largest cell) fails check A at dt >= 0.01 and loses cells in CONV-BOX; check C exceeds its limit at dt 0.001. dx 1 um (0.5 times the cell) passes A and CONV-BOX but check B overdraws one voxel by 0.035% of its content at dt 0.001. These two exceptions are taken as outside the validity of the framework, not as a reason to change the limits, which are kept as written.

| Use case | Voxel (um) | Voxel volume (um3) | Configured cell volume (um3) | Largest living cell measured (um3) | Voxel / largest cell |
|---|---|---|---|---|---|
| Unit test, dx 2 (default) | 2 x 2 x 2 | 8 | 1.3 | 2.29 (UT-BASE) | 3.5 |
| Unit test, dx 1 | 1 x 1 x 1 | 1 | 1.3 | 1.95 (CONV-BOX) | 0.5 |
| Unit test, dx 0.5 | 0.5^3 | 0.125 | 1.3 | 1.95 | 0.06 |
| Colony (shipped and pre-grown) | 2 x 2 x 2 | 8 | 1.3 | 1.49 (pre-grown, 10 min) | 5.4 |
| Acetate switch, shipped | 2 x 2 x 2 | 8 | 1.3 | 2.31 (earlier batch) | 3.5 |
| Acetate switch, ACE-REF design | 4 x 4 x 4 | 64 | 1.3 | 2.31 (earlier batch) | 28 |
| Cross-feeding | 1.25^3 | 1.95 | 1.3 / 1.8 | 2.2 / 2.9 (earlier run) | **0.67 (outside)** |
| Tumour | 15^3 | 3375 | 2494 | not measured | 1.35 (configured) |
| Motility | 5 x 5 x 10 | 250 | 8 | 8.0 (earlier run) | 31 |

Sources: voxel and configured volumes from the shipped configurations. Measured volumes come from the pinned unit-test and colony runs where they exist; the acetate, cross-feeding and motility values come from earlier runs (read only for cell size) and must be re-read from the definitive runs. The colony value is from the 10 min pre-grown design; a full-length run is likely to show larger cells (2.3 um3 in the unit test), still inside the zone. The tumour has no measurement yet. Condition 2 is not yet computed for the use cases other than the colony.

## Reading

The time step and the voxel size do not change the single-cell and shared-voxel results as long as the voxel is large enough for the cells and the initial cap number is below 0.3. The colony's own dependence on dt comes from the cap at its density (`../colony/RESULTS.md`).


## Analytical endpoint of the mass-conversion test (UT-BASE)

The endpoint is the initial biomass plus (mu / q_glc) times the glucose in the box, with mu and q_glc from FBA of the run's own model at the glucose (8 mmol/gDW/h) and oxygen (18) Vmax of its configuration (`scripts/make_figure.py`): mu 0.7333 1/h, q_glc 8.0, yield 0.0917 gDW per mmol glucose (0.509 gDW/g), 0.338 pg initial biomass, so 4.093 pg; the run ends at 4.093 pg (8 cells, glucose exhausted), a difference below 0.01%. The shipped `Ecoli_core.xml` has ATPM at 0; with ATPM at 8.39 the same calculation gives 3.874 pg. **The earlier figure and notebook used a hard-coded yield of 0.5197 gDW/g (`dfba_analysis.ipynb`, "matches FBA model"), which this model does not give (0.5092): it overstates the endpoint (4.170 pg) and shows a spurious -1.85% error.** The values in the submitted text and S4 Table 1 (4.08 / 4.09 and 3.87 / 3.88 pg) agree with the FBA-derived 4.093 and 3.874 pg to within 0.02 pg. The with-ATPM variant and the iML1515 repetition are not in the pinned runs (UT-ATPM, UT-IML1515 would have to be added to report them from this repository).

### Other models (UT-ATPM, UT-IML1515; 2026-11-03)

Same test (8 cells, closed 16 um box, dx 2, dt 0.001, pinned code), analytical endpoint from FBA of each run's own model (`scripts/endpoint.py`, `results/ut_endpoints.csv`):

| run | model | mu (1/h) | yield (gDW/g glc) | simulated (pg) | analytical (pg) | difference |
|---|---|---|---|---|---|---|
| UT-BASE | E. coli core (ATPM lb 0) | 0.7333 | 0.5092 | 4.093 | 4.093 | -0.002% |
| UT-ATPM | E. coli core, ATPM lb 8.39 | 0.6906 | 0.4796 | 3.869 | 3.874 | -0.135% |
| UT-IML1515 | iML1515 | 0.6964 | 0.4836 | 3.900 | 3.904 | -0.085% |

The endpoint agrees within 0.14% for the three models (expected 3.874 and 3.904 pg were fixed before the runs). With maintenance the yield drops by 6%, as FBA predicts, and the simulation follows it.

## Open

1. UT-DOMAIN (needs the addon hotfix).
2. Condition 2 (R) for the other use cases; the cross-feeding setting at dx 1.25 is outside condition 1 and its main-figure setting needs a decision.
3. Voxel/cell volumes of the acetate, cross-feeding, motility and tumour use cases to be re-read from definitive runs.

## UT-ELEM: elemental balance of the models (2026-11-06; `scripts/elem_balance.py`, `results/elem_balance.csv`, `results/elem_unbalanced.csv`)

cobra 0.32 `check_mass_balance` (elements and charge) on every reaction that is not an exchange and not a biomass reaction ("internal"); a reaction is checkable when all its metabolites have a formula. Files are the ones the runs read (sha256 in the table).

| Model | Internal reactions | Checkable | Balanced | Unbalanced | No formula |
|---|---|---|---|---|---|
| e_coli_core (unit test, acetate, colony, motility; also with ATPM 8.39) | 74 | 74 | 74 | 0 | 0 |
| iML1515 | 2373 | 2373 | 2372 | 1 (PUACGAMS, a generic R group in the formula) | 0 |
| C. beijerinckii iCB925 (trimmed) | 484 | 462 | 462 | 0 | 22 |
| M. barkeri iMG746 (trimmed) | 517 | 517 | 517 | 0 | 0 |
| MCF7 Zielinski 2017 core (trimmed) | 346 | 346 | 340 | 6 | 0 |

C. beijerinckii: 22 reactions involve 19 pseudo-metabolites without formula (DNA, RNA, protein, peptidoglycan, teichoic acid, carbohydrate, phospholipid classes, menaquinones) and cannot be checked. MCF7: ARTFR11, ARTPLM1 and ARTPLM2 carry generic acyl-chain symbols in their formulas (a formula artefact, not a missing atom balance), DHCR242r, PSERT and PSP_L differ by one H (protonation state). Not a check of the biomass reactions (their composition is lumped by construction). Biomass reactions (every reaction with "biomass" in its id) are excluded from the table.

## UT-TRIM: trimmed against full models (2026-11-06; `scripts/trim_check.py`, `results/trim_check.csv`, `results/trim_removed.csv`; full models in `models_full/`, origin in `models_full/ORIGIN.md`)

Medium: every exchange of the trimmed model gives its bounds to the same reaction of the full model (case-insensitive id, iCB925 uses `ex_` against `EX_`), all other exchanges of the full model are closed. Gurobi; FBA optimum and FVA (fraction 0.99999) of the exchanges of the trimmed model. Variant "full bounds": the trimmed topology with the internal bounds of the full model, to separate the removal of reactions from the tightening of bounds.

| Model | Reactions full / trimmed | Removed reactions that carry flux in the full model | Internal bounds changed | Optimum full / trimmed (full bounds) / shipped trimmed | Largest difference of an exchange limit: topology only / shipped |
|---|---|---|---|---|---|
| C. beijerinckii iCB925 | 938 / 501 | 0 of 437 | 485 | 0.332977 / 0.332977 / 0.332977 | 1e-6 / 1e-6 |
| M. barkeri iMG746 | 816 / 551 | 0 of 265 | 518 (314 upper bounds from 1000 to below 1) | 0.051696 / 0.051696 / **0.043597 (-15.7%)** | 0 / **20.0** |
| MCF7 Zielinski 2017 core | 381 / 374 | 0 of 9 | 1 | 0.025915 (biomass_NCI60) / 0.003776 / 0.003776 (biomass_reaction) | 0.28 / 0.28 (not comparable) |

Findings. (1) Removing reactions did not change the behaviour of any model: the reactions removed are all blocked in the medium (range 0), and the trimmed topology with the full bounds gives the full optimum and the full exchange ranges (iCB925, iMG746; 1e-6). (2) The shipped trimmed iCB925 also has different internal bounds (485 reactions) but gives the same optimum and exchange ranges. (3) **The shipped trimmed iMG746 has 518 internal bounds changed, 314 of them capped (upper bound 1000 to below 1), which lowers the optimum in the trimming medium by 15.7% and changes the exchange limits of the species (largest differences: H+ 20.5 against 0.46, H2O -6.5 against 8.2, H2S -10 against -0.03, SO3 10.0 against 0.0, CO2 -13.5 against -5.9, CH4 14.4 against 7.1 mmol/gDW/h)**: a restricted model, not only a trimmed one. The cross-feeding result therefore uses the capped M. barkeri; whether the caps bind in the medium of the runs (acetate and H2 supplied by C. beijerinckii) was not tested. The text of the paper should say "trimmed and with bounds from the preparation notebook" (to confirm with the model's author: the cap values look like flux ranges from an earlier FVA). (4) The MCF7 trimmed file is not a trimmed copy of the full file: it has another biomass reaction (biomass_reaction against biomass_NCI60), an added ATPM, 27 exchanges with other bounds (glucose -0.5 against -5 etc.), so the comparison is not valid and the optimum differs (0.0038 against 0.0259); it is the calibrated tumour model of S5, not a trimming of the full one. (5) The trimmed iMG746 of the repository (sha256 f324fd6c2a40) differs from the one the runs read (d219cf19bd06) but gives the same optimum and ranges; which one the manuscript cites is to be settled. Not done: iML1515 (the runs use the full model), and the effect of the iMG746 caps in the cross-feeding medium.
