# Cross-feeding (C. beijerinckii + M. barkeri): results synthesis, reduced-size series

Status 2026-11-03. Six of the eleven planned runs are done (group A of `runs_def.csv`); the dx 0.8 and dt 0.0005 runs (CF-DX dx 0.8, CF-DT dt 0.0005, CF-DXDT) are not run, so the series has no finest step yet and every comparison below is against CF-REF only. All runs are on the pinned code (`review-runs-v1`, commit 7f2b1d0a, MN5-built binary 1f60748dc4ee2e45; `provenance.csv`). Analysis: `scripts/analyse_cf.py` -> `results/cf_series.csv` (per run, snapshot, species), `results/cf_cap.csv` (cap per run, snapshot, species, substrate), `results/cf_runs.csv` (one row per run). Figure: `results/fig_crossfeeding_ref.png/.pdf` (manuscript Fig. 6 layout from CF-REF; `results/fig_crossfeeding_preview.*` is the earlier partial 840 min debug run, CF-FIG, which timed out at 2 h and is superseded).

## Design

Reduced size: box +-60 um (0.36 of the area of the 200 x 200 um paper box), half cluster of 111 cells (62 C. beijerinckii, 49 M. barkeri; same density), 2D (one layer, dz 1.25 um), dx = dy 1.25 um, dt_diffusion = intracellular_dt 0.001 min, 1800 min (30 h), snapshots every 30 min, Dirichlet on all faces (H2 = 0 at the box edge, so the box edge is a sink), ledger on, 56 threads, 2 runs per node, 5 h 12-17 min per run (dt 0.002: 2 h 29 min). Seeds 51 and 52 (CF-SEED), a same-seed repeat (CF-REP), dx 1.6 um (CF-DX) and dt 0.002 (CF-DT) are done. Criteria fixed before the runs: qualitative (both species grow; C. beijerinckii takes glucose and releases acetate, H2 and CO2; M. barkeri takes acetate and H2; M. barkeri growth follows C. beijerinckii growth) and quantitative (cell number per species at 1800 min within max(10%, 2 x seed spread of CF-SEED) of the finest step of the series).

## Reference run (CF-REF, seed 50)

C. beijerinckii: 62 to 3,945 cells (first new cell at 210 min), growth rate 0.164 /h at 3 h falling to 0.134 /h at 30 h as glucose falls (box mean 1.0 to 0.8 mM); uptake of glucose, release of acetate, H2 and CO2. M. barkeri: 49 to 65 cells (first new cell at 1590 min), growth rate rising from 0.005 /h at 3 h to 0.039 /h at 30 h as the H2 in its voxels rises (0.007 to 0.32 mM) and acetate to 0.25 mM; uptake of H2, acetate and CO2, release of methane. M. barkeri does not double within 30 h in the reduced system. In the figure the M. barkeri cluster is pushed to the edge of the expanding C. beijerinckii cluster at 30 h. Ledger: cumulative residual of the cell terms is at most 1.2e-16 of the gross flux for every substrate in every run.

## Findings

1. **Qualitative criteria are met in all six runs.** Exchange signs (median flux of the living cells at 1800 min) are as predicted in every run, M. barkeri grows in every run (60-65 cells, growth rate 0.038-0.039 /h), and its growth rate rises with the H2 supplied by C. beijerinckii (first new M. barkeri cell 1470-1680 min, 21-24 h after the first C. beijerinckii division at 210 min). The "sequential growth" pattern is therefore present, with M. barkeri growth starting late and slowly.

2. **Seed spread and repeat.** Cells at 1800 min, C. beijerinckii 3,945 / 4,004 / 3,987 (REF, s51, s52) and 3,935 (REP): seed spread 1.7%, REF vs REP 0.25%. M. barkeri 65 / 60 / 61 and 64 (REP): spread 8%, REF vs REP 1.6% (one cell is 1.6%). The threshold of the quantitative criterion is therefore 10% for C. beijerinckii and 16% for M. barkeri (2 x 8%), with the caveat that M. barkeri counts are 60-65 cells.

3. **dx 1.6 um and dt 0.002 against CF-REF (provisional, no finest step yet).**

   | run | C.b. cells 1800 | vs REF | M.b. cells 1800 | vs REF | C.b. growth 1800 (1/h) | M.b. first new cell (min) | H2 in cells' voxels, median (mM) |
   |---|---|---|---|---|---|---|---|
   | REF dx 1.25 dt 0.001 | 3,945 | - | 65 | - | 0.1339 | 1590 | 0.36 |
   | dx 1.6 dt 0.001 | 4,050 | +2.7% | 63 | -3.1% | 0.1363 | 1470 | 0.31 |
   | dx 1.25 dt 0.002 | 3,771 | -4.4% | 62 | -4.6% | 0.1249 | 1680 | 0.54 |

   Both are inside the thresholds (10% and 16%). The growth of C. beijerinckii differs by +1.8% (dx) and -6.7% (dt) in the growth rate at 30 h. The median H2 and acetate in the cells' voxels differ more than the cell numbers (H2 -14% for dx 1.6, +49% for dt 0.002; acetate -9% and +29%), so concentration observables are more sensitive to the step than the growth.

4. **Cap (R2.M1.5, R2.M1.6): M. barkeri is in the capped regime for the first half of the run.** Ratio of the recorded uptake bound to the Michaelis-Menten bound, median over living M. barkeri cells (REF): H2 0.39 / 0.78 / 0.96 at 600 / 1200 / 1800 min, with 100 / 96 / 46% of the cells below 0.95; acetate 0.59 / 1.0 / 1.0 (74 / 0 / 0% below 0.95). The cap number R of H2 has median 1.4-1.5 up to 480 min, 0.9 at 930 min and 0.5 at 1380 min, above the initial threshold of about 0.3 found in the unit test until about 20 h. The H2 and acetate uptake of M. barkeri during that time is therefore limited by the substrate present in the voxel in each step, not by the kinetics, and its growth during the first 20 h depends on dx and dt in the way the unit test predicts. C. beijerinckii is outside the cap (glucose ratio 0.988-0.993, R max 0.06) except in the dt 0.002 run (6.5% of cells below 0.95 at 1800 min). The M. barkeri growth rate at 600 min differs across runs by up to 33% (0.0115 at dt 0.002 to 0.0153), as expected where the cap acts; by 1800 min it agrees within 4% (0.0379-0.0393 /h).

5. **Reduced box against the paper box.** Not compared yet. The reduced box has the edge sink closer (H2 = 0 at +-60 um), so the H2 and acetate levels are lower than in the paper figure; the reduced reference is for pattern and for the dx / dt series, not for the figure values.

## CF2: larger box, Km 0.01, initial H2 / acetate 1 mM (2026-11-03)

CF2-BOX (box +-100 um, H2 and acetate initially 1 mM, boundaries zero except glucose and CO2 as shipped, dx 1.25, dt 0.001, 1800 min, seed 50) and CF2-KM (same, M. barkeri H2 Km 0.01 mM) on the pinned code (binary 1f60748dc4ee2e45, MN5 job 48447730, 5 h 15 min each). Analysis: `results/cf2_series.csv`, `cf2_cap.csv`, `cf2_runs.csv`; figures `results/fig_crossfeeding_CF2-BOX_s50.*`, `fig_crossfeeding_CF2-KM_s50.*` (same layout as the reference figure).

| | CF-REF (box +-60) | CF2-BOX | CF2-KM |
|---|---|---|---|
| C.b. cells 1800 min | 3,945 | 3,865 | 3,864 |
| C.b. growth rate 1800 min (1/h) | 0.134 | 0.127 | 0.127 |
| M.b. cells 1800 min | 65 | 63 | 63 |
| M.b. first new cell (min) | 1590 | 1590 | 1650 |
| M.b. growth rate 600 / 1200 / 1800 min (1/h) | 0.014 / 0.029 / 0.039 | 0.014 / 0.027 / 0.036 | 0.014 / 0.027 / 0.040 |
| box-mean H2 at 1800 min (mM) | 0.305 | 0.187 | 0.186 |
| box-mean acetate at 1800 min (mM) | 0.221 | 0.160 | 0.160 |
| M.b. H2 bound / Michaelis-Menten bound, median, 600 / 1200 / 1800 min | 0.39 / 0.78 / 0.96 | 0.36 / 0.76 / 0.93 | 0.15 / 0.44 / 0.98 |
| M.b. cells below 0.95 of the H2 bound, 600 min | 100% | 100% | 100% |

Findings:

1. **None of the three changes makes M. barkeri grow faster or changes the shape of the curves.** M. barkeri still does not double in 30 h (63 cells), its first new cell is at 1590-1650 min, and the curves of biomass, concentrations and net exports are smooth and monotonic, without the changes of concavity of the original figure. Lowering Km changes the cap ratio (0.15 against 0.36 at 600 min: the Michaelis-Menten bound is larger, the voxel limits more) but not the growth.
2. **The initial 1 mM of H2 and acetate does not matter.** The box mean is 1.0 mM at 0 and 0.003 mM at 30 min in both runs: with zero boundaries the 1 mM is removed within the first snapshot interval (diffusion to the edge and uptake) and the concentrations afterwards are set by C. beijerinckii production. The first 30 min are the only time M. barkeri sees mM concentrations.
3. **The larger box lowers the H2 and acetate levels** (0.19 against 0.31 mM and 0.16 against 0.22 mM at 30 h), because the cluster is the same size and the field is spread over a larger area, and C. beijerinckii grows 2% less. The cap acts on M. barkeri during the first 20 h as before (cap ratio 0.36 at 600 min, 81% of cells below 0.95 of the bound at 1800 min against 46% in CF-REF).
4. Ledger residuals are at most 2.4e-13 of the gross flux for every substrate in the three runs (largest: SO4 in CF2-BOX; glucose, H2, acetate and methane 1e-16 to 1e-18).

Conclusion: the slow, smooth M. barkeri growth is not explained by the box size, the Km or the initial concentrations. The remaining differences from the original figure are the old code (dead cap, dt 0.01 with unequal steps) and the sink at the boundary; whether a larger system or a no-flux H2 boundary reproduces the original behaviour is not tested.

## CF3: convergence series on the CF2-KM design, voxels not smaller than the cell (2026-11-04)

Design: box +-100 um, 111 cells, H2 and acetate initially 1 mM, M. barkeri Km(H2) 0.01 mM, boundaries zero except glucose and CO2 as shipped, 1800 min, 2D slab dz 1.25 um, intracellular_dt = dt_diffusion, 56 threads, pinned code (binary 1f60748dc4ee2e45). The largest living cells reach 2.26 (C. beijerinckii) and 2.98 um3 (M. barkeri), so the voxel volume dx x dx x 1.25 must be at least 2.98 um3: dx 1.6 (3.2 um3) and 2.0 (5.0 um3) are inside, dx 1.25 (1.95 um3) is not and dx 0.8 is not run. Runs (MN5 jobs 48468796-48468798; 2.98-10.5 h each): CF3-REF (dx 1.6, dt 0.001), CF3-SEED_s51 (seed 51), CF3-DX (dx 2.0, dt 0.001), CF3-DT dt 0.002 and dt 0.0005 (dx 1.6; dt 0.0005 is the finest step), plus CF2-KM (dx 1.25, dt 0.001) for reference. Analysis: `results/cf3_series.csv`, `cf3_cap.csv`, `cf3_runs.csv`; figures `results/fig_crossfeeding_CF3-REF_s50.*` and `fig_crossfeeding_CF3-DT_dx1.6_dt0.0005_s50.*` (manuscript Fig. 6 layout). Criterion fixed before the runs: cell number of each species and M. barkeri growth rate at 1800 min within max(10%, 2 x seed spread) of the finest step. Seed spread (CF3-REF against CF3-SEED_s51): C. beijerinckii cells 2.0%, M. barkeri cells 8%, M. barkeri growth rate at 1800 min 1.6%; thresholds 10%, 16% and 10%.

| run | C.b. cells 1800 | vs dt 0.0005 | M.b. cells 1800 | vs dt 0.0005 | M.b. growth 1800 (1/h) | vs dt 0.0005 | M.b. growth 600 (1/h) | vs dt 0.0005 | M.b. H2 bound / MM bound, median, 600 min | M.b. first new cell (min) |
|---|---|---|---|---|---|---|---|---|---|---|
| dt 0.0005 (dx 1.6) | 3,932 | - | 74 | - | 0.0412 | - | 0.0199 | - | 0.28 | 1320 |
| REF (dx 1.6, dt 0.001) | 3,870 | -1.6% | 64 | -13.5% | 0.0400 | -2.8% | 0.0155 | -22% | 0.17 | 1470 |
| seed 51 (dx 1.6, dt 0.001) | 3,795 | -3.5% | 59 | -20.3% | 0.0394 | -4.4% | 0.0149 | -25% | 0.16 | 1500 |
| dx 2.0, dt 0.001 | 3,893 | -1.0% | 66 | -10.8% | 0.0402 | -2.3% | 0.0155 | -22% | 0.17 | 1590 |
| dx 1.6, dt 0.002 | 3,772 | -4.1% | 55 | -25.7% | 0.0363 | -11.8% | 0.0115 | -42% | 0.12 | 1680 |
| CF2-KM (dx 1.25, dt 0.001) | 3,864 | -1.7% | 63 | -14.9% | 0.0399 | -3.2% | 0.0145 | -27% | 0.15 | 1650 |

Findings:

1. **Voxel size does not matter inside the safe zone.** dx 2.0 and dx 1.6 (and, outside the zone, dx 1.25) agree: C. beijerinckii cells within 1.1%, M. barkeri cells within 3-4% (inside the seed spread of 8%), M. barkeri growth rate at 1800 min within 1%.
2. **C. beijerinckii is converged in dt** (cells within 4.1% of the finest step at every dt, growth rate within 5%; glucose bound within 1-3% of the Michaelis-Menten bound).
3. **M. barkeri is not converged in dt, and the dependence is the cap.** Against dt 0.0005, the M. barkeri cell number at dt 0.001 is 13.5% lower (inside the 16% threshold) and at dt 0.002 25.7% lower (outside); the growth rate at 1800 min is 2.8% and 11.8% lower (the second outside the 10% threshold). At 600 min the growth rate is 22% (dt 0.001) and 42% (dt 0.002) lower, and the steps do not shrink with dt (0.0115, 0.0155, 0.0199 1/h at dt 0.002, 0.001, 0.0005), so no limit is reached at dt 0.0005 for the first 10-20 h. The ratio of the H2 bound to the Michaelis-Menten bound rises as dt falls (0.12, 0.17, 0.28 at 600 min) and 100% of the M. barkeri cells are below 0.95 of the bound at 600 min in every run: during the first 20 h M. barkeri is limited by the H2 in its voxel in each step (cap number above 0.3), and what it takes up per step depends on dt. The growth rate at 1800 min converges (differences 3.6e-3 and 1.2e-3 1/h between consecutive steps).
4. **Cumulative cell numbers:** M. barkeri reaches 55, 64 and 74 cells at 1800 min for dt 0.002, 0.001 and 0.0005 (seed spread of the dt 0.001 setting: 64 and 59). A finer step (dt 0.00025) was not run; it would cost about 21 h per run on 56 threads (more than the 13 h limit used here).
5. **Qualitative criteria are met in all runs** (both species grow, exchange signs as expected, M. barkeri grows after C. beijerinckii and its growth rate rises with the H2 it receives). The ledger residual is at most 1.8e-15 of the gross flux for H2, acetate, glucose and CO2 in the six runs.

What this means for the figure and the text: the C. beijerinckii results and the qualitative pattern are converged; the M. barkeri growth during the first 20 h and its final cell number depend on dt where the cap acts (the figure from dt 0.0005 shows the highest growth). The reference for the manuscript figure should be the dt 0.0005 run (finest step), and the limitation (M. barkeri growth and cell number are cap-limited and not converged in dt at dt 0.0005; the growth rate at 1800 min is) must be stated. Conclusions about the absolute time of M. barkeri growth should not be drawn from this system.

## Limits and open items

- No finest step: the dx 0.8 (dt 0.001), dt 0.0005 (dx 1.25) and dx 0.8 / 1.6 at dt 0.0005 runs are needed to say whether M. barkeri growth converges; dx 0.8 at dt 0.0005 costs about 5 x CF-REF.
- 30 h is short for M. barkeri (no doubling): the cell-number criterion rests on 60-65 cells, the growth rate and the time of the first new cell (30 min snapshot resolution) are better observables.
- Median fluxes are over cells with a nonzero flux; the exchange sign of M. barkeri CO2 is uptake in this model.
- Concentration observables are medians over the voxels of living cells and depend on the voxel size.
- CF-FIG (112 threads, 2 h debug queue) timed out at about 840 min; it has a different build directory and binary sha (7f4412de56929182, same commit).

## FVA replay (2026-11-06; `common/fva_replay.py`)

Method: for a sample of living cells of one snapshot the linear problem of the cell is rebuilt (SBML model of its cell type, lower bound of every recorded exchange set to the `_lb` column, objective bound `max_growth_rate` not binding), solved with Gurobi (gurobipy 13.0.3 through cobra 0.32.1), and flux variability analysis is run on the recorded exchanges at fraction_of_optimum 0.99999 and 0.9999 (GLPK gives the same range limits to 2e-5, checked on all 29 tables before they were redone with Gurobi). Checks: the replay objective equals the recorded `growth_rate` (median 1.2e-14, max 1.5e-10, relative); infeasible replays are the cells with zero recorded growth; the recorded flux lies inside its FVA range in 100% of the cells. A range is called undetermined when it exceeds 1e-2 mmol/gDW/h; `width_ratio` is the median of range(0.9999) / range(0.99999): 10 for a unique optimum (the range is the tolerance), about 1 for a true degeneracy. The replay is independent of the Clp solver of the run. Tables in `results/fva/` (`*_cells.csv` one row per sampled cell, fraction and exchange; `*_summary.csv`).

**CF-REF_s50 (30 h; 200 living C. beijerinckii and 65 living M. barkeri, each type sampled separately).** Acetate, CO2, H2, NH4+, Pi, SO4, glucose (CB) and CH4 (MB) are unique in all 265 cells (100% determined, widest range 3.2e-4 mmol/gDW/h, width_ratio 10); the objective equals the recorded growth rate (median 1.2e-14, max 1.5e-10, relative), the recorded flux is in the range in 100% of the cells. Same result as the earlier replay (H2 and acetate partitioning not a tie-break artefact). Only the dx 1.25 um / dt 0.002 reference was replayed, not the finer runs.
