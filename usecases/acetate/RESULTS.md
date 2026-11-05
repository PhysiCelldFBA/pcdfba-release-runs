# Acetate switch: results synthesis

Status 2026-10-29. All 42 planned runs (`runs_def.csv`) were run on the pinned code (`review-runs-v1`, commit 7f2b1d0a, local build, binary 13c39c4af8f3c0bb; `provenance.csv`: 42 ok), locally on an 8-core machine, 2-5 min each. Analysis: `scripts/analyse_ace.py` -> `results/ace_summary.csv` (one row per run), `results/ace_series.csv` (per snapshot), `results/ace_rel_to_ref.csv` (relative change against ACE-REF, PULSE excluded). ACE-FVA and ACE-SCALE are not run. The earlier batch (`results/regimes_final.csv`, binary db0d368a, before the pin) is superseded.

## Design

20 cells in a closed 3D box (+-18 um, dx 4, dt 0.01 min), 8 mM glucose, 5 mM O2 at t = 0, a 5 mM O2 pulse at 210 min, 600 min, 10 threads, snapshots every 5 min, ledger on.
Observables (`analyse_ace.py`): fermentation start / end = first / last time the median acetate exchange flux of the living cells is > 0.5; acetate reuse start = first time it is < -0.5 (5 min resolution); peak and floor of the box-mean acetate (floor = minimum after the peak); final biomass (g dry weight per litre of box); final cell number. The cap number R = q X dt / c and the fraction of cells whose recorded bound is below the Michaelis-Menten bound ("capped") are computed as in the colony, using only snapshots with t > 0 and a box-mean concentration > 0.1 mM (when a substrate is exhausted the voxel limit is the physical limit, not a discretisation effect).

## Reference run (ACE-REF, seed 50)

Fermentation (acetate secretion) from 75 to 185 min, acetate reuse from 215 min (5 min after the pulse), peak box-mean acetate 4.93 mM at 190 min, 80 cells (two doublings) and 0.479 g/L at the end. After the 5 mM O2 pulse the oxygen is used up (about 320 min) and acetate stops at 1.73 mM: with a 5 mM pulse the acetate floor is set by the oxygen supplied, not by the kinetics (pulses of 10-80 mM reach 0.652 mM, finding 5).

## Findings

1. **Seed, repeat and solver choice do not change the observables.** ACE-SEED (seeds 51, 52), ACE-REP (same seed repeated) and ACE-PFBA (parsimonious secondary objective) give the same regime times (75 / 185 / 215 min), cell number and biomass as ACE-REF; the largest difference in any observable is 0.1% (acetate floor, seed 52 and pFBA). The cells are placed from a table, so the seed enters only through the dynamics. REF and REP agree to the printed digits (no thread-scheduling difference visible here).

2. **Mechanics and phenotype steps do not matter.** ACE-DTM (dt_mechanics and dt_phenotype 0.001 instead of 0.1 and 6 min): all observables identical to ACE-REF.

3. **dx and dt (ACE-CONV): the regime times do not move; the magnitudes move by up to 10%.**

   | run | R max, O2 (glucose) | O2 cells capped (max, %) | regime times | peak acetate | final biomass | acetate floor |
   |---|---|---|---|---|---|---|
   | dt 0.001 | 0.02 (0.006) | 0 | same | -0.6% | 0.0% | -1.7% |
   | REF (dx 4, dt 0.01) | 0.25 (0.064) | 0 | - | - | - | - |
   | dx 6 | 0.07 (0.018) | 0 | same | -0.5% | 0.0% | -1.4% |
   | dt 0.05 | 2.5 (0.41) | 51 | same | +2.5% | +0.2% | +7.8% |
   | dx 2 | 21 (0.86) | 41 | same | +3.1% | +0.5% | +10.5% |

   R max is the maximum over cells and snapshots while the substrate is present (> 0.1 mM). The reference stays below the initial threshold of about 0.3 found in UT-CAPRULE; dt 0.05 and dx 2 exceed it, and these are the only runs with capped cells. Even there the regime times are identical and the peak acetate and biomass change by 0.2-3%; the acetate floor, the smallest number (about 1.7 mM), changes by 8-10%. dx 6 and dt 0.001 agree with the reference within 2%. All voxels (8, 64, 216 um3) are larger than the cell (a few um3), so the runs are inside the cell-size condition of the unit-test safe zone.

4. **Sensitivity (26 runs, factors 0.5 and 2, seed 50; plus the volume pair below; summary figure `results/fig_acetate_sens.png`, table `results/ace_sens_rel.csv`). Criterion fixed on 2026-10-28: all six observables (fermentation start, end, reuse start, peak acetate, final biomass, acetate floor) within 5% of ACE-REF.** 11 of 26 runs satisfy it.
   - *Within 5% on all six:* diffusion coefficients of glucose, acetate and oxygen (six runs, identical to ACE-REF), Km of glucose and oxygen (four runs), Vmax of oxygen x2.
   - *Regime times unchanged, but the amount of acetate left depends on acetate uptake:* Vmax acetate x0.5 (floor +118%, biomass -7%) and x2 (floor -33%, biomass +6%), Km acetate x0.5 / x2 (floor -7% / +12%). The pulse size x0.5 / x2 (2.5 / 10 mM O2) and the pulse time x2 leave the regime times unchanged (pulse x0.5: floor +92%, 56 cells; pulse x2: floor -62%).
   - *Timing depends:* Vmax glucose x0.5 (fermentation 145-360 min instead of 75-185), x2 (5-100 min); Vmax oxygen x0.5 (fermentation starts at 5 min); initial glucose x0.5 (ends at 110 min, biomass -22%) and x2 (ends at 285 min, biomass +56%, no acetate reuse by 600 min); number of cells x0.5 (fermentation 115-210 min, reuse from 230 min) and x2 (45-125 min); pulse time x0.5 (pulse at 105 min: no reuse by 600 min, peak acetate -39%).
   - Acetate reuse appears in 40 of 42 runs; it is absent by 600 min in two sensitivity runs (initial glucose x2, pulse at 105 min). The sequence fermentation then reuse is robust to diffusion, Km of glucose and oxygen and most kinetic changes; the timing and the acetate left depend on glucose and oxygen uptake capacity, initial glucose, cell number and the time of the pulse.

   - *Cells and volume scaled together (ACE-SENS_vol, added 2026-11-05, 28 runs in all):* 40 cells in 36 x 36 x 72 um and 9 cells in 36 x 36 x 16 um (18 um is not a multiple of dx 4; density 4.29e-4 and 4.34e-4 per um3 against 4.29e-4) give the regime times of the reference (75 / 185 / 215 min) and all six observables within 0.4% (peak acetate -0.02 / -0.12%, biomass -0.01 / +0.37%, floor -0.06 / -0.32%); 160 and 36 cells at the end against 80. So the dependence of the timing on the cell number in the fixed box is the biomass-to-substrate ratio, not the size of the system. 13 of 28 runs within 5% on all six observables; acetate reuse in 26 of 28.

5. **Pulse size (ACE-PULSE, 900 min).** Pulses of 10, 20, 40 and 80 mM O2 give identical regime times (75 / 185 / 215 min), 80 cells, 0.485 g/L and the same acetate floor (0.652 mM); the three seeds of the 10 mM pulse agree to 0.01 mM in acetate. Once the pulse is large enough to exhaust the acetate, the floor does not depend on its size; with 5 mM the oxygen runs out first (floor 1.73 mM).

6. **Cap in the reference box.** While a substrate is present (> 0.1 mM) no cell is capped in the reference, DTM, PFBA, PULSE, SEED or REP runs; capping appears only for dt 0.05 and dx 2 (finding 3). When a substrate is exhausted all cells are capped by definition (nothing left in the voxel): glucose at about 200 min and oxygen at about 320 min in the reference, with no change to the regimes.

## What can be said to the reviewers (R2.M1.5, R2.M1.6, R2.M1.8, R2.M3)

- The acetate switch (fermentation, then acetate reuse after re-oxygenation) occurs in 40 of 42 runs; its times are identical across three seeds, a same-seed repeat, pFBA, the mechanics / phenotype steps, dx 2, 4, 6 and dt 0.001-0.05 (5 min resolution).
- The magnitudes depend on the discretisation only where the uptake cap acts (R above about 0.3): up to 3% in peak acetate and biomass, up to 10% in the acetate floor.
- The quantitative timing depends on the kinetics and initial conditions (Vmax of glucose and oxygen, initial glucose, cell number, time of the pulse), as expected for a mechanistic model; diffusion and the Km of glucose and oxygen do not matter in this closed box.

## Limitations and open items

- Time resolution is 5 min, 6.7% of the fermentation start (75 min); the 5% criterion cannot resolve smaller changes in that time. "Same" for the regime times means identical snapshots.
- The box is small (20-80 cells) and closed; the conclusions about dx and dt hold at this density (R 0.003-0.25 in the reference). Denser or larger systems are not covered (see `../colony/RESULTS.md`).
- 80 cells (two doublings) is the end state of 39 of 42 runs; growth stops when glucose and oxygen are used up, so final-biomass comparisons mostly measure the amount of substrate converted.
- ACE-FVA (replay on the ACE-REF snapshots) and ACE-SCALE are not run.
- `ALL_CITATIONS.txt` in this directory is written by the binary at start-up; it is not a result.
