# Colony: results synthesis (pre-grown study)

Status 2026-10-29. Everything below comes from the 16 runs of `runs_def.csv` that have been run (one more, COL-NEC_x1_t40_s50, is listed and not run), all on the pinned code (`review-runs-v1`, commit 7f2b1d0a, MN5-built binary 226118d62b231202; see `provenance.csv`). Numbers are in `results/colony_summary.csv` (per run at 2 min, 6 min and the last snapshot), `results/pg_profiles.csv` (radial profiles) and `results/cap_check_*.csv`.
The full-length growth runs (COL-REF, COL-SEED, COL-REP; quarter colonies, dt 0.0005, 600 min) were run on 2026-11-02 and are reported in the section "Full-growth reference, seeds and repeat" at the end of the quarter-domain section. The pre-grown sections are not statements about the 600 min colony.

## Design

A disc of 87,685 cells (R0 200 um, hexagonal packing) is placed at t = 0 with 4 mM glucose and 5 mM O2 on the box faces (dx 2, dz 2, 2D), necrosis rate x10, 10 min, ledger on. This isolates the zonation (necrotic core, acetate-using ring) from the slow growth phase and costs 1 h per run instead of 10 h.
Observables: dead fraction, necrotic-core radius (outer edge of the contiguous region that is at least 50% dead, 10 um annuli), acetate ring (annuli where at least 20% of living cells take up acetate), and the fraction of living cells whose recorded uptake bound is below the Michaelis-Menten bound ("capped", ratio < 0.99). Cap number R = q X dt / c is defined in `scripts/cap_check.py`.

## Findings

1. **Seed and repeat spread are negligible.** Seeds 50, 51, 52 and the seed-50 repeat give 41.8% dead at 10 min in all four, core radius 130 um, ring 120-130 to 190 um. (The only variation is one ring edge moving by one annulus.) The cells are placed identically in this design, so this is the spread of the dynamics only.

2. **The structure depends strongly on dt at the shipped setting and converges only for dt <= 0.0005 min.**

   | dt (min) | dead % at 2 min | dead % at 6 min | core / ring at the last snapshot (um) | capped glucose / acetate (last snapshot, %) |
   |---|---|---|---|---|
   | 0.01 | 0.1 | 2.0 | none / 10-200 (10 min) | 96 / 95 |
   | 0.004 | 1.0 | 12.2 | 100 / 80-190 (10 min) | 79 / 92 |
   | 0.002 | 2.9 | 24.6 | 120 / 110-190 (10 min) | 55 / 89 |
   | 0.001 (shipped) | 6.5 | 37.3 | 130 / 120-190 (10 min) | 29 / 60 |
   | 0.0005 | 13.1 | 45.6 | 140 / 140-190 (10 min) | 2 / 2 |
   | 0.00025 | 23.0 | 48.5 | 140 / 140-190 (6 min) | 0 / 0 |
   | 0.0001 | 38.8 | - | 130 / 140-180 (2 min) | 0 / 0 |

   Between dt 0.0005 and 0.00025 the late state agrees (core 140 um in both, ring 140-190 um, dead fraction 45.6% vs 48.5% at 6 min, 6% relative). dt 0.001 is not converged (core 130 um, 41.8% dead at 10 min against 47.4% at dt 0.0005). dt 0.01 gives no necrotic core.

3. **What changes with dt is the uptake cap, not the kinetics.** The fraction of capped cells falls from about 95% (dt 0.01) to 2% (dt 0.0005) and 0% (dt <= 0.00025), and the core radius tracks it. With the cap inactive the median uptake equals the Michaelis-Menten bound, and the late structure no longer changes with dt. The cap number at dt 0.0005 is below 1 for most cells (median R 0.3 glucose, 0.7 acetate); at dt 0.001 it is 0.7 and 1.6-1.7.

4. **The time at which the core forms keeps shortening as dt falls, even with the cap inactive.** At 2 min the dead fraction is 13.1, 23.0 and 38.8% for dt 0.0005, 0.00025 and 0.0001, with no sign of levelling. dt 0.00025 is at its plateau by about 5 min and dt 0.0001 reaches 38.8% by 2 min, which suggests the same end state is reached earlier, but the dt 0.0001 run stops at 2 min and cannot confirm it. Timing of necrosis (when) is therefore not a converged quantity at any dt tested; the end-state structure (where) is.

5. **Voxel size: with the cap mostly inactive (dt 0.0005, 6 min) the structure no longer depends on dx.** dx 1, 2 and 4 give 45.2, 45.6 and 44.6% dead (within 2%), core radius 130, 140 and 130 um (one 10 um annulus apart), and the same ring edges (140-190 um). The acetate-using fraction of living cells varies more (40.4, 44.5 and 50.1%), following the acetate cap, which is slightly more active at dx 1 and 4 (7.2 and 6.4% of cells capped against 2.0% at dx 2). At dt 0.001, where the cap acts in 35-73% of cells, dx 1 and 4 gave 40.2 and 37.3% dead at 10 min (core 130 and 120 um) against 41.8% (130 um) at dx 2: that comparison mixes discretisation with the cap.

5b. **A smaller colony is not a substitute.** R0 100 um (21,931 cells), dt 0.0005, 6 min gives 4.0% dead, core 20 um and ring 20-80 um, against 45.6%, 140 um and 140-190 um for R0 200 um at the same settings. The cap is not the reason (1.0% glucose and 1.7% acetate capped); the small colony has a different structure at this age. The sensitivity experiments therefore keep R0 200 um.

6. **Necrosis rate x1 instead of x10 changes the state at 10 min** (22.2% dead, core 100 um at 10 min, against 41.8% and 130 um), as expected for a death rate ten times lower. The runs do not show whether x1 reaches the same core given more time; the x10 speed-up is therefore an acceleration whose equivalence is not demonstrated.

7. **Cost.** Wall time on one 112-core node for 87,685 cells: dt 0.01 6 min, 0.004 17 min, 0.002 32 min, 0.001 57-67 min, 0.0005 80 min (all for 10 simulated min); dt 0.00025 60 min for 6 simulated min; dt 0.0001 45 min for 2 simulated min. Cost scales as 1/dt.

## What can be said to the reviewers (R2.M1.5, R2.M1.6, R2.M1.8, R2.M3)

- The colony zonation (necrotic core surrounded by an acetate-using ring) is present at every dt up to 0.004 and disappears only at dt 0.01. It is a qualitative result of the model.
- Core radius, ring inner edge and dead fraction are quantitative results that depend on dt at the shipped setting through the uptake cap, and converge once the cap is inactive (dt <= 0.0005 in this colony). The shipped setting (dt 0.001) lies 7% below that limit in core radius and about 12% in dead fraction at 10 min.
- Stochastic spread in this design is below the resolution of the observables (annulus width 10 um).
- Voxel size does not change the structure once the cap is mostly inactive (dx 1, 2, 4 at dt 0.0005: dead fraction within 2%, same ring edges, core within one annulus).
- The cap number R (unit-test experiment UT-CAPRULE) predicts where the cap acts; the colony data are the in-situ check.

## Criterion

The 5% criterion fixed on 2026-10-28 (core radius and ring fraction against COL-PG) was relaxed on 2026-10-29, after seeing the data, at the user's decision: core and ring edges are compared to the resolution of the analysis (one 10 um annulus, 7% of the core radius) and the dead fraction within 10% relative. With this criterion dx 1, 2 and 4 at dt 0.0005 agree; dt 0.001 against dt 0.0005 (core 130 against 140 um, 41.8 against 47.4% dead at 10 min) is at the limit; dt >= 0.002 is outside it. The relaxation is stated here and must be stated in the supplement.

## Limitations observed

- **Crowding and cell number.** The more cells and the denser the packing, the less stable and the more step-dependent the simulation becomes: the cap fraction grows with the biomass per voxel, the structure changes with dt up to a 10-20 fold finer step than the shipped one, and the one crash seen in this use case (dx 4, dt 0.01, mechanics step too large for the density) happened at the densest settings. The cause of that crash was not confirmed. The small colony (21,931 cells) does not reproduce the 87,685-cell colony. This is an observation across these runs, not a demonstrated mechanism; it is stated as a limit of the framework at high cell density.
- The dt that converges the colony (0.0005 min) costs about 80 min per 10 simulated minutes for 87,685 cells on one 112-core node.

## Quarter-domain validation (COL-QTR, COL-QTRG, COL-QTRL, COL-QTRN, COL-QTRB)

Purpose: the full colony costs 15-18 h per run at dt 0.001 (legacy ledger run: 18 h 10 min on one 112-core node), too much for the reference, seed and repeat runs. By symmetry only the quadrant x >= 0, y >= 0 is simulated: domain from 0, no-flux planes at x = 0 and y = 0 (BioFVM default without Dirichlet), Dirichlet only on the outer faces, cells held by `virtual_wall_at_domain_edge`. Needs the option `quarter_colony` (PhysiCelldFBA branch `feature/colony-quarter-disc`, commit 23947c49 on the pin; local, not pushed; binary c846ce99 built on MN5). It is a validation set: none of these runs replaces a COL-REF run, and the legacy full run (`runs/colony_ref_s51`, 600 min, +-500 um box) is not on the pin.

| Experiment | Design | Result |
|---|---|---|
| COL-QTR | pre-grown colony (R0 200 um), dt 0.0005, 6 min, quarter vs COL-DT_dt0.0005 | structure matches (criterion of this file) |
| COL-QTRG | shipped start (radius 10 um), +-150 um box, 240 min, full vs quarter, seeds 50, 51 | quarter has about 1.5% higher growth rate, about 4% larger radius; axis and diagonal radii agree |
| COL-QTRL | as COL-QTRG to 340 min (22k cells) | cell count quarter x4 / full quadrant 1.09-1.10 from 40 min on, both seeds, constant in time (no divergence); R99 +3.7% (99.9 vs 96.3 um at 340 min); axis/diagonal radii within 1%; no dead cells in either setup |
| COL-QTRN | quarter only, +-250 um box, 500 min, seeds 50, 51 | first dead cells at 360 min (0.1%), 1.8-2.0% at 380 min, as in the legacy run (2.4% at 360); core 90 um and dead fraction 18% at 500 min; but cells and radius 22-50% above the legacy run at 360-480 min: the box was smaller than the legacy one (supply 250 vs 500 um away), a confounded comparison |
| COL-QTRB | quarter only, +-500 um (legacy box), 600 min, seeds 50, 51 | matches the legacy full run (table below) |

COL-QTRB against the legacy full run (`colony_ref_s51`; cells = quarter x 4; both seeds agree to 0.3%):

| t (min) | cells quarter / legacy | dead % | R99 (um) | core (um) | ring outer edge (um) |
|---|---|---|---|---|---|
| 360 | 24.6k / 23.2k | 3.1 / 2.4 | 108 / 106 | 0 / 10 | 80 / 80 |
| 480 | 49.8k / 47.9k | 29.4 / 28.5 | 160 / 158 | 80 / 80 | 150 / 140 |
| 600 | 96.0k / 93.6k | 40.0 / 39.7 | 213.6 / 211.9 | 140 / 130 | 200 / 200 |

At 600 min the quarter has 2.6% more cells, 0.8% larger radius, the same ring outer edge and a core one 10 um annulus larger (the relaxed criterion allows one annulus). Radii along the axes and along the diagonal agree within 0.5% at every time.

**Artifact at the symmetry planes.** Inside the necrotic core the legacy colony is 100% dead; in the quarter, 65 of 5,149 cells within 100 um are alive, all of them within 6 um of the x = 0 or y = 0 plane (some at x down to -0.4 um, i.e. slightly outside the domain), healthy (growth rate 0.79, glucose and oxygen uptake at full rate). The cause was not investigated. Consequences: the acetate ring appears to start at 0 um instead of 140 um in the quarter, so the ring inner edge is not measurable there, and zone statistics should exclude cells within about 6 um of the planes. Effect on the dead fraction: small (1.3% of the cells in the core).

Limits of the validation: the cell-count offset of +9-10% in the +-150 um box (COL-QTRL) is not seen in the legacy box at 600 min (+2.6%); the two differ in box and in the run used for reference, and cell counts at a single snapshot have a spread of several percent from synchronous division waves (QTRL quarter 5916 vs QTRN 5537 cells at 340 min with the same seed; thread count and box differ), so the offset is not established beyond that scatter. The comparison with the legacy run is indicative only.

## Full-growth reference, seeds and repeat (COL-REF, COL-SEED, COL-REP)

Four quarter-colony runs from the shipped start (colony_radius 10 um) to 600 min: COL-REF_s50, COL-SEED_s51, COL-SEED_s52, COL-REP_s50 (same seed as REF). Legacy box (+-500 um, quarter 0..500 um), dx 2, dt_diffusion = intracellular_dt 0.0005, dt_mechanics = dt_phenotype 0.001, 4 mM glucose / 5 mM O2, ledger on, 112 threads, 4 nodes; 5.0-5.6 h per run (job 48389037), rc 0, 31 snapshots each. Code: the pin (7f2b1d0a) plus the unmerged option `quarter_colony` (branch `feature/colony-quarter-disc`, commit 23947c49, binary c846ce99, MN5-built; the branch is not pushed). Numbers: `results/ref_summary.csv`, `results/ref_zones.csv`, `results/ref_profiles.csv`, `results/cap_check_COL-{REF_s50,SEED_s51,SEED_s52,REP_s50}.csv`. Cell counts are the quarter count x 4. Zone statistics (core, ring) exclude nothing by hand: the alive cells at the symmetry planes (see above) affect only the ring inner edge, which is not measurable (reported 0).

| run | cells x4 at 360 / 480 / 600 min | dead % at 600 | R99 at 600 (um) | core at 360 / 480 / 600 (um) | ring outer edge 360 / 480 / 600 (um) |
|---|---|---|---|---|---|
| REF s50 | 24,796 / 50,388 / 98,008 | 46.2 | 215.2 | 30 / 90 / 150 | 80 / 150 / 200 |
| SEED s51 | 24,784 / 50,380 / 98,020 | 46.2 | 215.3 | 30 / 90 / 150 | 80 / 150 / 200 |
| SEED s52 | 24,788 / 50,340 / 98,124 | 45.9 | 215.1 | 30 / 90 / 150 | 80 / 150 / 200 |
| REP s50 | 24,800 / 50,436 / 98,188 | 46.2 | 215.2 | 30 / 90 / 150 | 80 / 150 / 200 |

1. **Seed spread and repeat meet the criterion.** Cell number differs by 0.1% between seeds and 0.2% between REF and REP at 600 min; the dead fraction by at most 0.3 percentage points (0.6% relative); R99 by 0.2 um; core and ring edges are identical at the 10 um annulus resolution at all three times. REF and REP are not bit-identical (98,008 vs 98,188 cells): thread scheduling contributes noise of the same size as the seed spread.
2. **Cap fraction at dt 0.0005 is small.** Fraction of cells whose recorded bound is below the Michaelis-Menten bound (REF, 360 / 480 / 600 min): glucose 4.3 / 6.4 / 7.4%, acetate 4.1 / 14.2 / 8.4%, oxygen about 0.01%; median R 0.13-0.42 (glucose) and 0.49-0.87 (acetate). In COL-QTRB at dt 0.001 the fractions were 24-36% (glucose) and 43-66% (acetate), median R 0.3-0.9 and 1.0-1.7. The other three runs agree with REF within 1.5 percentage points. Acetate is still partly capped at its median R (0.7-0.9, above the initial threshold of about 0.3 only in part of the cells): the cap acts on a minority.
3. **Ledger (R2.M1.7).** The cumulative residual of the cell sink / source terms is at most 5e-11 of the gross flux for oxygen and glucose and 2e-14 for acetate and CO2 in all four runs; the largest single-step residual is 4.8e-9 mmol (oxygen) against a cumulative gross flux of 2.4e8 mmol at 600 min. Diffusion-side residuals were not analysed.
4. **Against COL-QTRB (dt 0.001) and the legacy full run (colony_ref_s51, seed 51, not on the pin).**

   | | cells x4 at 360 / 480 / 600 | dead % at 360 / 480 / 600 | R99 at 600 (um) | core at 360 / 480 / 600 (um) | ring outer edge at 600 (um) |
   |---|---|---|---|---|---|
   | legacy full | 23,160 / 47,938 / 93,632 (full count) | 2.4 / 28.5 / 39.7 | 211.9 | 10 / 80 / 130 | 200 |
   | QTRB s50 (dt 0.001) | 24,636 / 49,812 / 96,048 | 3.1 / 29.3 / 40.0 | 213.6 | 0 / 80 / 140 | 200 |
   | REF s50 (dt 0.0005) | 24,796 / 50,388 / 98,008 | 6.7 / 36.0 / 46.2 | 215.2 | 30 / 90 / 150 | 200 |

   REF has +4.7% cells and +1.6% radius against the legacy run at 600 min, and the ring outer edge is identical. The dead fraction is 6.5 percentage points higher (16% relative) and the core 20 um (two annuli) larger; the same offset is present at 360 and 480 min (dead 6.7 vs 2.4%, 36.0 vs 28.5%). **Against the legacy run, the relaxed criterion (core and ring edges within one 10 um annulus, dead fraction within 10% relative) is not met for the core and the dead fraction.** QTRB, at the legacy dt, reproduces the legacy run (dead 40.0 vs 39.7%), so the difference is the time step (0.0005 vs 0.001) and not the quarter method. This agrees with the pre-grown finding that structure converges only for dt <= 0.0005: dt 0.001 results, including the legacy colony, understate the necrotic fraction by about 6 percentage points at 600 min.

Limits: convergence at dt 0.0005 at 600 min is addressed in the next section (COL-REFDT); the ring inner edge is not measurable in the quarter runs; the legacy run is a different code build (not the pin), seed 51, one run; seeds 50-52 are three samples.

## Time-step check at 600 min (COL-REFDT)

COL-REFDT_dt0.00025_s50: COL-REF with dt_diffusion = intracellular_dt 0.00025 (half), everything else identical (seed 50, dt_mechanics = dt_phenotype 0.001, quarter colony, binary c846ce99, 112 threads, one node, job 48410453, 7 h 39 min, rc 0, 31 snapshots). Numbers: `results/refdt_zones.csv`, `results/refdt_profiles.csv`, `results/cap_check_COL-REFDT_dt0.00025_s50.csv`.

| | cells x4 at 360 / 480 / 600 | dead % at 360 / 480 / 600 | R99 at 600 (um) | core at 360 / 480 / 600 (um) | ring outer edge at 360 / 480 / 600 (um) |
|---|---|---|---|---|---|
| REF s50 (dt 0.0005) | 24,796 / 50,388 / 98,008 | 6.7 / 36.0 / 46.2 | 215.2 | 30 / 90 / 150 | 80 / 150 / 200 |
| REFDT s50 (dt 0.00025) | 24,876 / 50,700 / 99,260 | 10.2 / 40.1 / 48.0 | 216.2 | 30 / 100 / 150 | 80 / 150 / 200 |

1. **At 600 min the relaxed criterion is met.** Cells +1.3%, R99 +1 um, core and ring edges identical at 600 min (core 150 um, ring 200 um), dead fraction 48.0 vs 46.2% (+3.9% relative, +1.8 percentage points; the seed spread of REF is 0.3 points, so the difference is real but inside the 10% criterion).
2. **Earlier times are not converged.** The dead fraction at 480 min is 40.1 vs 36.0% (+11% relative, core 100 vs 90 um, one annulus), and at 360 min 10.2 vs 6.7%: necrosis starts earlier as dt falls, as in the pre-grown study. The structure at 600 min is converged for the criterion; the timing of necrosis is not.
3. **The dead fraction at 600 min against dt:** 39.7% (legacy, dt 0.001), 46.2% (0.0005), 48.0% (0.00025). The steps are 6.5 and 1.8 points, a ratio of about 3.6, consistent with convergence from below. A Richardson-type extrapolation would put the limit near 48.5-49% (assumes first order; not tested).
4. **Cap fraction at dt 0.00025 is negligible:** glucose 0.008%, acetate 0.008%, oxygen 0.008% at 600 min (median R 0.09 glucose, 0.38 acetate), against 7.4 / 8.4 / 0.01% at dt 0.0005. The residual difference between dt 0.0005 and 0.00025 is therefore partly the cap acting on a minority of cells at dt 0.0005.
5. **Ledger:** cumulative residual of the cell terms at 600 min: oxygen 1.3e-6 mmol against 2.4e8 mmol gross (5e-15), as for COL-REF.

Limits: one run, seed 50; REF and REFDT differ by thread scheduling noise too (REF vs REP: 0.2% in cells); dt 0.0001 was not run at 600 min.

## Time-step series with the necrosis increase scaled with dt (COL-NECS)

The necrosis rate of a cell grows by `death_rate_increase` at each dFBA update, so the ramp per minute is increase/dt and the dead fraction can depend on dt through the model. COL-NECS repeats the dt series with increase = 1.67e-5 x dt/0.001 (the ramp is then 1.67e-5 per min2 at every dt): dt 0.001 (increase 1.67e-5), 0.0005 (8.35e-6) and 0.00025 (4.175e-6), otherwise COL-REF (seed 50, quarter colony, 600 min, 112 threads, one node each, jobs 48490287-48490289, 31 snapshots). Numbers: `results/necs_summary.csv` (computed with `pg_analyse.analyse`).

| at 600 min | dt 0.001 | dt 0.0005 | dt 0.00025 |
|---|---|---|---|
| dead %, increase scaled with dt | 40.1 | 45.0 | 47.0 |
| dead %, increase not scaled (COL-QTRB / COL-REF / COL-REFDT) | 39.7 | 46.2 | 48.0 |
| core radius (um) | 140 | 150 | 150 |
| ring outer edge (um) | 200 | 200 | 200 |
| R99 (um) | 214.6 | 215.8 | 217.0 |
| quarter cells | 23,968 | 24,517 | 24,829 |

1. **The scaling changes the dead fraction by at most 1.2 points** (1.2 at dt 0.0005, 1.0 at 0.00025). The dependence of the dead fraction on dt is therefore in the state of the cells and not in the necrosis ramp.
2. **Convergence:** the dead fraction changes by 4.9 points (11% relative) from 0.001 to 0.0005 and by 2.0 points (4.4%) from 0.0005 to 0.00025. The relaxed criterion is met at dt 0.0005 (and, unscaled, 3.8%). The steps shrink by a factor of 2.5, so the limit dt -> 0 is about 48-49% on the scaled series (first-order assumption, not tested). Core and ring edges agree within one annulus, cells within 1.3% between 0.0005 and 0.00025.
3. **Earlier times** are not converged, as before (dead % at 360 min: 3.1 / 5.6 / 8.1 scaled; 6.7 / 10.2 unscaled at 0.0005 / 0.00025).
4. **Decision:** the reference stays COL-REF_s50 (dt 0.0005, 1.67e-5 per update); Table S4 states the increase per dFBA update; S7.3 gets the paragraph C5 (`Colony_params_edits.tex`).

Limits: one run per dt, seed 50; the cause of the remaining dt dependence (cap number, trigger timing or other) was not tested.

## Open items and replan

| Item | Status | Proposal |
|---|---|---|
| COL-PG (seeds + repeat) | done | none |
| COL-DT | done for dt 0.0001-0.01; dt 0.00025 stops at 6 min and dt 0.0001 at 2 min (2 h debug limit) | optional: dt 0.0001 to 6 min in the standard queue (about 135 min), to confirm the plateau |
| COL-DX | done at dt 0.001 (10 min) and at dt 0.0005 (6 min) | none |
| COL-NEC | done at dt 0.001, 10 min only | in `runs_def.csv`, not run: x1 at dt 0.001 to 40 min in the standard queue (about 4 h), to see whether x1 reaches the x10 core; otherwise state x10 as an acceleration |
| COL-NECS | done 2026-11-04 (dt 0.001, 0.0005, 0.00025 with the increase scaled with dt) | reference unchanged; see the section above |
| COL-REF / SEED / REP | done 2026-11-02 (quarter colonies, dt 0.0005, 600 min) | COL-REFDT (dt 0.00025) done 2026-11-03: structure at 600 min within the relaxed criterion; timing of necrosis (360-480 min) still shifts earlier with smaller dt |
| COL-R0 | done: R0 100 um is not equivalent | none; sensitivity runs keep R0 200 um |
| COL-BND, COL-KIN, COL-PFBA | not run | run at R0 200 um, dt 0.0005, 6 min (about 50-60 min per run on one node, so the 32 kinetic sets fit one 32-node debug job at R0 200 um instead of the planned 350 um); whether 6 min is long enough for the kinetic variations to show is untested |
| COL-FVA | not run | the pre-grown snapshots on the pin already carry the `_lb` columns; the replay can start on COL-PG_s50 and COL-DT_dt0.0005_s50 now |
| COL-SCALE | wall times above | threads scaling not run |
| Ledger (R2.M1.7) | recorded in all runs | residuals not yet analysed |
