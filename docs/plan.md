# Definitive experiments: plan for a clean rerun (draft, 2026-10-27)

Why: a bug was caught and fixed (parallel races and locks, merged in `development` as `a2519f06`), configs then changed several times, and some "seed" runs turned out to be same-seed repeats or wrote into another run's folder. Everything that goes into the paper is regenerated from one commit, with one provenance record. Nothing is deleted; old results are moved to `dfba_runs/_archive_pre_rerun/` (decision O1).

## 0. Provenance rules (apply to every run)

1. **One build.** A tagged commit of `development` (after PR #18 and any other open fixes), built once; the commit hash goes in every run's record. Colony-size runs use the same hash on MN5 and locally.
2. **One generator per use case.** Run configs are produced by a script from the shipped sample config (`make_run_config*.py`), never edited by hand. The script changes only: output folder, threads, snapshot interval, `max_time`, recording of `_lb`, and the factor under test (dx, dt, seed, a parameter). The script and its arguments are stored with the run.
3. **Seed = both values.** `<options>` seed and user `random_seed` are set together (the sample's `SeedRandom` reads the user parameter). A repeat with the same seed is named `_rep`, never `_s`.
4. **Unique output folder, checked.** The generator writes the folder name from the run name and the launcher refuses a non-empty folder. After each run, a check reads the saved `PhysiCell_settings.xml` in the output and compares seed, dx, dt, folder and commit with the run record.
5. **Run record.** One `runs_def.csv` (planned) and one `provenance.csv` (what happened: commit, config hash, seed, dx, dt, threads, machine, start/end, status, output path) per use case. Figure scripts read only run ids from these files (this is also the reproducibility index of R1.M2).
6. **Folder layout.** `dfba_runs/<use_case>/<experiment>/<run_id>/`. Experiments are named as in this document (E1 ... E6).

## 1. Experiments

### E1. Unit tests (`dfba_unit_test`, PR #18 plus merged fixes)
- Grid: dx {0.8, 1.6, 3.2} x dt {0.005, 0.01, 0.02}; inactive and capped series; ledger residual per metabolite; carbon balance; elemental balance of the model (checks A-E in S11).
- Cells per voxel: 1, 4, 8, and mixed (consumers and non-consumers), single thread and multi-thread, 10 update orders (existing voxel-sharing and update-order tests).
- Elemental balance of every model used (Table R10: e_coli_core, iML1515, C. beijerinckii, M. barkeri, MCF7 trimmed) and the trimmed-SBML before/after check (R2.M6).
- Repeat the baseline against itself (the test is not reproducible run-to-run) to set the noise floor.
- Cost: minutes. First thing to run.

### E2. Degeneracy (alternative optima)
- FVA replay at 0.99999 of the optimum for every use case (colony, acetate, cross-feeding, motility, tumour) on the final reference runs: where flux ranges are non-zero for the exchanges that set the gradients. Expected: only the tumour (and 40 min in the acetate onset) is affected.
- pFBA run for the tumour at the shipped point; FBA vs pFBA for each use case (one pair each), plus the tolerance check.
- Open question to settle: the FBA vs pFBA difference in tumour lactate (~30%) is untested. Decide whether it gets a dedicated pair (O3).
- Depends on E4 reference runs (the FVA is a replay of those).

### E3. Sensitivity and seeds
- **Seeds:** the reference config of every use case at seeds 50, 51, 52 (both seed values set), 3 independent runs each. Plus one `_rep` of seed 50 per use case for the thread-scheduling noise floor.
- **Sensitivity per use case** (design from the existing work, rerun on the final build):
  - colony: pre-grown colony, boundary glucose/O2 grid, kinetic Sobol (Km, Vmax, D);
  - tumour: vessel O2 x glucose grid (5x5), 16-point Sobol on Km, Vmax, D at two O2 levels; necrosis threshold variation;
  - cross-feeding: one-at-a-time scan of Km, Vmax, D, boundary;
  - acetate switch: one-at-a-time scan (not done so far);
  - motility: regime map over Km and boundary.
- Cost: the largest block. Use the surrogate/earlier results only to choose the points; every number in the paper comes from the new runs.

### E4. Harmonised reference runs (the main figures)
- Rule rather than identical values: **dx and dt are set per use case by the biology and the geometry, but the uptake cap must not be the strongest limit**, tested with the cap-vs-Michaelis-Menten ratio (`cap_vs_mm.py`). `intracellular_dt = dt_diffusion` everywhere.
- Current state, to be unified:

| Use case | dx (um) | dt_diffusion (min) | dx^2/dt | Note |
|---|---|---|---|---|
| Colony | 2 | 0.001 | 4000 | dt_mechanics, dt_phenotype also 0.001 |
| Acetate switch | 1.0 -> 2 | 0.01 -> 0.001 | 100 -> 4000 | decided 2026-10-27 to match the colony; test runs on 2026-10-27 |
| Cross-feeding | 0.8 -> 1.25 | 0.01 -> 0.001 | 64 -> 1562 | decided D2; cap weak at 1562, still acting at 640 (dx 0.8) |
| Motility | 5 | 0.001 | 25000 | harmonised in PR #16 |
| Tumour | 15 | 0.01 | 22500 | cell scale is tissue; dx/dt convergence untested |

- **2D slab thickness (found 2026-10-27).** In 2D PhysiCell sets the z range to +-0.5*dz (`modules/PhysiCell_settings.cpp:323`), so the medium depth per unit area is dz, and the shipped 2D configs tie dz to dx. Changing dx therefore also changes how much substrate lies under each cell. Evidence: acetate switch at dx 2 / dz 2 (test runs, partial) has 1600 cells against 800 at dx 1 / dz 1, fermentation starts at 120 min against 40, peak acetate flux 7.0 against 9.4. A dx test in 2D must either hold dz fixed (convergence in the discretisation only) or state that the slab thickness changes (a different system). To decide per use case in E4, and to report thickness in every parameter table: colony dz 2, acetate dz 1, motility dz 5 (to be checked).
- **Acetate switch: the shipped timings were set by the cap (found 2026-10-27, `acetate_switch/cap_vs_mm.py`).** Old setting (dx 1, dz 1, dt 0.01, dx^2/dt = 100): the recorded O2 uptake bound is 0.90 / 0.82 / 0.26 / 0.03 of the Michaelis-Menten rate at 40 / 60 / 100 / 150 min while median O2 at the cells is still 1.45 / 0.81 / 0.17 / 0.025 mM (Km 0.002 mM). Test run (dx 2, dz 2, dt 0.001): ratio 1.00 until 200 min, O2 exhausted by kinetics (1e-11 mM at 150 min). So the old fermentation onset (40 min) and acetate peak were cap-limited, and the use case has to be re-established at a setting where the cap is not the strongest limit: dx and dt from the cap rule, dz held at 1, then initial O2 (shipped 5 mM, far above air-saturated ~0.2 mM) and glucose (8 mM) chosen from a cheap well-mixed dFBA surrogate to give the sequence aerobic -> O2-limited overflow -> acetate re-use -> plateau, and checked with a few PhysiCell runs on MN5.
- Also harmonise: simulated time per use case (cross-feeding 1800 vs 3000 min), snapshot interval, `dt_mechanics`/`dt_phenotype` (decide whether to match the colony's 0.001 or leave at 0.1/6), recording of `_lb`.
- Convergence table (R2.M1.8): each use case at the reference and one finer/coarser dx and dt, criterion 5% fixed in advance (cross-feeding already has the series; the other four need it).
- Output feeds: Figs 2-6, S10, the summary table.

### E5. Items from the reviews that are not experiments of the above
- R2.M5 scaling study (cells vs wall time, threads) and the hardware table.
- R3.M1 benchmark against COMETS/BacArena: decide run or state as limitation (O4).
- R2.M6 archival release: tag, commit, DOI, Zenodo; models repo version.
- Fig. 4: vessel O2 value checked against a literature dissolved-O2 value (text, no run).
- Unit check of the cross-feeding Fig. 6 panel D (fmol/s vs per-gDW flux) against the original notebook.
- Cap-ratio diagnostic as addon output (D5, deferred).
- Colony boundary values for the main figure (D1).

## 2. Order of work

1. Freeze: merge PR #18, tag the build, write the generators and the run-record/check script (provenance rules). Move old results to the archive.
2. E1 (minutes), then E4 reference runs (these define the settings), with the cap-ratio check.
3. E2 on the E4 runs.
4. E3 seeds on E4 settings; then sensitivity blocks per use case (MN5 for colony/cross-feeding/tumour, local for acetate/motility).
5. Regenerate figures and tables only from the run ids in `runs_def.csv`/`provenance.csv`; then fill the letter and manuscript placeholders.

## 3. Items you may not have listed (check)

- Convergence/dx-dt for the **non-cross-feeding** use cases (part of E4, easily forgotten).
- The 5% criterion and the comparison metric fixed **before** the runs.
- Thread-scheduling noise floor per use case (`_rep`).
- Mass-conservation ledger on **use-case runs**, not only the unit test (R2.M1.7): which binary produced the earlier numbers is unknown.
- The `vox20` crash (dx 2.0, cross-feeding): cause unknown; either find it or state it.
- Reproducibility index and the hardware/wall-time table (R1.M2, R3 minor 2).
- Developers' email and the `fix/parallel-race-conditions` PR are done or user-owned; nothing to run.

## 4. Decisions (2026-10-27)

- **O1. Settled.** Archive everything and rerun everything from E1. Old results are reference only. Nothing is deleted; the archive starts once the two acetate dx 2 / dt 0.001 test runs (started 2026-10-27) have finished, because they still write into `dfba_runs/acetate_switch/run`.
- **Repository. Settled.** New repository `pcdfba-review-runs`. One pinned copy of the code as the submodule `code/` (PhysiCelldFBA, branch development, currently `a2519f06`; replaced by a tag before the definitive runs), `build/build_all.sh` builds one binary per use case from it in scratch trees, and `usecases/<uc>/` holds config, bin, scripts, runs (untracked, symlink to scratch) and results. Rules in the README; the check is `common/check_run.py` (saved settings == config rebuilt from the tables).
- **O2. Open, decided by data.** `dt_mechanics`/`dt_phenotype`: one run per use case in E4 with the shipped values and with 0.001; keep the shipped values if the observables agree within the 5% criterion.
- **O3. Settled.** The tumour FBA versus pFBA pair is computed in E2 and the result decides the wording.
- **O4. Parked.** COMETS/BacArena benchmark is the last item, after all use cases are harmonised.
- **O5. Settled.** Cross-feeding final figure: 50 h from t = 0. Seeds (3, from t = 0) and sensitivity use shorter runs (1800 min unless E4 shows a reason to change). Pre-grown cross-feeding (restart from a snapshot of the reference run) is used only for sensitivity runs, with the history-fixed caveat stated; first test in E4 that positions, volumes and the microenvironment can be reloaded and that the cells reach the reference state within a few steps (the addon has no `save_dFBA()`, so each cell's LP state is rebuilt).
- **O6. Open.** Whether all five sensitivity blocks of E3 are run in full or some are cut to one-at-a-time scans; decided after E4 fixes the reference runs.
- **Ledger. Proposed, not started.** Move the conservation ledger from `dfba_unit_test` into the addon (opt-in, off by default, `feature/ledger` from development) so every use case can write `conservation_ledger.csv`; it answers R2.M1.7 on the use-case runs and needs a per-step cost check. Must be merged before the code tag used for the definitive runs.

## 5. Acetate switch: reference and 0-D surrogate (2026-10-27)

Reference: Mahadevan, Edwards & Doyle 2002 (notes in `usecases/acetate/reference/mahadevan2002_notes.md`): well-mixed batch, glucose 10.8 mM, acetate 0.4 mM, O2 0.21 mM with gas-phase resupply kLa = 7.5 1/h, X0 = 0.001 g/L, 10 h; O2 exhausted at 4.6 h, glucose gone at 6.9 h, acetate peak then re-use, acetate back to 0 at 9.5 h, X about 1 g/L.
Surrogate: `usecases/acetate/scripts/surrogate_0d.py` (our E. coli core model, our bound rule, explicit Euler, 6 s per run) and `scan_0d.py` (`results/scan_0d.csv`).
- With the shipped kinetics (glucose Vmax 8, Km 0.02; O2 Vmax 18, Km 0.002; acetate Vmax 5, Km 1) and kLa 7.5, X0 = 0.005 g/L, glucose 8 mM: overflow starts 4.1 h, glucose gone 7.35 h, re-use from 7.3 h, acetate peak 1.85 mM (reference about 4), X 0.67 g/L. The timeline is qualitatively that of the reference.
- With the paper's Vmax (glucose 10, O2 15) our core model overflows from t = 0 (it needs about 2.2 O2 per glucose): the reference's phase 1 is lost. Keep the shipped Vmax.
- X0 sets the timing (X0 0.001 gives glucose still left at 10 h). The shipped PhysiCell setup has 200 cells (3.4e-13 gDW each: 1.3 um^3, fluid 0.75, 1.04 g/mL) in 700 x 700 x 1 um = 0.138 g/L, 28 times the X0 that fits and 138 times the reference.
- To match X0 in PhysiCell the medium per cell must increase: slab thickness dz about 28 um for 200 cells in 700 x 700 um, or fewer cells in a smaller domain. A thick slab also makes the cap weak (cap scales with the voxel volume), so the shipped dt 0.01 would be adequate.
- Missing in the PhysiCell sample: oxygen resupply. BioFVM has bulk supply (`supply_rates`, `supply_target_densities_times_supply_rates`, `bulk_supply_rate_function`); the sample does not use it. kLa = 7.5 1/h would be a supply rate of 7.5/60 per min toward 0.21 mM.

## Acetate decisions (2026-10-27)

- Batch 1 (10 runs, MN5 job 48145642, all rc=0): the pulse at dx 4 / dt 0.01 gives the three-phase switch with regime times independent of dx (2, 4, 6), dt (0.001-0.05) and seed (50-52); the continuous O2 supply gives no clean switch and is dropped (`feature/o2-supply` stays local, not to be merged). See `usecases/acetate/results/regimes_3d.csv`.
- Acetate figure: pulse, dx = dy = dz = 4, dt 0.01, 3D cube.
- Cell data renamed to the SBML ids (`R_EX_*` and `R_EX_*_lb`) in branch `feature/acetate-rex-fluxes` (worktree `pcdfba-acetate-rex`), to harmonise with the other use cases. The final acetate binary needs this branch merged and the code tag moved.
- TO DO LATER: scale the model up (more cells, larger domain) to reproduce the density maps of the original figure.
