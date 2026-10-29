# Unit test: results synthesis

Status 2026-10-29. All 70 runs of `runs_def.csv` were run on the pinned code (`review-runs-v1`, commit 7f2b1d0a; binary f397aa14e87ace45 built from it) and all finished (rc 0); `provenance.csv` has one row per run. They were run **locally** (8-core machine, single-thread runs 6 at a time), not on MN5; the 16-thread run of UT-THREADS was oversubscribed on 8 cores, so its timing means nothing (its result does). The full printout of the checks is `results/checks_local_pin.txt`, the CONV-BOX comparison is `results/ut_convbox.csv`, BASE and THREADS are in `results/ut_base_threads.csv`.
Not done: UT-DOMAIN (needs the addon hotfix), UT-TRIM, UT-CAPRULE on these runs, UT-ELEM for the four other models.

## Results

- **dx/dt invariance with the cap inactive (check A, 45 comparisons against the plateau dx 4, dt 0.001, limit 0.5%):** 40 pass. The five failures are all dx 0.5 um at dt 0.01 and 0.02 (volume 12% and 49% from the plateau, glucose taken 16% and 63%; at dt 0.02 two of four cells remain). dx 1, 2 and 4 pass at every dt from 0.001 to 0.02 (volume within 0.054%, glucose taken within 0.041%). dx 0.5 is kept as the failure case.
- **Cap active (glucose D = 0, 1 mM; check B):** 31 of 32 pass; the voxel of the cell is emptied to within 0.035%. One failure: dx 1, dt 0.001, lowest glucose -3.5e-4 mM against a limit of -2.7e-5 mM (2 mu_max dt c0).
- **Ledger residual (check C):** 48 inactive and 47 capped checks pass; one failure: dx 0.5, dt 0.001, oxygen residual 2.3e-4 of gross against a limit of 2.7e-5 (the next largest is 4.6e-5, within its limit).
- **Carbon balance (check D):** 16 of 16 pass (largest difference 0.05%, limit 0.5%). **Elemental balance (E):** e_coli_core, 0 of 74 internal reactions unbalanced.
- **UT-CONV-BOX (8 cells growing to 32 in a 16 um box, 50 mM glucose, 26 mM left at 120 min):** compared with the plateau at 30, 60, 90 and 120 min, volume agrees within 0.03% for every dx >= 1 and every dt, and for dx 0.5 at dt 0.001 and 0.005 (0.046%); cell counts are identical. The six comparisons outside the 5% criterion (of 48 at 30, 60 and 90 min) are all dx 0.5 with dt 0.01 and 0.02: volume -9% and -20% at 30 min, -23% and -37% at 60 min, -24% and -51% at 90 min, with cells lost at 60 and 90 min.
- **Voxel sharing (1, 4, 8 cells, mixed):** all pass; glucose stays above -8e-6 mM; the 4-cell and mixed cases draw identical glucose (relative difference 0), so cells that cannot take a substrate do not dilute the others.
- **UT-ORDER (12 cells in one voxel, 5 creation orders, 1 and 4 threads):** all agree with the identity order at 1 thread to at most 1.1e-15 (criterion 1e-12); the lowest concentration is -3.4e-6 mM (solver tolerance).
- **UT-THREADS (1, 4, 16 threads):** identical in every printed digit (volume 45.0398, 32 cells, mean growth rate 0.732744).
- **UT-BASE (5 identical single-thread runs, default test):** identical results. The default test uses up its glucose (about 3e-61 at the end), so its end-state growth rate is 0.

## Reading

The time step and the voxel size do not change the single-cell and shared-voxel results as long as the voxel is not too small for the step: at dx 0.5 um the cap acts at dt >= 0.01 and the cells starve or are lost. For the colony, dx >= 1 um was safe here at every dt tested, but the colony's own dependence on dt comes from the cap at its density (`../colony/RESULTS.md`).

## Open

1. Check B at dx 1, dt 0.001 and check C at dx 0.5, dt 0.001 exceed their limits (numbers above). To decide: keep the limits and report these as exceptions, or look at the cause.
2. UT-CAPRULE (cap number R against the deviation) is not computed on these runs.
3. UT-DOMAIN, UT-TRIM, UT-ELEM for four models.
