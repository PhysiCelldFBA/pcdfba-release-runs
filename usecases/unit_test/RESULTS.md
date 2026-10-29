# Unit test: results synthesis

Status 2026-10-29. **The definitive runs on the pinned code (`review-runs-v1`, 7f2b1d0a) have not been run.** What exists is a first trial of the whole table (70 runs, 70 finished) made with an earlier build (binary a1bad17583ed09cc, before the pin), kept in `runs_trial1_pre_pin/` and checked with the scripts of `dfba_unit_test`; the full printout is `results/trial1_checks.txt`. The statements below are indicative until the runs are repeated on the pin, and they are not for the manuscript yet.

## Observed in the trial

- **dx/dt invariance with the cap inactive (check A, 45 comparisons against the finest run, limit 0.5%):** 40 pass. All five failures are dx 0.5 um at dt 0.01 and 0.02 (volume 12% and 49% from the reference, one run loses half its cells). dx 1, 2 and 4 pass at every dt from 0.001 to 0.02 (volume within 0.054%, glucose taken within 0.041%). dx 0.5 is kept as the failure case: the likely cause is the cap (a 0.125 um3 voxel holds little substrate per step), but the cap number was not computed for these runs.
- **Cap active (glucose D = 0, check B):** 31 of 32 checks pass; the cell's voxel is emptied to within 0.035%. The one failure is dx 1, dt 0.001: lowest glucose -3.5e-4 mM against a limit of -2.7e-5 mM.
- **Ledger residual (check C):** 47 capped and 48 inactive checks pass; one failure (dx 0.5, dt 0.001, oxygen residual 2.3e-4 of gross against 2.7e-5).
- **Carbon balance (check D):** 16 of 16 pass. **Elemental balance (E):** e_coli_core, 74 internal reactions, 0 unbalanced.
- **Voxel sharing (1, 4, 8 cells, mixed):** glucose stays non-negative to -8e-6 mM (solver tolerance); the 4-cell and mixed cases draw identical glucose (relative difference 0), so cells that cannot take a substrate do not dilute the others.
- **UT-BASE (5 identical single-thread runs):** identical to the printed digits. **UT-THREADS:** volume agrees to 0.001% between 1, 4 and 16 threads.
- **UT-CONV-BOX (16 cells, shared voxels, growth and division):** volume agrees with the finest run to 0.023% for dx >= 2 at every dt, and for dx 1 up to dt 0.01; dx 1 at dt 0.02 differs by 13.7%, and dx 0.5 loses half its cells at dt >= 0.005.

## Problems the trial found (to fix before the definitive runs)

1. **UT-ORDER failed**: the update-order check cannot find `R_EX_o2_e` in the trial binary. The label is missing from that build's output; rerun on the pin.
2. **UT-CONV-BOX exhausts the glucose** in the box before 120 min, so the remaining glucose is zero and the glucose comparison is void (the percentages in `results/ut_convbox.csv` for glucose are meaningless). Use more glucose or a shorter time so substrate remains at the end.
3. The trial binary is not in `binaries.tsv`; the definitive runs will be.
4. Check C has no tolerance margin at dx 0.5 and check B has none at dx 1, dt 0.001; decide whether these limits stay as written or the failure is reported.

## Not done

UT-DOMAIN (needs the addon hotfix), UT-TRIM, UT-CAPRULE on the pin (`scripts/cap_rule.py` exists), UT-ELEM for the other four models.

## Link to the other use cases

The cap number R = q X dt / c relates the failures here to the colony dt dependence: see `../colony/RESULTS.md`.
