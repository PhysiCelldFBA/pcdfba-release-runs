# Acetate switch: results synthesis

Status 2026-10-29. **No definitive run on the pinned code (`review-runs-v1`, 7f2b1d0a) exists for this use case.** The runs in `runs/` are an earlier batch (binary db0d368a012b0ae1, before the pin) with the old design (dx 4, dt 0.01, 20 cells in the ±18 um box, seeds 50-52, pulse at 210 min). `results/regimes_final.csv` is their analysis. The numbers below are indicative; they do not replace `ACE-REF` and the experiments that depend on it.

## Observed in the earlier batch (9 runs: seeds 50, 51, 52 at the standard pulse, plus pulse variants)

- **Regime times are identical for the three seeds** at the 5 min snapshot resolution: fermentation from 75 to 185 min, acetate reuse starting at 215 min (the pulse is at 210 min). 80 cells at the end in every run.
- **Peak acetate varies by 3%** between seeds (7.29, 7.49, 7.51 mM); the other recorded extremes (glucose and oxygen exchange minima) are identical.
- **The regime times do not depend on the size of the pulse** in the runs tried (runs named o210, o220, o240, o280, plus the standard pulse; seed 50, three seeds for the standard and o210): same 75/185/215 min.
- The discretisation was dx 4, dt 0.01 throughout; no dx or dt comparison has been made on the pin.

## What the colony study suggests for this use case

The colony shows that dt dependence enters through the uptake cap (cap number R = q X dt / c): the structure converges once R is below 1 for most cells (`../colony/RESULTS.md`). The cap number for this box (20 cells in 36^3 um3, dx 4) has not been computed or measured. The assumption that the acetate switch is insensitive to dx and dt must be tested (ACE-CONV, ACE-DTM) with the cap fraction recorded.

## Not done

ACE-REF, SEED, REP, CONV, DTM, PFBA, FVA, PULSE and SENS on the pin (configs are generated from the tables: `runs_def.csv`, `overrides_def.csv`); the cap-fraction measurement for this box; ACE-SCALE.
