# Benchmark against COMETS and BacArena (Supplementary Section S12, R3.M1)

Two scenarios of the manuscript reconstructed in COMETS 2.12.5 and BacArena 1.8.2, with the same biology as the
PhysiCelldFBA runs they are compared with (models, uptake kinetics, medium, boundary, initial biomass, duration):

| run_id | tool | scenario | PhysiCelldFBA run | wall time |
|---|---|---|---|---|
| BM-ACE-COMETS | COMETS 2.12.5 | acetate switch (ACE-REF) | ACE-REF_s50 | 3.4 s |
| BM-ACE-BACARENA | BacArena 1.8.2 | acetate switch (ACE-REF) | ACE-REF_s50 | 258 s |
| BM-CF-COMETS | COMETS 2.12.5 (Convection 2D) | consortium of Fig. 6 | CF3-NC50_s50 | 4.6 h |
| BM-CF-BACARENA | BacArena 1.8.2 (closed arena) | consortium of Fig. 6 | CF3-NC50_s50 | 1.6 h |

These four runs are the ones in the three figures (S23 to S25) and the three tables (S24 to S26) of S12. The trial runs that
led to them are not part of this branch; what they taught (how each scenario had to be adapted to each tool) is in S12.1 and in the
header comment of each script. All four ran on the workstation mpws (bsclife006, Ubuntu 24.04.5, 8 cores), one thread
each. Settings, inputs, end time and wall time of each run: `runs.csv`.

## Layout

    runs.csv                 one row per run: script, settings (env variables), inputs, host, end time, wall time, folder on mpws
    scripts/                 the scripts as they ran (sha256 in config/SHA256SUMS = the copies in the run folders on mpws)
      comets_ecoli_aceref_v2.py         BM-ACE-COMETS
      run_bacarena_aceref_v6.R          BM-ACE-BACARENA
      comets_cfref_10um.py, extract_comets_cf_qual.py           BM-CF-COMETS (run_comets_10um.sh: the launcher used on mpws)
      run_bacarena_fig6_80x80.R, extract_bacarena_fig6_qual.R   BM-CF-BACARENA
      rerun_mpws.sh                     reruns one of the four on mpws (below)
      make_benchmark.py                 tables and figures from data/ and the PhysiCelldFBA results of usecases/
      extract_comets_fields.py, extract_bacarena_fields.R, extract_physicelldfba_fields.py
                                        snapshot fields and cells of the three consortium runs (added 2026-11-20)
      make_microenv.py                  Supplementary Figure S25 (fields and cells, 10-30 h) from those snapshots
    config/                  inputs that are not elsewhere in the repository, and SHA256SUMS of every input
    env/                     the software environment on mpws, captured 2026-11-19
    data/<run folder>/       the small outputs of each run that make_benchmark.py reads
    results/                 benchmark_summary.csv (every number of S12), benchmark_fig6_criteria.csv, the criteria fixed
                             before the consortium runs (cf_fig6_qualitative_criteria.md) and the three figures
    runs/                    raw outputs of a rerun (not tracked)

Raw outputs are not stored in the repository or on Zenodo: the BacArena consortium state (`.rds`, 134 MB), the COMETS
`media_out` / `flux_out` / `biomass` files and the logs remain in the run folders on mpws named in `runs.csv`
(about 210 MB in all).

## Redraw the tables and the figures

From `benchmark_frameworks/`:

    python3 scripts/make_benchmark.py
    python3 scripts/make_microenv.py

It reads `data/` and, for PhysiCelldFBA, `usecases/acetate/results/ace_series.csv`, `usecases/crossfeeding/results/cf3nc50_series.csv`
and the two `provenance.csv` (all tracked), and writes `results/`. On 2026-11-19 the two PNGs it wrote were byte-identical
to `figures/supp_fig_benchmark_{aceref,crossfeeding}.png` of the manuscript, and the two CSVs to the ones S12 was written from.

### Microenvironment figure (S25, 2026-11-20)

`make_microenv.py` reads the snapshot fields extracted from the raw outputs (in `data/`, about 5 MB):

- `data/fig6_comets/comets_fig6_conv_fields.npz`: `extract_comets_fields.py` on the BM-CF-COMETS run folder. COMETS writes
  `media_out` with 1-based coordinates (the static ring at 1 and N) and `biomass` with 0-based ones; the script shifts the
  media by one. `extract_comets_cf_qual.py`, as it ran, does not; checked on 2026-11-20 that this changes none of the numbers
  of S12 (the H2 at M. barkeri is identical, acetate within 0.05%); only the COMETS box means of
  `comets_cfref_fig6_conv_media_means.csv`, which S12 does not use, include the ring and are about 11% low.
- `data/fig6_bacarena/bacarena_fig6_closed_nopool_{media,cells}.csv.gz`: `extract_bacarena_fields.R` on the saved state of
  BM-CF-BACARENA (every hourly step).
- `data/fig6_physicelldfba/physicelldfba_CF3-NC-CTL_s50_fields.npz`: `extract_physicelldfba_fields.py` on CF3-NC-CTL_s50
  (MN5, snapshots at 600, 1200 and 1800 min): the configuration of Fig. 6 run to 30 h. The raw output of CF3-NC50_s50 is
  not available, so the figure stops at 30 h.

In COMETS the static ring is restored only once per step, after the diffusion substeps (`applyStaticMedia` after
`diffuseMediaFick`; the `dirichlet` mask of `Utility.diffuseFick` is never set by the loader of 2.12.5), and at 10 um boxes
and 0.25 min the whole domain mixes within a step. The ring therefore exchanges at most 84/484 of the mixed volume per step
(test on a 7 x 7 grid, 2026-11-20), the fields in the domain are uniform and H2 is higher than in PhysiCelldFBA. S12.1
describes this.

## Rerun on mpws

Environment (`env/`):

- COMETS 2.12.5 Linux distribution in `$BENCHMARK_TOOLS/comets/comets_linux/comets_2.12.5` (jar sha256 `57b26959...`,
  OR-Tools 9.4.1874, unmodified `comets_scr`: `env/mpws_comets_launcher.txt`), OpenJDK 21.0.13.
- Python 3.12.3 venv `$BENCHMARK_TOOLS/venv` (cobra 0.32.1; full `pip freeze` in `env/mpws_comets_python.txt`). The
  cometspy that is imported is the one of the COMETS distribution (`lib/cometspy-master`, version 0.5.2), put on
  `PYTHONPATH`, ahead of the cometspy 0.6.3 of the venv.
- R 4.3.3 with `$BENCHMARK_TOOLS/R_libs`: BacArena 1.8.2, sybil 2.2.0, glpkAPI 1.3.4.1 on GLPK 5.0
  (`env/mpws_R_sessionInfo.txt`). The scripts set this library path themselves.

`BENCHMARK_TOOLS` defaults to `/home/ohayoun/benchmark_tools`. With Steps 3 and 4 of the top README done (they provide
`usecases/acetate/config/shipped/Ecoli_core.xml`, `usecases/crossfeeding/config/cells_CF3-REF_s50.csv`), from
`benchmark_frameworks/`:

    bash scripts/rerun_mpws.sh BM-CF-BACARENA runs/BM-CF-BACARENA

The script stages the inputs into the run folder, checks each one against `config/SHA256SUMS`, and runs the tool with the
settings of `runs.csv`. It then runs the extraction script. Copy the files listed in `data/<run folder>/` from the run
folder to compare or replace them. BacArena's own `Ec_core` (shipped with the package) is the model of BM-ACE-BACARENA.

### Inputs in `config/`

- `cb_modelorg.rds`, `mb_modelorg.rds`: the sybil modelorg objects of the two trimmed models
  (`usecases/crossfeeding/config/shipped/sbml_models/*_trimmed.xml`), read with sybilSBML on 2026-11-09. sybilSBML is no longer
  installed on mpws, so the objects are kept. Checked on 2026-11-19: same 501 / 551 reactions and bounds as the SBML files.
- `nocaps_bounds_{cb,mb}.csv`: the bounds that remove the FVA caps (485 / 522 reactions). These are exactly the reactions
  whose bounds differ between the trimmed models and `usecases/unit_test/models_full/*_trimmed_nocaps.xml` (checked
  2026-11-19). `run_bacarena_fig6_80x80.R` applies them to the modelorgs; the COMETS run reads the `_nocaps.xml` files
  directly.

## Branch

Branch `benchmark-comets-bacarena`, started from `bd51a5e` of `cf-acetate-kinetics-check`, because
`cf3nc50_series.csv` and the crossfeeding provenance it reads are on that branch. Merge it after that branch.
