# pcdfba-review-runs

Definitive experiments for the PhysiCelldFBA revision. The plan and the settled decisions are in `docs/plan.md`.

**One code version.** `code/` is a submodule of PhysiCelldFBA pinned to the commit `7f2b1d0a` (tag `review-runs-v1`: `development` after PRs #20, #24, #25, #27).
`build/build_all.sh` builds one binary per use case from it, each in its own scratch tree; `build/build_record.tsv` records commit and binary hash.

**One folder per use case** (`usecases/<uc>/`), described by tables:

| file | content |
|---|---|
| `experiments.csv` | one row per experiment: question, reviewer comments it answers, what it feeds, observables, pass criterion, base config |
| `runs_def.csv` | the planned runs: `run_id`, `experiment_id`, `seed` |
| `overrides.csv`, `overrides_def.csv` | `scope, xpath, value` on top of the base config; scope is `*` (all runs), an experiment id or a run id |
| `provenance.csv` | what actually happened per run (commit, binary, config hash, machine, times, status), written by `common/make_provenance.py` |
| `binaries.tsv` | binary hash to commit |

`scripts/design_def.py` writes `runs_def.csv` and `overrides_def.csv`; nothing is edited by hand in a generated config.

    git submodule update --init
    PCDFBA_EXT=/path/to/ext build/build_all.sh [usecase ...]       # binaries and the shipped base configs
    python3 common/build_configs.py usecases/<uc> [RUN_ID ...]      # configs (and cells files) from base + overrides
    python3 common/check_run.py usecases/<uc> RUN_ID                # saved settings of a run == the rebuilt config
    python3 common/make_provenance.py usecases/<uc>

The seed is `<options><random_seed>` only. Run outputs (`runs/`) and binaries (`bin/`) are not tracked; figures, tables and the provenance records are.

**Running on MareNostrum 5.** Each use case is run from a directory that holds `tree/` (the exported pinned code with the sample copied over and the vendored `addons/dFBA/ext` linked in) and `uc/` (a copy of the use-case directory with the generated `config/`, `scripts/` and an empty `runs/`).
Submit from that directory: `sbatch uc/scripts/<name>.sbatch`. The scripts take that directory from `PCDFBA_MN5_DIR`, or from where the job was submitted. The account `bsc08` is the project group code.
The Python tools need numpy, scipy, pandas, lxml and matplotlib (and cobra for the model checks); set `PYTHON` where a script uses it.
