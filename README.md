# pcdfba-release-runs

Definitive experiments for the PhysiCelldFBA revision. The plan and the settled decisions are in `docs/plan.md`.

**One code version.** `code/` is a submodule of PhysiCelldFBA pinned to the commit `7f2b1d0a` (tag `review-runs-v1`: `development` after PRs #20, #24, #25, #27).
`build/build_all.sh` builds one binary per use case from it, each in its own scratch tree; `build/build_record.tsv` records commit and binary hash.

**One folder per use case** (`usecases/<uc>/`), described by tables:

| file | content |
|---|---|
| `experiments.csv` | one row per experiment: question, reviewer comments it answers, what it feeds, observables, pass criterion, base config |
| `runs_def.csv` | the planned runs: `run_id`, `experiment_id`, `seed` (other `runs_*.csv` files are read as well) |
| `overrides.csv`, `overrides_def.csv` | `scope, xpath, value` on top of the base config; scope is `*` (all runs), an experiment id or a run id |
| `provenance.csv` | what actually happened per run (commit, binary, config hash, machine, times, status), written by `common/make_provenance.py` |
| `binaries.tsv` | binary hash to commit |

`scripts/design_def.py` writes `runs_def.csv` and `overrides_def.csv`; nothing is edited by hand in a generated config.
The six use cases are `unit_test`, `colony`, `acetate`, `crossfeeding`, `motility` and `tumour`.

## Reproducing a run from a fresh clone

The folders `usecases/<uc>/bin/`, `config/` and `runs/` are **empty after a clone on purpose**: binaries, generated configs and outputs are not tracked, they are produced by the steps below. What is tracked is the code pin, the tables that describe the experiments, the scripts, and the small result tables and figures.

### Step 0. System packages (Ubuntu / Debian)

The addon links libSBML and Clp statically (step 2 downloads them) and needs a few system libraries and a compiler:

    sudo apt install build-essential git python3 python3-pip python3-venv libxml2-dev libbz2-dev liblapack-dev zlib1g-dev

`libxml2-dev`, `libbz2-dev` and `liblapack-dev` are the ones that are not there on a new machine; the link line of the samples is `-lsbml-static -lxml2 -lbz2 -lz -llapack` and fails without them. The generators and analysis scripts need Python with `numpy scipy pandas lxml matplotlib` (and `cobra` for the model checks):

    python3 -m venv ~/.venvs/pcdfba && . ~/.venvs/pcdfba/bin/activate
    pip install numpy scipy pandas lxml matplotlib cobra

### Step 1. Clone and fetch the pinned code

    git clone git@github.com:PhysiCelldFBA/pcdfba-release-runs.git
    cd pcdfba-release-runs
    git submodule update --init

`code/` is not part of this repository: it is a pointer to one commit of PhysiCelldFBA. The first command gives an empty `code/` folder, the second downloads that commit into it (`git submodule status` must show `7f2b1d0a...` without a `+` or `-` in front). Stay on `main`: a branch of this repository that predates the use cases has most of the folders empty or missing.

### Step 2. Vendored libraries (`ext`) and the variable `PCDFBA_EXT`

The addon needs libSBML 5.20.4 and Clp 1.17.10. They are not in git; the pinned code has a script that downloads them (network needed) into `code/addons/dFBA/ext`, which is git-ignored there:

    cd code
    python3 beta/setup_fba.py          # downloads libSBML and Clp into addons/dFBA/ext (skipped if already there)
    export PCDFBA_EXT=$PWD/addons/dFBA/ext
    cd ..

`PCDFBA_EXT` is read by `build/build_all.sh` and has no default: it must be the directory that contains `coin-or/` and `libsbml/` (check that `$PCDFBA_EXT/coin-or/include/coin/CoinPackedMatrix.hpp` exists: the script only prints "No package information found ... Skipping" when there is no download for your OS or architecture, and does not fail). If you already have these libraries from another PhysiCelldFBA checkout, point the variable at that directory instead and skip the download.

### Step 3. Build the binaries and the shipped base configs

    build/build_all.sh                 # all use cases listed in build/usecases.tsv
    build/build_all.sh tumour colony   # or only some of them (names as in the first column of build/usecases.tsv)
    JOBS=16 build/build_all.sh tumour  # make parallelism, default 8

For each use case the script (i) exports the pinned code with `git archive` into its own scratch tree `build/trees/<uc>/`, so `code/` is never modified (it refuses to run if `code/` has uncommitted changes); (ii) links `$PCDFBA_EXT` into that tree; (iii) copies the sample project of the use case over the root of the export, which is what the staging targets of PhysiCelldFBA do; (iv) runs `make` there (log in `build/trees/<uc>.build.log`); (v) copies the binary to `usecases/<uc>/bin/`; (vi) copies the sample's own `config/` to `usecases/<uc>/config/shipped/`, the base that the generator below starts from; (vii) appends commit and binary hash to `build/build_record.tsv`. The sample directory and binary name of each use case are in `build/usecases.tsv`.

The colony runs that use a quarter of the disc need the binary rebuilt with `build/patches/colony_quarter.patch` (see `build/patches/README.md`); the plain `build_all.sh colony` binary is the unpatched one.

### Step 4. Generate the configs from the tables

    python3 common/build_configs.py usecases/tumour                    # every run of the use case
    python3 common/build_configs.py usecases/tumour TUM5-REF_s50       # only the given run ids

For every run in the `runs*.csv` files this reads the base config (`config/shipped/...`, named in `experiments.csv`), applies the rows of `overrides*.csv` that match the run (scope `*`, then its experiment, then its run id) and writes `usecases/<uc>/config/<run_id>.xml` (plus `cells_<run_id>.csv` for use cases with a cell generator). Every xpath must match exactly one element, so a typo or a changed base fails loudly and nothing is written silently wrong. The seed is set from the `seed` column and the output folder to `./runs/<run_id>`.

### Step 5. Run

From the use-case folder, one run is the binary with the generated config:

    cd usecases/tumour
    OMP_NUM_THREADS=56 bin/cancer-core-metabolism config/TUM5-REF_s50.xml

Output goes to `runs/<run_id>/` (not tracked, never delete it). Each use case has `scripts/run_task.sh RUN_ID`, which does the same, refuses a non-empty `runs/<run_id>` and writes `run.log` with start and end lines; the thread count is taken from the config. A full-length run takes minutes to hours on 56 cores, so the larger sets were run on MareNostrum 5 (below).

### Step 6. Check and analyse

    python3 common/check_run.py usecases/<uc> RUN_ID    # saved settings of the run == the rebuilt config (exit 1 on any difference)
    python3 common/make_provenance.py usecases/<uc>     # writes provenance.csv from the runs' own output
    python3 usecases/<uc>/scripts/analyse_*.py ...      # per use case; see usecases/<uc>/RESULTS.md and the docstring of each script

`RESULTS.md` of a use case lists the experiments, the numbers and the open items; each script starts with a docstring that gives its usage.

### Step 7. Redraw the figures of the manuscript

`paper_figures_manifest.csv` has one row per figure of the manuscript (main and supplementary): the manuscript file name, the use case, the runs it is drawn from (`source_runs`), the command that draws it (`command`, run from `usecases/<use_case>/`) and the file it writes (`source_file`). For example, Figure 2:

    cd usecases/unit_test && python3 scripts/make_figure.py runs/UT-BASE_s50 results/fig_unit_test_biomass

The runs it names must be in `runs/` (Step 5, or the simulation outputs archived on Zenodo). Then, from the repository root,

    python3 common/make_paper_figures.py --paper-figures <folder of the manuscript figures>

copies every figure under its manuscript name into `paper_figures/` and writes `paper_figures/MANIFEST.csv` with the SHA-256 of each file and whether it differs from the manuscript's copy. Rows with status `external` (Figure 1, the PhysiCell-Studio screenshot, the workflow drawing and the graphical abstract) are not produced by a simulation and have no command.

## Running on MareNostrum 5

Each use case is run from a directory that holds `tree/` (the exported pinned code with the sample copied over and the vendored `addons/dFBA/ext` linked in) and `uc/` (a copy of the use-case directory with the generated `config/`, `scripts/` and an empty `runs/`).
Submit from that directory: `sbatch uc/scripts/<name>.sbatch`. The scripts take that directory from `PCDFBA_MN5_DIR`, or from where the job was submitted. The account `bsc08` is the project group code.

## Licence

The scripts and code in this repository are distributed under the BSD 3-Clause licence (`LICENSE`). Data, figures and tables
(the run outputs, provenance tables and figure files) are distributed under the Creative Commons Attribution 4.0 licence
(CC-BY-4.0). The PhysiCelldFBA source in `code/` keeps its own licence (BSD 3-Clause; PhysiCell and BioFVM notices in its
`licenses/` folder), and the metabolic model files keep the terms of their sources (see `licenses/THIRD_PARTY.md` of
PhysiCelldFBA).
