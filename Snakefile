# Reproducibility workflow of pcdfba-review-runs (work in progress: the figure-collection step is the first one written).
# Planned order (a later rule feeds the earlier ones): configs (common/build_configs.py) -> runs (usecases/<uc>/scripts/run_task.sh, MN5 launchers) -> analysis -> figures -> paper_figures -> provenance workbooks -> Zenodo package.
# Existing outputs are inputs: no rule deletes or regenerates a run output.
# usage: snakemake -c1 paper_figures            (PAPER_FIGURES=<manuscript figures folder> for the figures that are not outputs of the runs)
import os
PAPER_FIGURES = os.environ.get("PAPER_FIGURES", "../pcdfba_paper/figures")

rule paper_figures:
    """The figures of the manuscript under their manuscript names (paper_figures/) and paper_figures/MANIFEST.csv"""
    input: "paper_figures_manifest.csv", "common/make_paper_figures.py"
    output: "paper_figures/MANIFEST.csv"
    shell: "python common/make_paper_figures.py --paper-figures {PAPER_FIGURES} --out paper_figures"
