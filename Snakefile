# Reproducibility workflow of pcdfba-release-runs (work in progress: the figure-collection step is the first one written).
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

rule provenance_workbooks:
    """One Excel workbook per use case (provenance/<uc>_provenance.xlsx) with a CSV copy of every sheet and provenance/all_runs.csv; the provenance.csv files are written first by common/make_provenance.py"""
    input: expand("usecases/{uc}/provenance.csv", uc=["unit_test", "colony", "acetate", "crossfeeding", "motility", "tumour"]), "common/make_provenance_workbooks.py"
    output: "provenance/all_runs.csv"
    shell: "python common/make_provenance_workbooks.py"
