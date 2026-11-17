# Supplementary Data workbook

`build_workbook.py` builds `Supplementary_Data.xlsx` (the numerical tables of the convergence, sensitivity and cost
analyses) from the LaTeX source of the tables in `tables_source.tex` (taken from the supplement of the manuscript).
`table_map.py` maps each supplementary table (number before the workbook, LaTeX label) to its workbook sheet and its
number in the PDF.

    python build_workbook.py            # reads tables_source.tex, writes Supplementary_Data.xlsx

Needs openpyxl. The values are those of the manuscript tables; the runs behind them are listed in
`../../paper_figures_manifest.csv` and the `RESULTS.md` file of each use case.
