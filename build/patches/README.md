# Patches on top of the pinned code

`colony_quarter.patch` (from PhysiCelldFBA commit 23947c49, branch `feature/colony-quarter-disc`, 3 files, about 40 lines): a change to the **colony sample project only** (`bacterial_colony/custom_modules/custom.cpp`, `custom.h` and the user parameter `quarter_colony` in its settings file). With `quarter_colony` true the initial disc is created only in the quadrant x >= 0, y >= 0 (same rows, spacing and cells as the full disc). The symmetry planes are set in the run configuration (domain from 0; Dirichlet only on the outer faces), see `usecases/colony/overrides_def.csv`. The addon and the PhysiCell core are not touched.

Used by the full-length colony runs and the quarter-domain validation (binary c846ce99, `usecases/colony/binaries.tsv`). To rebuild that binary: export the pinned code, apply the patch (`git apply build/patches/colony_quarter.patch` in the export), then build the colony sample as `build/build_all.sh` does.
