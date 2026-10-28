#!/bin/bash
# Build every use-case binary from the one pinned copy of the code (code/ submodule), each in its own scratch tree.
# usage: build/build_all.sh [usecase ...]      (default: all in build/usecases.tsv)
# env:   PCDFBA_EXT   directory holding the vendored libSBML/Clp (the content of addons/dFBA/ext)  [required]
#        JOBS         make parallelism (default 8)
# Nothing in code/ is modified: the tree is exported with git archive, the sample is copied over the root of the export
# (this is what the staging targets do, without their Makefile-backup trap), and make runs there.
set -euo pipefail
HERE=$(cd "$(dirname "$0")/.." && pwd)
: "${PCDFBA_EXT:?set PCDFBA_EXT to the vendored ext directory}"
JOBS=${JOBS:-8}
COMMIT=$(git -C "$HERE/code" rev-parse HEAD)
git -C "$HERE/code" diff --quiet && git -C "$HERE/code" diff --cached --quiet || { echo "code/ has uncommitted changes: refusing"; exit 1; }
WANT=("$@"); [ ${#WANT[@]} -eq 0 ] && mapfile -t WANT < <(grep -v '^#' "$HERE/build/usecases.tsv" | cut -f1)
for uc in "${WANT[@]}"; do
  read -r _ sample binary _ < <(grep -P "^$uc\t" "$HERE/build/usecases.tsv")
  tree="$HERE/build/trees/$uc"
  rm -rf "$tree"; mkdir -p "$tree"
  git -C "$HERE/code" archive "$COMMIT" | tar -x -C "$tree"
  rm -rf "$tree/addons/dFBA/ext"; ln -s "$PCDFBA_EXT" "$tree/addons/dFBA/ext"
  S="$tree/sample_projects_intracellular/fba/$sample"
  mkdir -p "$tree/custom_modules" "$tree/config" "$tree/scripts"
  cp -r "$S/custom_modules/." "$tree/custom_modules/"; cp "$S/main.cpp" "$tree/main.cpp"; cp "$S/Makefile" "$tree/Makefile"
  [ -d "$S/config" ] && cp -r "$S/config/." "$tree/config/"
  ( cd "$tree" && make -j"$JOBS" PHYSICELL_CPP=g++ PROGRAM_NAME="$binary" > "$HERE/build/trees/$uc.build.log" 2>&1 ) || { echo "BUILD FAILED: $uc (see build/trees/$uc.build.log)"; exit 1; }
  cp "$tree/$binary" "$HERE/usecases/$uc/bin/$binary"
  # the shipped sample config is kept next to the binary as the base for the generators
  mkdir -p "$HERE/usecases/$uc/config/shipped"; cp -r "$S/config/." "$HERE/usecases/$uc/config/shipped/"
  printf '%s\t%s\t%s\t%s\n' "$uc" "$binary" "$COMMIT" "$(sha256sum "$HERE/usecases/$uc/bin/$binary" | cut -c1-16)" >> "$HERE/build/build_record.tsv"
  echo "built $uc -> usecases/$uc/bin/$binary @ ${COMMIT:0:8}"
done
