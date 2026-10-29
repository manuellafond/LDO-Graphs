#!/bin/bash
# Install the software used for the results: OrthoFinder 3 (GitHub main, pinned commit) with the opt-in
# normalisation patch, in a Python virtual environment that also runs all the pipeline scripts
# (numpy, which the scripts need, comes with OrthoFinder).
# BLAST+ (blastn, makeblastdb) must be installed separately: OrthoFinder searches DNA with blastn.
#
#   bash setup.sh [ENV_DIR]          (default: ./env)
# then:   source ENV_DIR/env.sh
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
ENV=$(realpath -m "${1:-$HERE/env}")
COMMIT=c85f0034679a9812959655f1cd60c74a30597cf1      # OrthoFinder main, 2026-05-13 (reports v3.1.4.post1.dev1)

mkdir -p "$ENV"
if [ ! -d "$ENV/OrthoFinder" ]; then
  git clone https://github.com/OrthoFinder/OrthoFinder.git "$ENV/OrthoFinder"
fi
cd "$ENV/OrthoFinder"
git checkout -q "$COMMIT"
if git apply --check "$HERE/orthofinder_waterfall_norm.patch" 2>/dev/null; then
  git apply "$HERE/orthofinder_waterfall_norm.patch"
  echo "normalisation patch applied (inactive unless OF_NORM_* variables are set)"
else
  git apply --reverse --check "$HERE/orthofinder_waterfall_norm.patch" && echo "patch already applied"
fi
python3 -m venv "$ENV/venv"
"$ENV/venv/bin/pip" install -q "$ENV/OrthoFinder"

cd /    # not from inside the clone, whose orthofinder.py would shadow the installed package
SITE=$("$ENV/venv/bin/python" -I -c "import orthofinder, os; print(os.path.dirname(orthofinder.__file__))")
cat > "$ENV/env.sh" <<EOF
# source this file before swisstree/reproduce.sh or simphy/simphy_run.sh
export PY="$ENV/venv/bin/python"
export OF="$ENV/venv/bin/orthofinder"
export MAFFT="$SITE/bin/mafft/bin/mafft"
export MAFFT_BINARIES="$SITE/bin/mafft/libexec"
EOF
"$ENV/venv/bin/orthofinder" --version | sed 's/\x1b\[[0-9;]*m//g' | tail -1
if command -v blastn > /dev/null && command -v makeblastdb > /dev/null; then
  echo "BLAST+: $(blastn -version | head -1)"
else
  echo "WARNING: blastn / makeblastdb not found. Install BLAST+ (e.g. apt-get install ncbi-blast+, or"
  echo "         conda install -c bioconda blast) before running OrthoFinder on DNA."
fi
echo "done; now: source $ENV/env.sh"
