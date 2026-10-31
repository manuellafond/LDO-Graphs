#!/bin/bash
# Install the software used for the results: OrthoFinder 3 (GitHub main, pinned commit) with the opt-in
# normalisation patch, in a Python virtual environment that also runs all the pipeline scripts
# (numpy, which the scripts need, comes with OrthoFinder).
# BLAST+ (blastn, makeblastdb) must be installed separately: OrthoFinder searches DNA with blastn.
#
#   bash setup.sh [ENV_DIR]          (default: ./env)
#   PYTHON=python3.11 bash setup.sh   to choose the Python used for the environment (default python3); OrthoFinder at
#                                    this commit needs Python >= 3.11
# then:   source ENV_DIR/env.sh
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
ENV=$(realpath -m "${1:-$HERE/env}")
COMMIT=c85f0034679a9812959655f1cd60c74a30597cf1      # OrthoFinder main, 2026-05-13 (reports v3.1.4.post1.dev1)
PYTHON=${PYTHON:-python3}

# OrthoFinder at this commit requires Python >= 3.11: check before cloning anything.
if ! "$PYTHON" -c 'import sys; sys.exit(sys.version_info < (3, 11))' 2>/dev/null; then
  if command -v "$PYTHON" > /dev/null; then FOUND=$("$PYTHON" -V 2>&1); else FOUND="not found"; fi
  echo "ERROR: $PYTHON: $FOUND; OrthoFinder needs Python >= 3.11." >&2
  echo "       Install one (e.g. apt-get install python3.11 python3.11-venv, conda create -n py311 python=3.11," >&2
  echo "       or uv python install 3.11) and run:  PYTHON=python3.11 bash setup.sh" >&2
  exit 1
fi

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
# a venv left by an earlier run with another Python is rebuilt
if [ -x "$ENV/venv/bin/python" ] && ! "$ENV/venv/bin/python" -c 'import sys; sys.exit(sys.version_info < (3, 11))'; then
  echo "removing $ENV/venv (built with $("$ENV/venv/bin/python" -V 2>&1))"; rm -rf "$ENV/venv"
fi
"$PYTHON" -m venv "$ENV/venv"
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
