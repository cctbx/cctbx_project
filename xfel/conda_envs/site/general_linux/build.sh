#!/bin/bash
#
# Build cctbx.xfel with psana in a conda environment, following the "General
# build" procedure in xfel/conda_envs/README.md.
# Run it, do not source it.
#
# modules/cctbx_project is your clone of this repo; build/, mc3/, conda_base/,
# opt/ and bootstrap.py are created next to modules/.  Clone where you want
# the build to live:
#
#   git clone https://github.com/cctbx/cctbx_project $PREFIX/modules/cctbx_project
#   $PREFIX/modules/cctbx_project/xfel/conda_envs/site/general_linux/build.sh

set -eo pipefail

[[ -z ${CONDA_PREFIX:-} ]] || { echo "conda env active; run 'conda deactivate' first" >&2; exit 1; }

SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
BUILD_PREFIX=$(cd "$SCRIPT_DIR/../../../../../.." && pwd)
CONDA_BASE=$BUILD_PREFIX/conda_base
NPROC=$(nproc)

MINIFORGE_VERSION=26.3.2-3
MINIFORGE=Miniforge3-$MINIFORGE_VERSION-Linux-x86_64.sh
MINIFORGE_URL=https://github.com/conda-forge/miniforge/releases/download/$MINIFORGE_VERSION/$MINIFORGE

step () { printf '\n===== %s -- %s\n' "$*" "$(date '+%F %T')"; }

cd "$BUILD_PREFIX"
mkdir -p opt

step "copy bootstrap.py and download Miniforge3"
cp modules/cctbx_project/libtbx/auto_build/bootstrap.py bootstrap.py
wget -O "opt/$MINIFORGE" "$MINIFORGE_URL"

step "install Miniforge3 into $BUILD_PREFIX/mc3"
bash "opt/$MINIFORGE" -b -p "$BUILD_PREFIX/mc3"

source "$BUILD_PREFIX/mc3/etc/profile.d/conda.sh"

step "conda env create (several minutes)"
conda env create -f modules/cctbx_project/xfel/conda_envs/psana_environment.yml -p "$CONDA_BASE"

step "install gemmi"
conda install -y -p "$CONDA_BASE" -c conda-forge --no-deps gemmi

conda activate "$CONDA_BASE"

step "bootstrap.py hot update build"
python bootstrap.py --builder=xfel --use-conda="$CONDA_BASE" --nproc="$NPROC" \
    --no-boost-src hot update build

printf '\nDone.  Source %s/build/conda_setpaths.sh to use this installation.\n' "$BUILD_PREFIX"
