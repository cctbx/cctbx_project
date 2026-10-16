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

case "${PSANA_CHANNEL:-}" in
    lcls-i|lcls-ii) ;;
    *) echo "set PSANA_CHANNEL=lcls-i (psana1) or lcls-ii (psana2)" >&2; exit 1 ;;
esac

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
bash "opt/$MINIFORGE" -b -u -p "$BUILD_PREFIX/mc3"   # -b: no prompts, no conda init.  -u: re-runnable.

source "$BUILD_PREFIX/mc3/etc/profile.d/conda.sh"

step "conda env create (several minutes)"
case "$PSANA_CHANNEL" in
    lcls-i)  PSANA_SPEC='psana=4.0.66' ;;
    lcls-ii) PSANA_SPEC='psana' ;;
esac
sed -e "s/__PSANA_CHANNEL__/$PSANA_CHANNEL/" -e "s/__PSANA__/$PSANA_SPEC/" \
    "$SCRIPT_DIR/psana_environment.yml" > "$BUILD_PREFIX/opt/environment.yml"
grep -q '__PSANA' "$BUILD_PREFIX/opt/environment.yml" \
    && { echo "placeholder substitution into opt/environment.yml did not fully apply" >&2; exit 1; }
conda env create -f "$BUILD_PREFIX/opt/environment.yml" -p "$CONDA_BASE" --yes -v

step "install gemmi"
conda install -y -p "$CONDA_BASE" -c conda-forge --no-deps gemmi

conda activate "$CONDA_BASE"

step "bootstrap.py hot update build"
python bootstrap.py --builder=xfel --use-conda="$CONDA_BASE" --nproc="$NPROC" \
    --no-boost-src hot update build

printf '\nDone.  Source %s/build/conda_setpaths.sh to use this installation.\n' "$BUILD_PREFIX"
