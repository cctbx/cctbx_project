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
MINIFORGE_PREFIX=$BUILD_PREFIX/mc3
ENV_PREFIX=$BUILD_PREFIX/conda_base
NPROC=$(nproc)

SKIP_CONDA=
if [[ -e $ENV_PREFIX ]]; then
    case "${CCTBX_RECREATE_CONDA:-}" in
        1) ;;
        0) SKIP_CONDA=1 ;;
        *) echo "$ENV_PREFIX exists; set CCTBX_RECREATE_CONDA=1 (recreate) or 0 (reuse)" >&2; exit 1 ;;
    esac
fi

MINIFORGE_VERSION=26.3.2-3
MINIFORGE=Miniforge3-$MINIFORGE_VERSION-Linux-x86_64.sh
MINIFORGE_URL=https://github.com/conda-forge/miniforge/releases/download/$MINIFORGE_VERSION/$MINIFORGE

step () { printf '\n===== %s -- %s\n' "$*" "$(date '+%F %T')"; }

cd "$BUILD_PREFIX"
mkdir -p opt

step "copy bootstrap.py and download Miniforge3"
cp modules/cctbx_project/libtbx/auto_build/bootstrap.py bootstrap.py
wget -O "opt/$MINIFORGE" "$MINIFORGE_URL"

step "install Miniforge3 into $MINIFORGE_PREFIX"
bash "opt/$MINIFORGE" -b -u -p "$MINIFORGE_PREFIX"   # -b: no prompts, no conda init.  -u: re-runnable.

source "$MINIFORGE_PREFIX/etc/profile.d/conda.sh"

if [[ -z $SKIP_CONDA ]]; then
    step "conda env create (several minutes)"
    case "$PSANA_CHANNEL" in
        lcls-i)  PSANA_SPEC='psana=4.0.66' ;;
        lcls-ii) PSANA_SPEC='psana' ;;
    esac
    sed -e "s/__PSANA_CHANNEL__/$PSANA_CHANNEL/" -e "s/__PSANA__/$PSANA_SPEC/" \
        "$SCRIPT_DIR/psana_environment.yml" > "$BUILD_PREFIX/opt/environment.yml"
    grep -q '__PSANA' "$BUILD_PREFIX/opt/environment.yml" \
        && { echo "placeholder substitution into opt/environment.yml did not fully apply" >&2; exit 1; }
    conda env create -f "$BUILD_PREFIX/opt/environment.yml" -p "$ENV_PREFIX" --yes -v

    step "install gemmi"
    conda install -y -p "$ENV_PREFIX" -c conda-forge --no-deps gemmi
fi

conda activate "$ENV_PREFIX"

step "bootstrap.py hot update build"
python bootstrap.py --builder=xfel --use-conda="$ENV_PREFIX" --nproc="$NPROC" \
    --no-boost-src hot update build

printf '\nDone.  Source %s/build/conda_setpaths.sh to use this installation.\n' "$BUILD_PREFIX"
