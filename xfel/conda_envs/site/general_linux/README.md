# cctbx-general-linux-build

Builds cctbx.xfel with psana in a conda environment (installed to
`conda_base/` from `psana_environment.yml`), following the "General build"
procedure in [`xfel/conda_envs/README.md`](../../README.md).

## Running it

Clone cctbx_project into `modules/` under wherever you want the build to
live (the parent of `modules/` is the install prefix), then run `build.sh`
from a clean shell with no conda environment active:

    git clone https://github.com/cctbx/cctbx_project $PREFIX/modules/cctbx_project
    cd $PREFIX/modules/cctbx_project
    git checkout <branch>   # optional
    PSANA_CHANNEL=lcls-i $PREFIX/modules/cctbx_project/xfel/conda_envs/site/general_linux/build.sh

`PSANA_CHANNEL` is required: `lcls-i` builds against psana1 (pinned to
`psana=4.0.66`), `lcls-ii` builds against psana2 (unpinned, from the
`lcls-ii` channel instead). `build.sh` fills these into the env file it
generates from `psana_environment.yml`.

The conda solve takes several minutes; `bootstrap.py` can take 1-3 hours
depending on core count.

## Using the result

    source $PREFIX/build/conda_setpaths.sh

Re-running `build.sh` is safe but redoes everything. Paths are baked in, so
moving the checkout means rebuilding.

## Notes

* `build.sh` copies `bootstrap.py` from your `modules/cctbx_project`
  checkout, so it owns the repo list and configure flags.
* This directory's `psana_environment.yml` is a copy of
  `xfel/conda_envs/psana_environment.yml` with the channel and psana lines
  replaced by `__PSANA_CHANNEL__`/`__PSANA__` placeholders that `build.sh`
  fills in. See `xfel/conda_envs/README.md` for other variants (e.g.
  building with conda compilers).
* `NPROC` is auto-detected from `nproc`; edit the script if that isn't
  appropriate for your machine.
