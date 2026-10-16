# cctbx-perlmutter-build

Builds cctbx/DIALS on NERSC Perlmutter with psana, CUDA and Kokkos, in a conda
environment (installed to `conda_base/` from `perlmutter_environment.yml`)
with mpi4py built against Cray MPICH.

## Running it

Clone cctbx_project into `modules/` under wherever you want the build to live
on `$CFS` or `$SCRATCH` (the parent of `modules/` is the install prefix), then
run `build.sh` from a login node with no conda environment active:

    git clone https://github.com/cctbx/cctbx_project $SCRATCH/cctbx/modules/cctbx_project
    cd $SCRATCH/cctbx/modules/cctbx_project
    git checkout <branch>   # optional, e.g. cuda13-kokkos-4.7 -- see Open item below
    cd $SCRATCH/cctbx
    PSANA_CHANNEL=lcls-i script -q -e -c modules/cctbx_project/xfel/conda_envs/site/nersc_perlmutter/build.sh build.log

`PSANA_CHANNEL` is required: `lcls-i` builds against psana1 (pinned to
`psana=4.0.66`), `lcls-ii` builds against psana2 (unpinned, from the
`lcls-ii` channel instead). `build.sh` fills these into the env file it
generates from `perlmutter_environment.yml`.

The conda solve takes 5-10 minutes and `bootstrap.py` 1-3 hours. A GPU is not
needed to build. Use `script` rather than `| tee` for the log so progress bars
keep working; `col -b < build.log` strips the control characters.

## Using the result

    source $SCRATCH/cctbx/activate.sh

`activate.sh` also defines `mk-cctbx` for rebuilding without redoing the conda
environment:

    mk-cctbx build              # recompile
    mk-cctbx build hot update   # re-fetch sources, then recompile

Paths are baked in, so moving the checkout means rebuilding. If
`conda_base/` already exists, `build.sh` requires `CCTBX_RECREATE_CONDA` to be
set: `1` recreates the conda environment (including the mpi4py-against-Cray-MPICH
rebuild), `0` reuses it as-is and skips straight to the `bootstrap.py` build.

## Notes

* `build.sh` copies `bootstrap.py` from your `modules/cctbx_project` checkout,
  so it owns the repo list, the Kokkos pin and the configure flags. `bootstrap.py`'s
  `hot update` runs `git pull --rebase` on an existing `modules/cctbx_project`
  checkout, so whatever branch you checked out survives the build.
* Last checked on Perlmutter 2026-08-17: `PrgEnv-gnu/8.6.0`, `cudatoolkit/12.9`
  and `cpe-cuda/25.09` exist; `cc -shared -target-accel=nvidia80 -lmpi -lgdrapi`
  works. `SIT_ROOT=/reg/g/psdm` in `activate.sh` is an SLAC path and unverified.
* "could not clone ksdssp" is expected (private repo) and not fatal.

## Open item: CUDA 13

As of 2026-08-26 the unversioned `cudatoolkit` module gives CUDA 13, which the
Kokkos 4.2.00 pinned in cctbx master and cctbx's CUDA arch flags do not
support. Check out the cctbx branch `cuda13-kokkos-4.7` (Kokkos 4.7.03 plus
CUDA >= 13 arch flags) before running `build.sh`, until that branch is merged
upstream.
