# cctbx-perlmutter-build

Builds cctbx/DIALS on NERSC Perlmutter with psana, CUDA and Kokkos, in a conda
environment (`psana_env`, from `perlmutter_environment.yml`) with mpi4py built
against Cray MPICH.

## Running it

Clone on `$CFS` or `$SCRATCH` (the checkout directory is the install prefix)
and run from a login node with no conda environment active:

    git clone <this repo> $SCRATCH/cctbx-perlmutter-build
    cd $SCRATCH/cctbx-perlmutter-build
    script -q -e -c ./build.sh build.log

The conda solve takes 5-10 minutes and `bootstrap.py` 1-3 hours. A GPU is not
needed to build. Use `script` rather than `| tee` for the log so progress bars
keep working; `col -b < build.log` strips the control characters.

## Using the result

    source $SCRATCH/cctbx-perlmutter-build/activate.sh

`activate.sh` also defines `mk-cctbx` for rebuilding without redoing the conda
environment:

    mk-cctbx build              # recompile
    mk-cctbx build hot update   # re-fetch sources, then recompile

Re-running `build.sh` is safe but redoes everything. Paths are baked in, so
moving the checkout means rebuilding.

## Notes

* `bootstrap.py` is pinned to a cctbx commit (`BOOTSTRAP_COMMIT` in `build.sh`);
  it owns the repo list, the Kokkos pin and the configure flags. Record the SHA
  of whichever build succeeds.
* Last checked on Perlmutter 2026-08-17: `PrgEnv-gnu/8.6.0`, `cudatoolkit/12.9`
  and `cpe-cuda/25.09` exist; `cc -shared -target-accel=nvidia80 -lmpi -lgdrapi`
  works. `SIT_ROOT=/reg/g/psdm` in `activate.sh` is an SLAC path and unverified.
* "could not clone ksdssp" is expected (private repo) and not fatal.

## Open item: CUDA 13

As of 2026-08-26 the unversioned `cudatoolkit` module gives CUDA 13, which the
Kokkos 4.2.00 pinned in cctbx master and cctbx's CUDA arch flags do not
support. `build.sh` therefore checks out the cctbx branch `cuda13-kokkos-4.7`
(Kokkos 4.7.03 plus CUDA >= 13 arch flags). Revert that commit once the branch
is merged upstream.
