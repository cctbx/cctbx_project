# cctbx-xfel build scripts

This directory collects build scripts for various environments. The general_linux
build is equivalent to the procedure in xfel/conda_envs/README.txt. Other recipes
are suitable for specific clusters and supercomputing centers such as the
Perlmutter supercomputer at NERSC.

## Procedure

We recommend cloning the cctbx_project repository and running the checked-out
build script. This provides for simpler version control in the cases where you
might need to build from a cctbx_project branch.

```bash
PREFIX=/path/to/install/location
git clone https://github.com/cctbx/cctbx_project $PREFIX/modules/cctbx_project
cd $PREFIX/modules/cctbx_project
git checkout <branch> # optional
PSANA_CHANNEL=lcls-i $PREFIX/modules/cctbx_project/xfel/conda_envs/site/general_linux/build.sh
# or
PSANA_CHANNEL=lcls-ii $PREFIX/modules/cctbx_project/xfel/conda_envs/site/nersc_perlmutter/build.sh
```

## psana version
We support two LCLS data formats, psana1 and psana2. This choice is controlled
by installing the `psana` conda package from channels `lcls-i` or `lcls-ii`
respectively. The build scripts will take your choice from the `PSANA_CHANNEL`
environment variable at build time.

## Activating the build

The build scripts write an activation script which should be sourced:

```bash
source $PREFIX/activate.sh
```

The activation script also defines a function `mk-cctbx` which enables partial
rebuilds. It takes the same verbs as `bootstrap.py`:

```bash
mk-cctbx build                  # Recompile only
mk-cctbx hot update build       # Re-fetch sources, build
mk-cctbx                        # Recompile only. Equivalent to mk-cctbx build
```

The conda environment is created only by `build.sh`; `mk-cctbx` never touches
it, so to rebuild it you need to start the build over.
