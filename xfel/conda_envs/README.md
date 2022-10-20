# cctbx.xfel Build Instructions

The conda build environment psana_environment.yml is suitable
for general use and contains the usual CCTBX dependencies plus psana and its
dependencies.

The build steps below were tested on Oct 19, 2022. They should be done in a clean environment: start
a new shell before proceeding.

## Prerequisite: Install Miniconda3 and add mamba

If needed, visit: https://docs.conda.io/en/latest/miniconda.html and install the correct Miniconda
for your platform. Activate your `base` environment and do: `conda install mamba -c conda-forge`.
Mamba is a much faster C++ implementation of conda.

## General build

These steps were tested on a CentOS 7.9 machine with 64 cores. In the
bootstrap.py step you should adjust nproc to suit your environment.

```
$ mkdir cctbx; cd cctbx
$ wget https://raw.githubusercontent.com/cctbx/cctbx_project/master/libtbx/auto_build/bootstrap.py
$ wget https://raw.githubusercontent.com/cctbx/cctbx_project/master/xfel/conda_envs/psana_environment.yml
$ mamba env create -f psana_environment.yml -p $PWD/conda_base
$ conda activate `pwd`/conda_base
$ python bootstrap.py --builder=xfel --use-conda=$PWD/conda_base --nproc=48 \
    --python=39 --no-boost-src hot update build
$ echo $PWD/build/conda_setpaths.sh
```
To activate the cctbx environment, `source` the script that was printed in the final step.

## LCLS build

Since the `psana` compute nodes do not have internet access, we use `psexport` for everything except `build`.
```
dwpaley@pslogin02:~
$ ssh psexport
[...]
$ cd /reg/d/psdm/<experiment>/scratch/dwpaley
$ mkdir cctbx; cd cctbx
$ wget https://raw.githubusercontent.com/cctbx/cctbx_project/master/libtbx/auto_build/bootstrap.py
$ wget https://raw.githubusercontent.com/cctbx/cctbx_project/master/xfel/conda_envs/psana_environment.yml
$ mamba env create -f psana_environment.yml -p $PWD/conda_base
$ conda activate $PWD/conda_base
$ python bootstrap.py --builder=xfel --use-conda=$PWD/conda_base --nproc=48 \
    --python=39 --no-boost-src hot update
$ exit
$ ssh psana
[...]
$ cd /reg/d/psdm/<experiment>/scratch/dwpaley/cctbx
$ conda activate $PWD/conda_base
$ python bootstrap.py --builder=xfel --use-conda=$PWD/conda_base --nproc=12 \
    --python=39 --no-boost-src build
$ echo $PWD/build/conda_setpaths.sh
```

## Build with conda compilers

On machines that are old or weird, the system compilers may not work correctly. We can have `bootstrap.py` build
with standardized compilers from conda instead. Replace the step `python bootstrap.py ... build` with the following:
```
$ python bootstrap.py --builder=xfel --use-conda=$PWD/conda_base \
  --config-flags="--compiler=conda" --config-flags="--use_environment_flags" \
  --nproc=10 --python=39 --no-boost-src build
```

# cctbx.xfel tests

The cctbx.xfel regression tests include tests from several repositories.  The below instructions reproduce what we do nightly. If psana is configured, it will be tested as well.

```
$ cd modules
$ conda install -c conda-forge git-lfs
$ git clone https://gitlab.com/cctbx/xfel_regression.git
$ git clone https://github.com/nksauter/LS49.git
$ git clone https://gitlab.com/cctbx/ls49_big_data.git
$ cd xfel_regression
$ git lfs install --local
$ git lfs pull
$ cd ../uc_metrics
$ git lfs install --local
$ git lfs pull
$ cd ../ls49_big_data
$ git lfs install --local
$ git lfs pull
$ cd ../../
$ mkdir test; cd test
$ libtbx.configure xfel_regression LS49 ls49_big_data
$ export OMP_NUM_THREADS=4
$ libtbx.run_tests_parallel module=uc_metrics module=simtbx module=xfel_regression module=LS49 nproc=64
```

Note, bootstrap.py has several 'builders' available that set up which packages are cloned and configured.  The xfel builder will clone uc_metrics for you, but for reference, here's how to get it standalone if needed:

```
$ git clone https://gitlab.com/cctbx/uc_metrics.git
$ libtbx.configure uc_metrics
$ cd `libtbx.show_build_path`; make
```

