# cctbx.xfel Build Instructions

The conda build environment psana_environment.yml is suitable
for general use and contains the usual CCTBX dependencies plus psana and its
dependencies.

The build steps below were tested on Aug 17, 2022. They should be done in a
clean environment (no cctbx installation has been activated, etc.)

Creating the conda environment (step `base`) may take up to ~20 min with no
obvious progress occurring.

## General build

These steps were tested on a CentOS 7.9 machine with 12 cores. In the
bootstrap.py step you should adjust nproc to suit your environment.

```
$ mkdir cctbx; cd cctbx
$ wget https://raw.githubusercontent.com/cctbx/cctbx_project/master/libtbx/auto_build/bootstrap.py
$ wget https://raw.githubusercontent.com/cctbx/cctbx_project/master/xfel/conda_envs/psana_environment.yml
$ python bootstrap.py --builder=xfel --use-conda=psana_environment.yml --nproc=12 --python=39 --no-boost-src hot update base
$ conda activate `pwd`/conda_base # if no conda is availble, first source mc3/etc/profile.d/conda.sh
$ python bootstrap.py --builder=xfel --use-conda=psana_environment.yml --nproc=12 --python=39 build
$ source build/conda_setpaths.sh
$ libtbx.python -c "import psana" # Should exit with no output
```

## LCLS build

Follow the general build proceedure but first run the steps that require internet access using
an ssh connection to pslogin.slac.stanford.edu (up through the hot update base step). Use the psana
or psana-ffb nodes for the build step.

## Alternative build using mamba and the conda compilers

The `bootstrap.py base` step above takes >1 hr, mainly to create the conda environment. This can be improved to ~15 min using
Mamba, a C++ implementation of Conda. You need a base Conda environment with Mamba installed (`conda install mamba -c conda-forge`).

Additionally, it is sometimes needed to use a standardized set of compilers, such as the ones distributed by conda. These instructions
demonstrate these two use cases.

```
$ cd $INSTALL; mkdir cctbx; cd cctbx
$ wget https://raw.githubusercontent.com/cctbx/cctbx_project/master/libtbx/auto_build/bootstrap.py
$ wget https://raw.githubusercontent.com/cctbx/cctbx_project/master/xfel/conda_envs/psana_environment.yml
$ python bootstrap.py --builder=xfel --use-conda=psana_environment.yml \
  --no-boost-src --python=39 hot update
$ source ~/miniconda3/etc/profile.d/conda.sh # modify as needed
$ conda activate base
$ mamba env create -f psana_environment.yml -p `pwd`/conda_base
$ conda activate `pwd`/conda_base
$ python bootstrap.py --builder=xfel --use-conda=psana_environment.yml \
  --config-flags="--compiler=conda" --config-flags="--use_environment_flags" \
  --config-flags="--no_bin_python"             \
  --no-boost-src --python=39 --nproc=10 build
$ source build/conda_setpaths.sh
$ source modules/cctbx_project/xfel/conda_envs/test_psana_lcls.sh
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

