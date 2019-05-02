# simtbx
Simulation Toolbox - for simulating X-ray diffraction images ab
initio.

Note to developers.  This package requires additional testing with an external github repo, LS49.  The test should be performed on a Linux system with at least 16 cores available.

Instructions are as follows (2 May 2019). Use a development build of cctbx + dials:
```
cd ${DEVELOPER_DIR} # the normal location for bootstrap.py and the modules directory for source code
python bootstrap.py hot update --builder=dials
cd modules
git clone git@github.com:nksauter/LS49.git

cd ${BUILD_LS49} # create and change to an empty build directory; use conda if appropriate
python ../modules/cctbx_project/libtbx/configure.py LS49 prime iota --enable_openmp_if_possible=True --use_conda
source setpaths.sh
make
export OMP_NUM_THREADS=16 # or higher

cd ${DATA_DIR} # create an empty directory and download big data file from Nick Sauter
tar xzf ls49.tgz
export LS49_BIG_DATA=${DATA_DIR}/ls49_big_data

cd ${TEST} # create an empty directory for test execution
rm *
libtbx.run_tests_parallel module=LS49 nproc=15
rm *
libtbx.run_tests_parallel module=simtbx nproc=5
```
