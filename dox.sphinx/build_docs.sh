#!/bin/bash
# This script automatically creates the document tree for most of the cctbx modules.
# To keep your src directory clean it's best to copy the dox.sphinx directory into
# the build directory first:
# $> cd build
# $> cp -arf ../sources/cctbx_project/dox.sphinx .
# also it's necessary to call setpaths_all.sh before running this script:
# $> source setpaths_all.sh

MODULES="cctbx chiltbx cma_es crys3d fable fftw3tbx gltbx iotbx mmtbx omptbx rstbx scitbx smtbx spotfinder wxtbx xfel" #left out: cudatbx libtbx
mkdir $MODULES
for MODULE in $MODULES
do
  ./generate_modules.py --doc-header $MODULE -s rst -d ./$MODULE ../../sources/cctbx_project/$MODULE
done
make html -b coverage
