#! /bin/csh -f
if ("`gcc --version`" !~ 2.95.*) then
  echo "FATAL CONFIGURATION PROBLEM:"
  echo "  Available compiler: gcc version" `gcc --version`
  echo "  Required compiler:  gcc version 2.95.x, with x >= 2"
  echo "  See http://cctbx.sourceforge.net/page_compilers.html for more information."
  exit 1
endif
set PACKAGES="$0"
set PACKAGES="$PACKAGES:h"
set PACKAGES="$PACKAGES:h"
set PACKAGES="$PACKAGES:h"
mkdir linux_gcc
cd linux_gcc
mkdir boost
cd boost
cp $PACKAGES/boost/libs/python/build/linux_gcc.mak Makefile
make ROOT=$PACKAGES cp # with softlinks make test fails on some machines
make ROOT=$PACKAGES
make test
cd ..
mkdir cctbx
cd cctbx
cp $PACKAGES/cctbx/build/configuration_linux_gcc configuration
python $PACKAGES/cctbx/build/boot.py
python make.py softlinks
python make.py compile_all
python test.py
examples/cpp/getting_started
source setpythonpath.csh
python $PACKAGES/cctbx/examples/python/getting_started.py
echo "Done."
