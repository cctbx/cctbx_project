#! /bin/csh -f
set try_python = "`python -c 'print __name__'`"
if ($status != 0 || "$try_python" !~ __main__) then
  echo "FATAL CONFIGURATION PROBLEM:"
  echo "  python not available."
  echo "  See http://cctbx.sourceforge.net/page_installation.html for more information."
  exit 1
endif
set PACKAGES="$0"
set PACKAGES="$PACKAGES:h"
set PACKAGES="$PACKAGES:h"
set PACKAGES="$PACKAGES:h"
mkdir tru64_cxx
cd tru64_cxx
mkdir boost
cd boost
cp $PACKAGES/boost/libs/python/build/tru64_cxx.mak Makefile
make ROOT=$PACKAGES softlinks
make ROOT=$PACKAGES
make test
cd ..
mkdir cctbx
cd cctbx
cp $PACKAGES/cctbx/build/configuration_tru64_cxx configuration
python $PACKAGES/cctbx/build/boot.py
python make.py softlinks
python make.py compile_all
python test.py
examples/cpp/getting_started
source setpythonpath.csh
python $PACKAGES/cctbx/examples/python/getting_started.py
echo "Done."
