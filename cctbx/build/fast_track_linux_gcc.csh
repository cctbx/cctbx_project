#! /bin/csh -f
if ("`gcc --version`" !~ 2.95.* && "`gcc --version`" !~ 3.0) then
  if ("`gcc --version`" == 2.96) then
    echo "WARNING:"
    echo "  If you are running RedHat 7.0 then you need to update"
    echo "  glibc to version 2.2 and gcc to version 2.96-69."
    echo "  See http://www.redhat.com/support/errata/rh7-errata-bugfixes.html"
    echo "  for more information."
    echo ""
    echo "The installation will proceed in 10 seconds."
    echo ""
    sleep 10
  else
    echo "FATAL CONFIGURATION PROBLEM:"
    echo "  Available compiler: gcc version" `gcc --version`
    echo "  Required compiler:  gcc version 2.95.x, with x >= 2, or 2.96"
    echo "  See http://cctbx.sourceforge.net/page_compilers.html for more information."
    exit 1
  endif
endif
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
