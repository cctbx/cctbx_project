#!/bin/csh -f
rm -f fortran_lbfgs.so
echo '*** creating lbfgs.o'
f77 -O2 -fast -c lbfgs.f
echo '*** creating liblbfgs.a'
ld -r -o liblbfgs.a lbfgs.o
echo '*** calling pyfort'
python pyfort.py -b -L . -l lbfgs fortran_lbfgs.pyf
find . -name fortran_lbfgs.so -exec mv {} $LIBTBX_BUILD/libtbx \;
rm -rf so_locations fortran_lbfgs.txt lbfgs.o liblbfgs.a build
