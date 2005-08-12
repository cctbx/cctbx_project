#! /bin/csh -f
set verbose
rm -f fortran_lbfgs.so
f77 -O2 -fast -c lbfgs.f
ld -r -o liblbfgs.a lbfgs.o
python pyfort.py -b -L . -l lbfgs fortran_lbfgs.pyf
find . -name fortran_lbfgs.so -exec mv {} "`libtbx.show_lib_path`" \;
rm -rf so_locations fortran_lbfgs.txt lbfgs.o liblbfgs.a build
