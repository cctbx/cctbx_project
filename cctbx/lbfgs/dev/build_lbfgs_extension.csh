#!/bin/csh -f
rm -f lbfgs.so
echo '*** creating lbfgs.o'
f77 -O2 -c lbfgs.f
echo '*** creating liblbfgs.a'
ld -r -o liblbfgs.a lbfgs.o
echo '*** calling pyfort'
python pyfort.py -b -L . -l lbfgs lbfgs.pyf
find . -name lbfgs.so -exec mv {} . \;
