#! /bin/csh -f
tst3d fftw $1 $2 $3 $4 -1 > zw
tst3d fftpack $1 $2 $3 $4 -1 > zx
paste zw zx | std-linear-correlation
