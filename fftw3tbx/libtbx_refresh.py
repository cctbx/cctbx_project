import fftw3tbx
from libtbx import easy_run
import sys, os

if (self.env.is_ready_for_build()):
  source = self.env.under_build("base/lib/"+fftw3tbx.libfftw3)
  if (os.path.isfile(source)):
    target = self.env.under_build("lib/"+fftw3tbx.libfftw3)
    print "Copying:", fftw3tbx.libfftw3
    open(target, "wb").write(open(source, "rb").read())
