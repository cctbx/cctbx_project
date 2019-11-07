from __future__ import absolute_import, division, print_function
import fftw3tbx
import os

if self.env.is_ready_for_build():
  for libfftw3 in [fftw3tbx.libfftw3, fftw3tbx.libfftw3f]:
    source = self.env.under_build("base/lib/"+libfftw3)
    if os.path.isfile(source):
      target = self.env.under_build("lib/"+libfftw3)
      print("  Copying:", libfftw3)
      with open(target, "wb") as fh:
        with open(source, "rb") as rh:
          fh.write(rh.read())
