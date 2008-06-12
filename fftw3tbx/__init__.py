try:
  import scitbx.array_family.flex
  import boost.python
  ext = boost.python.import_ext("fftw3tbx_ext")
except ImportError:
  ext = None
if (ext is not None):
  from fftw3tbx_ext import *

import sys

fftw3_h = "fftw3.h"

if (sys.platform.startswith("darwin")):
  libfftw3 = "libfftw3.dylib"
else:
  libfftw3 = "libfftw3.so.3"
