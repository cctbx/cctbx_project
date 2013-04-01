
from __future__ import division
from libtbx.auto_build import create_windows_exe
import sys

if (__name__ == "__main__") :
  sys.exit(create_windows_exe.run(sys.argv[1:]))
