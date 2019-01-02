
from __future__ import absolute_import, division, print_function
import sys

if (__name__ == "__main__"):
  from libtbx.auto_build import create_windows_exe
  sys.exit(create_windows_exe.run(sys.argv[1:]))
