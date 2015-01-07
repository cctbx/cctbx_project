from __future__ import division
import sys

if (__name__ == "__main__") :
  from libtbx.auto_build import create_installer
  sys.exit(create_installer.run(sys.argv[1:]))
