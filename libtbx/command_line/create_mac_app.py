
from __future__ import division
import sys

if (__name__ == "__main__") :
  from libtbx.auto_build import create_mac_app
  sys.exit(create_mac_app.run(sys.argv[1:]))
