from __future__ import absolute_import, division, print_function
# XXX This is a workaround for the Boost floating-point error that is
# triggered when importing numpy (used in a variety of modules in CCTBX).
# Importing Numpy before any of the boost extensions avoids a crash.
try :
  import numpy
except ImportError :
  pass
import pydoc

def run():
  pydoc.cli()

if (__name__ == "__main__"):
  run()
