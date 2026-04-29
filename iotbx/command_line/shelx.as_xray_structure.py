"""Create xray_structure from a SHELX file and dump as pickle"""
from __future__ import absolute_import, division, print_function
from cctbx import xray
from libtbx import easy_pickle
import sys, os

for f in sys.argv[1:]:
  try:
    xs = xray.structure.from_shelx(filename=f, strictly_shelxl=False)
  except KeyboardInterrupt:
    raise
  except Exception:
    print("%s is not a .ins or a .res file" % f)
    continue
  r, _ = os.path.splitext(f)
  easy_pickle.dump(r + '.pickle', xs)
