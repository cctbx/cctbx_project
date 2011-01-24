from cctbx import xray
from libtbx import easy_pickle
import sys, os

for f in sys.argv[1:]:
  try:
    xs = xray.structure.from_cif(file_path=f)
  except KeyboardInterrupt:
    raise
  except:
    print "%s is not a .cif file" % f
    continue
  r, _ = os.path.splitext(f)
  easy_pickle.dump(r + '.pickle', xs)
