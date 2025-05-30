"""Convert SHELX to mmCIF format"""
from __future__ import absolute_import, division, print_function
from cctbx import xray
import sys, os
op = os.path

def run(args):
  for f in args:
    try:
      xs = xray.structure.from_shelx(filename=f, strictly_shelxl=False)
    except KeyboardInterrupt:
      raise
    except Exception:
      print("%s is not a .ins or a .res file" % f)
      continue
    r, _ = op.splitext(op.basename(f))
    xs.as_cif_simple(out=open(r + '.cif', 'w'))

if (__name__ == "__main__"):
  run(args=sys.argv[1:])
