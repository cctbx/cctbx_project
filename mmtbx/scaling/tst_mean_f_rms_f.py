
from __future__ import absolute_import, division, print_function
import sys
#from libtbx.test_utils import approx_equal
from libtbx.utils import null_out

def exercise(args):
  if 'verbose' in args:
    out=sys.stdout
  else:
    out=null_out()

  print("Testing mean_f_rms_f....", end=' ')
  from mmtbx.scaling.mean_f_rms_f import test
  test(out=out)

if (__name__ == "__main__"):
  exercise(sys.argv[1:])
  print("OK")
