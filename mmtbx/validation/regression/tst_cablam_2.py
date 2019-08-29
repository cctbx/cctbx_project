from __future__ import absolute_import, division, print_function
from libtbx import easy_run
from time import time
import libtbx.load_env

def exercise_01():
  easy_run.call("phenix.fetch_pdb 5o61 -c")
  assert not easy_run.call("phenix.cablam 5o61.cif"), "For Chris to fix"

if (__name__ == "__main__"):
  t0 = time()
  if (not libtbx.env.has_module(name="phenix")):
    print("Skipping: probe not configured")
  else:
    exercise_01()
  print("Time: %.2f" % (time() - t0))
  print("OK")
