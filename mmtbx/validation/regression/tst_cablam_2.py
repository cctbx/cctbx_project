from __future__ import absolute_import, division, print_function
from libtbx import easy_run
from iotbx.pdb.fetch import fetch_and_write
from time import time
import libtbx.load_env

def exercise_01():
  fetch_and_write(id='5o61', entity='model_cif')
  fb = easy_run.fully_buffered("phenix.cablam 5o61.cif")
  assert fb.return_code == 0, fb.return_code

if (__name__ == "__main__"):
  t0 = time()
  if (not libtbx.env.has_module(name="phenix")):
    print("Skipping: probe not configured")
  else:
    exercise_01()
  print("Time: %.2f" % (time() - t0))
  print("OK")
