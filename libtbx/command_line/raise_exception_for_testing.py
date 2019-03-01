from __future__ import absolute_import, division, print_function
from libtbx.utils import show_exception_info_if_full_testing
import libtbx.load_env
import sys

def run(args):
  assert args in [[], ["silent"]]
  if (len(args) == 0):
    libtbx.env.full_testing = True
  else:
    libtbx.env.full_testing = False
  try:
    raise RuntimeError("Just for testing.")
  except RuntimeError:
    show_exception_info_if_full_testing()

if (__name__ == "__main__"):
  run(sys.argv[1:])
