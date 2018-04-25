from __future__ import absolute_import, division, print_function

import os
import sys

def run(args):
  try: __file__
  except NameError: d = sys.path[0]
  else: d = os.path.dirname(__file__)
  d,b = os.path.split(os.path.dirname(d))
  if b == "libtbx" and os.path.isdir(d):
    sys.path.insert(0, d)
  import libtbx.introspection
  n = libtbx.introspection.number_of_processors(return_value_if_unknown=None)
  if (n is not None or len(args) == 0):
    print(n)
  else:
    print(" ".join(args))

if __name__ == "__main__":
  run(args=sys.argv[1:])
