import sys

def run(args):
  import os.path as op
  dn = op.dirname
  try: __file__
  except NameError: d = sys.path[0]
  else: d = dn(__file__)
  d,b = op.split(dn(d))
  if (b == "libtbx" and op.isdir(d)):
    sys.path.insert(0, d)
  import libtbx.introspection
  n = libtbx.introspection.number_of_processors(return_value_if_unknown=None)
  if (n is not None or len(args) == 0):
    print n
  else:
    print " ".join(args)

if (__name__ == "__main__"):
  run(args=sys.argv[1:])
