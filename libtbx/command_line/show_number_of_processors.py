import libtbx.introspection
import sys

def run(args):
  n = libtbx.introspection.number_of_processors(return_value_if_unknown=None)
  if (n is not None or len(args) == 0):
    print n
  else:
    print " ".join(args)

if (__name__ == "__main__"):
  run(args=sys.argv[1:])
