from boost.rational import gcd
from libtbx.utils import Usage
import libtbx.load_env
import sys

def run(args, func=gcd):
  if (len(args) == 0):
    raise Usage("%s integer [integer...]" % libtbx.env.dispatcher_name)
  values = [int(arg) for arg in args]
  result = values[0]
  for value in values[1:]:
    result = func(result, value)
  print result

if (__name__ == "__main__"):
  run(sys.argv[1:])
