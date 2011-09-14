from libtbx.math_utils import prime_factors_of
import sys

def run(args):
  for arg in args:
    n = int(arg)
    assert n > 0
    print "prime factors of %d:" % n, prime_factors_of(n)

if (__name__ == "__main__"):
  run(args=sys.argv[1:])
