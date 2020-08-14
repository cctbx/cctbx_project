from __future__ import absolute_import, division, print_function
import gcd
from boost_adaptbx.boost.rational import lcm
import sys

def run(args):
  gcd.run(args=args, func=lcm)

if (__name__ == "__main__"):
  run(sys.argv[1:])
