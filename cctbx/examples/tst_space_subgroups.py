from __future__ import absolute_import, division, print_function
from cctbx.examples import space_subgroups
import sys

def exercise(args):
  assert len(args) == 0
  space_subgroups.run(args=["1", "P1"])
  space_subgroups.run(args=["2", "P1"])
  space_subgroups.run(args=["2", "Pnnn:1"])
  space_subgroups.run(args=["2", "Fmmm"])
  print("OK")

if (__name__ == "__main__"):
  exercise(sys.argv[1:])
