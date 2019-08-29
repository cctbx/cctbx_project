from __future__ import absolute_import, division, print_function
from cctbx import sgtbx
import sys

def run(args):
  for symbol in args:
    space_group_info = sgtbx.space_group_info(symbol=symbol)
    space_group_info.show_summary()
    for m in space_group_info.group():
      rt = m.as_rational().as_float()
      print(" ", list(rt.r) + list(rt.t))
    print()

if (__name__ == "__main__"):
  run(sys.argv[1:])
