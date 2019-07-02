from __future__ import absolute_import, division, print_function
from scitbx.array_family import flex
import sys

def run(args):
  assert len(args) == 0
  print("Size of an empty container.")
  print("Some of those containers store a pointer to a raw storage:")
  print("the 'cumulative' entry refers to the total memory footprint")
  print("whereas the other entry does not include the raw storage.")
  print()
  flex.show_sizes_int()
  print()
  flex.show_sizes_double()
  print()

if (__name__ == "__main__"):
  run(args=sys.argv[1:])
