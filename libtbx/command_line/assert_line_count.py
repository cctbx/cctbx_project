from __future__ import absolute_import, division, print_function
import sys

def run(args):
  assert len(args) == 1
  expected_count = int(args[0])
  assert len(sys.stdin.read().splitlines()) == expected_count

if (__name__ == "__main__"):
  run(sys.argv[1:])
