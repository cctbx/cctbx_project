#! /usr/bin/env python

from iotbx import crystal_symmetry_from_any
from iotbx import format
import sys

def run():
  assert len(sys.argv) == 2
  crystal_symmetry = crystal_symmetry_from_any.extract_from(sys.argv[1])
  if (crystal_symmetry is None):
    raise RuntimeError, \
      "Unknown file format or unit cell and/or space group missing from file."
  format.crystal_symmetry(crystal_symmetry)

if (__name__ == "__main__"):
  run()
