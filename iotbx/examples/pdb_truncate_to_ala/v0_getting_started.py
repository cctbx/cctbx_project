from __future__ import absolute_import, division, print_function
import iotbx.pdb
import sys

def run(args):
  if (len(args) == 0):
    raise RuntimeError("Please specify one or more pdb file names.")
  for file_name in args:
    pdb_obj = iotbx.pdb.input(file_name=file_name)
    pdb_obj.construct_hierarchy().overall_counts().show()

if (__name__ == "__main__"):
  run(sys.argv[1:])
