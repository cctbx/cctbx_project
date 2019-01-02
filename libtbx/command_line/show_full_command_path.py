from __future__ import absolute_import, division, print_function
import libtbx.path
import sys

def run(args):
  assert len(args) == 1
  print(libtbx.path.full_command_path(command=args[0]))

if __name__ == "__main__" :
  run(sys.argv[1:])
