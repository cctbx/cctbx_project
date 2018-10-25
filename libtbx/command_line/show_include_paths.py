from __future__ import absolute_import, division, print_function
from libtbx.utils import Sorry
from libtbx.str_utils import show_string
import libtbx.load_env
import sys

def read_include_paths(file_name="include_paths"):
  result = []
  for line in open(libtbx.env.under_build(file_name)).read().splitlines():
    flds = line.split(None, 1)
    assert len(flds) == 2
    result.append(flds)
  return result

def run(args):
  remaining_args = []
  prohibit_white_space = False
  for arg in args:
    if (arg == "--prohibit-white-space"):
      prohibit_white_space = True
    else:
      remaining_args.append(arg)
  include_paths = read_include_paths()
  if (len(remaining_args) == 0):
    for key,path in include_paths:
      print(key)
  else:
    for target_key in remaining_args:
      n_hits = 0
      for key,path in include_paths:
        if (key == target_key):
          if (prohibit_white_space and len(path.split()) != 1):
            raise Sorry(
              "Include path contains white-space: %s" % show_string(path))
          print(path)
          n_hits += 1
      if (n_hits == 0):
        raise Sorry("No such include path: %s" % show_string(arg))

if (__name__ == "__main__"):
  run(sys.argv[1:])
