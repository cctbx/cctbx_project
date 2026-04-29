"""Show paths to all components of Phenix"""
from __future__ import absolute_import, division, print_function
from six.moves import range

def run(args):
  import libtbx.load_env
  remaining_args = []
  dirname_count = 0
  for arg in args:
    if (arg == "--dirname"):
      dirname_count += 1
    else:
      remaining_args.append(arg)
  def show(path):
    if (path is not None):
      from os.path import dirname
      for _ in range(dirname_count):
        path = dirname(path)
    print(path)
  if remaining_args:
    for arg in remaining_args:
      show(libtbx.env.dist_path(module_name=arg, default=None))
  else:
    for path in libtbx.env.dist_paths():
      show(path)

if __name__ == "__main__":
  import sys
  run(args=sys.argv[1:])

