from libtbx.utils import Sorry
from libtbx.str_utils import show_string
import libtbx.load_env
import sys

def run(args):
  remaining_args = []
  prohibit_white_space = False
  for arg in args:
    if (arg == "--prohibit-white-space"):
      prohibit_white_space = True
    else:
      remaining_args.append(arg)
  args = remaining_args
  def process_key(target_key):
    n_hits = 0
    for line in open(libtbx.env.under_build(
                  "include_paths")).read().splitlines():
      flds = line.split(None, 1)
      assert len(flds) == 2
      source_key, path = flds
      if (target_key is None or source_key == target_key):
        if (prohibit_white_space and len(path.split()) != 1):
          raise Sorry(
            "Include path contains white-space: %s" % show_string(path))
        print flds[1]
        n_hits += 1
    if (n_hits == 0):
      raise Sorry("No such include path: %s" % show_string(arg))
  if (len(args) == 0):
    process_key(target_key=None)
  else:
    for arg in args:
      process_key(target_key=arg)

if (__name__ == "__main__"):
  run(sys.argv[1:])
