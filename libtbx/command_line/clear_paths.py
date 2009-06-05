from libtbx.clear_paths \
  import remove_or_rename_files_and_directories_if_possible
import sys

def run(args):
  remaining = remove_or_rename_files_and_directories_if_possible(paths=args)
  for path in remaining:
    "WARNING: unable to remove or rename:", path

if (__name__ == "__main__"):
  run(args=sys.argv[1:])
