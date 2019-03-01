from __future__ import absolute_import, division, print_function
import libtbx.load_env

def run():
  print("# Listing of commands in: %s" % abs(libtbx.env.bin_path))
  file_names = sorted(libtbx.env.bin_path.listdir())
  print("# Number of commands: %d" % len(file_names))
  for file_name in file_names:
    print(file_name)

if (__name__ == "__main__"):
  run()
