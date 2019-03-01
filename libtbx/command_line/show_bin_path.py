from __future__ import absolute_import, division, print_function
import libtbx.load_env

def run():
  print(abs(libtbx.env.bin_path))

if (__name__ == "__main__"):
  run()
