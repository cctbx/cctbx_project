from __future__ import division
from __future__ import print_function
import libtbx.load_env

def run():
  print(abs(libtbx.env.bin_path))

if (__name__ == "__main__"):
  run()
