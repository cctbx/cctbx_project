from __future__ import absolute_import, division, print_function
from shutil import rmtree
import sys, os

def run():
  for arg in sys.argv[1:]:
    if (os.path.isfile(arg)):
      os.remove(arg)
    elif (os.path.isdir(arg)):
      rmtree(arg)

if (__name__ == "__main__"):
  run()
