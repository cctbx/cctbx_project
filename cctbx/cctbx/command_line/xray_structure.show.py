from cctbx import xray
from libtbx import easy_pickle
import sys

def run(args):
  for file_name in args:
    structures = easy_pickle.load(file_name)
    if (isinstance(structures, xray.structure)):
      structures = [structures]
    for structure in structures:
      structure.show_summary().show_scatterers()
      print

if (__name__ == "__main__"):
  run(sys.argv[1:])
