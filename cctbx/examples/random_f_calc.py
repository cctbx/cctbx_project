"""
usage: cctbx.python random_f_calc.py space_group_symbol

Example:
  cctbx.python random_f_calc.py P212121

This will write a file map_coeff.pickle that can be used to
view the electron density map with PyMOL (see view_fft_map.py).
"""
from __future__ import absolute_import, division, print_function

from cctbx.development import random_structure
from cctbx import sgtbx
from libtbx import easy_pickle
import sys

def generate_random_f_calc(space_group_info, n_elements=10, d_min=1.5):
  structure = random_structure.xray_structure(
    space_group_info,
    elements=["Si"]*n_elements,
    volume_per_atom=1000,
    min_distance=3.,
    general_positions_only=False)
  structure.show_summary().show_scatterers()
  print()
  f_calc = structure.structure_factors(
    d_min=d_min, anomalous_flag=False).f_calc()
  f_calc.show_summary()
  print()
  print("Writing file: map_coeff.pickle")
  easy_pickle.dump("map_coeff.pickle", f_calc)
  print()

def run():
  if (len(sys.argv) != 2):
    print(__doc__)
    return
  generate_random_f_calc(sgtbx.space_group_info(sys.argv[1]))

if (__name__ == "__main__"):
  run()
