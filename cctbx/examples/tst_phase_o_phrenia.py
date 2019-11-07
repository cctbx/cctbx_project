from __future__ import absolute_import, division, print_function
from cctbx.examples import phase_o_phrenia
from cctbx.development import random_structure
from cctbx.development import debug_utils
import sys

def exercise(space_group_info, n_scatterers=1, d_min=2, verbose=0):
  structure = random_structure.xray_structure(
    space_group_info,
    elements=["Hg"]*n_scatterers,
    volume_per_atom=500,
    min_distance=2.,
    general_positions_only=True)
  if (1 or verbose):
    structure.show_summary().show_scatterers()
  reduced_peaks = phase_o_phrenia.calculate_exp_i_two_phi_peaks(
    structure, d_min, min_peak_distance=3, max_reduced_peaks=20)
  for peak in reduced_peaks:
    print("%.6g" % peak.height, "%8.5f %8.5f %8.5f" % peak.site)

def run_call_back(flags, space_group_info):
  exercise(space_group_info, verbose=flags.Verbose)

def run():
  debug_utils.parse_options_loop_space_groups(sys.argv[1:], run_call_back)

if (__name__ == "__main__"):
  run()
