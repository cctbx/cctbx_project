from cctbx import find_centre_of_inversion
from cctbx import xray
from cctbx.development import random_structure
from cctbx.development import debug_utils
import math
import random
import sys

def exercise_core(structure, t_centre_of_inversion):
  refined_matches = find_centre_of_inversion.find_matches(structure)
  assert len(refined_matches) > 0
  assert len(refined_matches[0].pairs) == structure.scatterers().size()/2
  t_ctr = refined_matches[0].t_ctr
  for i,o in zip(t_centre_of_inversion, t_ctr):
    d = math.fmod(i-o, 1)
    tolerance = 1.e-6
    if (d < 0): d += 1
    if (d > tolerance): d -= 1
    assert d <= tolerance

def exercise(space_group_info, n_scatterers=8, verbose=0):
  structure = random_structure.xray_structure(
    space_group_info,
    n_scatterers=n_scatterers,
    volume_per_atom=500,
    min_distance=2.,
    general_positions_only=0001)
  t_centre_of_inversion = tuple([random.random() for i in xrange(3)])
  if (0 or verbose):
    print "t_centre_of_inversion: %.4f %.4f %.4f" % t_centre_of_inversion
  structure.build_scatterers(
    ["Hg"]*n_scatterers,
    t_centre_of_inversion=t_centre_of_inversion)
  if (0 or verbose):
    structure.show_summary().show_scatterers()
  exercise_core(structure, t_centre_of_inversion)
  mod_structure = xray.structure(
    special_position_settings=structure,
    scatterers=structure.scatterers()[:-1])
  if (0 or verbose):
    mod_structure.show_summary().show_scatterers()
  exercise_core(mod_structure, t_centre_of_inversion)
  mod_structure = structure.random_modify_parmeters("site", gauss_sigma=0.1)
  if (0 or verbose):
    mod_structure.show_summary().show_scatterers()
  refined_matches = find_centre_of_inversion.find_matches(mod_structure)
  assert len(refined_matches) > 0
  assert len(refined_matches[0].pairs) == n_scatterers/2

def run_call_back(flags, space_group_info):
  exercise(
    space_group_info,
    verbose=flags.Verbose)

def run():
  debug_utils.parse_options_loop_space_groups(sys.argv[1:], run_call_back)

if (__name__ == "__main__"):
  run()
