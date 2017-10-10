
from __future__ import division
from __future__ import print_function
from builtins import zip
from cctbx.development import random_structure
from cctbx.development import debug_utils
from cctbx.sgtbx import space_group_info
from scitbx.array_family import flex
from libtbx.test_utils import approx_equal, Exception_expected
from libtbx.utils import Sorry
import random
import math
import sys

if (1): # fixed random seed to avoid rare failures
  random.seed(0)
  flex.set_random_seed(0)

def exercise_ellipsoidal_truncation(space_group_info, n_sites=100, d_min=1.5):
  xrs = random_structure.xray_structure(
    space_group_info=space_group_info,
    elements=(("O","N","C")*(n_sites//3+1))[:n_sites],
    volume_per_atom=50,
    min_distance=1.5)
  f_obs = abs(xrs.structure_factors(d_min = d_min).f_calc())
  # exercise reciprocal_space_vector()
  for mi, d in zip(f_obs.indices(), f_obs.d_spacings().data()):
    rsv = flex.double(f_obs.unit_cell().reciprocal_space_vector(mi))
    assert approx_equal(d, 1./math.sqrt(rsv.dot(rsv)))
  ##
  print(f_obs.unit_cell())
  f = flex.random_double(f_obs.data().size())*flex.mean(f_obs.data())/10
  #
  f_obs1 = f_obs.customized_copy(data = f_obs.data(), sigmas = f_obs.data()*f)
  print("datat in:",f_obs1.data().size())
  r = f_obs1.ellipsoidal_truncation_by_sigma(sigma_cutoff=1)
  print("data left:",r.data().size())
  r.miller_indices_as_pdb_file(file_name="indices1.pdb", expand_to_p1=False)
  r.miller_indices_as_pdb_file(file_name="indices2.pdb", expand_to_p1=True)
  #
  f_obs.miller_indices_as_pdb_file(file_name="indices3.pdb", expand_to_p1=False)
  f_obs.miller_indices_as_pdb_file(file_name="indices4.pdb", expand_to_p1=True)
  print("*"*25)

def run_call_back(flags, space_group_info):
  exercise_ellipsoidal_truncation(space_group_info)

# TODO ideally this should loop over all possible space groups and all
# possible twin laws
def exercise_twinning () :
  xrs = random_structure.xray_structure(
    unit_cell=(12,5,12,90,90,90),
    space_group_symbol="P1",
    n_scatterers=12,
    elements="random")
  fc = abs(xrs.structure_factors(d_min=1.5).f_calc())
  fc = fc.set_observation_type_xray_amplitude()
  fc_twin = fc.twin_data("l,-k,h", 0.3)
  ic = fc.f_as_f_sq()
  fc_twin_2 = ic.twin_data("l,-k,h", 0.3).f_sq_as_f()
  assert (fc_twin.data().all_eq(fc_twin_2.data()))
  fc_tmp, fc_twin_tmp = fc.common_sets(other=fc_twin)
  # XXX in this particular crystal symmetry, a subset of reflections where h==l
  # will have the same value in the twinned and untwinned data - need to check
  # whether this is correct
  assert not fc_twin_tmp.data().all_approx_equal(fc_tmp.data())
  try :
    fc_twin = fc.twin_data("k,h,l", 0.5)
  except Sorry as s:
    pass
  else :
    raise Exception_expected
  fc_detwin = fc_twin.detwin_data("l,-k,h", 0.3)
  fc_detwin, fc = fc_detwin.common_sets(other=fc)
  assert fc_detwin.data().all_approx_equal(fc.data())
  # derived from PDB 3hfg; this confirms that the change in unit cell
  # parameters does not crash the routines
  xrs = random_structure.xray_structure(
    unit_cell=(56.438, 152.670, 74.203, 90.00, 92.41, 90.00),
    space_group_symbol="P21",
    n_scatterers=100,
    elements="random")
  fc = abs(xrs.structure_factors(d_min=1.5).f_calc())
  fc = fc.set_observation_type_xray_amplitude()
  fc_twin = fc.twin_data("h,-k,-l", 0.3)
  ic = fc.f_as_f_sq()
  fc_twin_2 = ic.twin_data("h,-k,-l", 0.3).f_sq_as_f()
  assert (fc_twin.data().all_eq(fc_twin_2.data()))
  fc_detwin = fc_twin.detwin_data("h,-k,-l", 0.3)
  fc_detwin, fc = fc_detwin.common_sets(other=fc)
  assert fc_detwin.data().all_approx_equal(fc.data())

if (__name__ == "__main__"):
  exercise_twinning()
  debug_utils.parse_options_loop_space_groups(sys.argv[1:], run_call_back)
