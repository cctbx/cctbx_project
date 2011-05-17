import sys
from scitbx.array_family import flex
from cctbx.sgtbx import space_group_info
from cctbx.development import random_structure
from libtbx.test_utils import approx_equal
from cctbx.development import debug_utils
import random
import math
import scitbx.matrix

if (1): # fixed random seed to avoid rare failures
  random.seed(0)
  flex.set_random_seed(0)

def exercise(space_group_info, n_sites=100, d_min=1.5):
  from cctbx import maptbx
  from cctbx.masks import vdw_radii_from_xray_structure
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
  print f_obs.unit_cell()
  f = flex.random_double(f_obs.data().size())*flex.mean(f_obs.data())/10
  #
  f_obs1 = f_obs.customized_copy(data = f_obs.data(), sigmas = f_obs.data()*f)
  print "datat in:",f_obs1.data().size()
  r = f_obs1.ellipsoidal_truncation_by_sigma(sigma_cutoff=1)
  print "data left:",r.data().size()
  r.miller_indices_as_pdb_file(file_name="indices1.pdb", expand_to_p1=False)
  r.miller_indices_as_pdb_file(file_name="indices2.pdb", expand_to_p1=True)
  #
  f_obs.miller_indices_as_pdb_file(file_name="indices3.pdb", expand_to_p1=False)
  f_obs.miller_indices_as_pdb_file(file_name="indices4.pdb", expand_to_p1=True)
  print "*"*25

def run_call_back(flags, space_group_info):
  exercise(space_group_info)

if (__name__ == "__main__"):
  debug_utils.parse_options_loop_space_groups(sys.argv[1:], run_call_back)
