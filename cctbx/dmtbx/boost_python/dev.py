from cctbx import dmtbx
from cctbx import miller
from cctbx.array_family import flex
from cctbx.development import random_structure
from cctbx.development import debug_utils
import math
import sys

def exercise(space_group_info, n_scatterers=8, d_min=2, verbose=0,
             e_min=1.5):
  structure = random_structure.xray_structure(
    space_group_info,
    elements=["const"]*n_scatterers,
    volume_per_atom=100,
    min_distance=2.,
    general_positions_only=0001,
    u_iso=0.0)
  if (0 or verbose):
    structure.show_summary().show_scatterers()
  f_calc = structure.structure_factors(
    d_min=d_min, anomalous_flag=00000).f_calc()
  f_obs = abs(f_calc)
  q_obs = miller.array(
    miller_set=f_obs,
    data=f_obs.data()
        / math.sqrt(f_obs.space_group().order_p() * n_scatterers)
        / f_obs.space_group().n_ltr())
  q_obs = q_obs.sort(by_value="abs")
  q_obs.setup_binner(auto_binning=True)
  n_obs = q_obs.normalize_structure_factors(quasi=0001)
  r = flex.linear_regression(q_obs.data(), n_obs.data())
  if (0 or verbose):
    r.show_summary()
  assert r.is_well_defined()
  assert abs(r.y_intercept()) < 0.1
  assert abs(r.slope() - 1) < 0.3
  q_large = q_obs.apply_selection(
    q_obs.quasi_normalized_as_normalized().data() > e_min)
  print "Number of e-values > %.6g: %d" % (e_min, q_large.size())
  print "slow:"
  slow_tprs = dmtbx.triplet_invariants(
    q_large.space_group_info().type(), q_large.indices(), 0001, 0001)
  #slow_tprs.dump_triplets(q_large.indices())
  print "fast:"
  fast_tprs = dmtbx.fast_triplets(
    q_large.space_group_info().type(), q_large.indices())
  #fast_tprs.dump_triplets(q_large.indices())
  slow_packed = slow_tprs.pack_triplets()
  assert slow_packed.size() % 7 == 0
  slow_packed.resize(flex.grid((slow_packed.size()/7,7)))
  print "slow grid:", slow_packed.all()
  fast_packed = fast_tprs.pack_triplets()
  assert fast_packed.size() % 7 == 0
  fast_packed.resize(flex.grid((fast_packed.size()/7,7)))
  print "fast grid:", fast_packed.all()
  if (slow_packed.all()[0] != fast_packed.all()[0]):
    print "LOOK Different number of triplets", q_large.space_group_info()
  qi = q_large.indices()
  for i in xrange(min(slow_packed.all()[0], fast_packed.all()[0])):
    sp = [slow_packed[(i,j)] for j in xrange(7)]
    #print "s", sp
    fp = [fast_packed[(i,j)] for j in xrange(7)]
    #print "f", fp
    if (sp[:6] != fp[:6]):
      print "LOOK Mismatch", q_large.space_group_info()
      break
    #if (sp[1] == sp[3]):
    #  print "EQUAL", qi[sp[0]], qi[sp[1]], qi[sp[3]], sp[6], fp[6]
    if (sp[6] != fp[6]):
      print "LOOK weights", q_large.space_group_info()
    assert sp == fp

def run_call_back(flags, space_group_info):
  e_min = 1.5
  if (flags.e_min != 00000):
    e_min = float(eval(flags.e_min))
  exercise(space_group_info, verbose=flags.Verbose,
    e_min=e_min)

def run():
  debug_utils.parse_options_loop_space_groups(sys.argv[1:], run_call_back, (
    "e_min",
  ))

if (__name__ == "__main__"):
  run()
