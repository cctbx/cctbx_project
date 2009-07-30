import scitbx.rigid_body.essence.tst_tardy
from scitbx.graph import tst_tardy_pdb
from scitbx.array_family import flex
from scitbx import matrix
from libtbx.test_utils import \
  approx_equal, is_above_limit, is_below_limit, Exception_expected
from libtbx.utils import format_cpu_times
import random
import sys

def compare_essence_and_fast_tardy_models(etm):
  ftm = scitbx.rigid_body.tardy_model(
    labels=etm.labels,
    sites=flex.vec3_double(etm.sites),
    masses=flex.double(etm.masses),
    tardy_tree=etm.tardy_tree,
    potential_obj=etm.potential_obj,
    near_singular_hinges_angular_tolerance_deg=
      etm.near_singular_hinges_angular_tolerance_deg)
  assert ftm.bodies_size() == len(etm.bodies)
  assert ftm.number_of_trees == etm.number_of_trees
  assert ftm.degrees_of_freedom == etm.degrees_of_freedom
  assert ftm.packed_q_size == etm.packed_q_size
  assert ftm.packed_qd_size == etm.packed_qd_size
  ftm.flag_positions_as_changed()
  ftm.flag_velocities_as_changed()
  #
  assert list(ftm.root_indices()) == etm.root_indices()
  #
  packed_q_orig = etm.pack_q()
  packed_qd_orig = etm.pack_qd()
  #
  def check_packed():
    e = etm.pack_q()
    f = ftm.pack_q()
    assert approx_equal(e, f)
    e = etm.pack_qd()
    f = ftm.pack_qd()
    assert approx_equal(e, f)
  check_packed()
  mt = flex.mersenne_twister(seed=0)
  packed_q_rand = mt.random_double(size=packed_q_orig.size())*2-1
  packed_qd_rand = mt.random_double(size=packed_qd_orig.size())*2-1
  for tm in [etm, ftm]: tm.unpack_q(packed_q=packed_q_rand)
  check_packed()
  for tm in [etm, ftm]: tm.unpack_q(packed_q=packed_q_orig)
  check_packed()
  for tm in [etm, ftm]: tm.unpack_qd(packed_qd=packed_qd_rand)
  check_packed()
  for tm in [etm, ftm]: tm.unpack_qd(packed_qd=packed_qd_orig)
  check_packed()
  for tm in [etm, ftm]: tm.unpack_q(packed_q=packed_q_rand)
  check_packed()
  for tm in [etm, ftm]: tm.unpack_qd(packed_qd=packed_qd_rand)
  check_packed()
  #
  fnosiet = ftm.number_of_sites_in_each_tree()
  assert list(fnosiet) == etm.number_of_sites_in_each_tree()
  e = etm.sum_of_masses_in_each_tree()
  f = ftm.sum_of_masses_in_each_tree()
  assert approx_equal(e, f)
  #
  def xxx_check_spatial_inertia():
    e = [body.i_spatial for body in etm.bodies]
    #print "e:", [sum(ei) for ei in e]
    f = ftm.xxx_spatial_inertia()
    #print "f:", [sum(fi) for fi in f]
    for ei,fi in zip(e, f):
      assert approx_equal(ei, fi)
  xxx_check_spatial_inertia()
  #
  def check_mean_linear_velocity():
    e = etm.mean_linear_velocity(number_of_sites_in_each_tree=None)
    f = ftm.mean_linear_velocity(number_of_sites_in_each_tree=None)
    assert approx_equal(e, f)
    f = ftm.mean_linear_velocity(number_of_sites_in_each_tree=fnosiet)
    assert approx_equal(e, f)
  check_mean_linear_velocity()
  value = matrix.col(mt.random_double(size=3)*2-1)
  for tm in [etm, ftm]:
    tm.subtract_from_linear_velocities(
      number_of_sites_in_each_tree=None,
      value=value)
  check_mean_linear_velocity()
  for tm in [etm, ftm]:
    tm.subtract_from_linear_velocities(
      number_of_sites_in_each_tree=fnosiet,
      value=value)
  check_mean_linear_velocity()
  #
  e = etm.sites_moved()
  f = ftm.sites_moved()
  assert approx_equal(e, f)
  e = etm.e_pot()
  f = ftm.e_pot()
  assert approx_equal(e, f)
  e = etm.d_e_pot_d_sites()
  f = ftm.d_e_pot_d_sites()
  assert approx_equal(e, f)
  e = etm.d_e_pot_d_q_packed()
  f = ftm.d_e_pot_d_q_packed()
  e_max_abs = flex.max(flex.abs(e))
  if (etm.potential_obj is None):
    assert is_below_limit(value=e_max_abs, limit=1e-10)
  else:
    assert is_above_limit(value=e_max_abs, limit=1) # random generator dependent
  assert approx_equal(e, f)
  #
  e = etm.e_kin()
  f = ftm.e_kin()
  assert approx_equal(e, f)
  e = etm.e_tot()
  f = ftm.e_tot()
  assert approx_equal(e, f)
  etm.reset_e_kin(e_kin_target=1.234)
  ftm.reset_e_kin(e_kin_target=1.234)
  e = etm.e_kin()
  assert approx_equal(e, 1.234)
  f = ftm.e_kin()
  assert approx_equal(e, f)
  try:
    ftm.reset_e_kin(e_kin_target=1, e_kin_epsilon=0)
  except RuntimeError, e:
    assert str(e).find("e_kin_epsilon > 0") > 0
  else: raise Exception_expected
  etm.assign_zero_velocities()
  ftm.assign_zero_velocities()
  assert etm.e_kin() == 0
  assert ftm.e_kin() == 0
  random.seed(0)
  etm.assign_random_velocities()
  random.seed(0)
  ftm.assign_random_velocities()
  e = etm.featherstone_system_model().spatial_velocities()
  f = ftm.xxx_spatial_velocities()
  assert approx_equal(e, f)
  e = etm.e_kin()
  f = ftm.e_kin()
  xxx_check_spatial_inertia()
  assert approx_equal(e, f)
  #
  etm.unpack_q(packed_q=packed_q_orig)
  etm.unpack_qd(packed_qd=packed_qd_orig)

def run(args):
  assert len(args) == 0
  for tc in tst_tardy_pdb.test_cases:
    tt = tc.tardy_tree_construct()
    masses = [1.0]*len(tc.sites)
    etm = scitbx.rigid_body.essence.tardy.model(
      labels=tc.labels,
      sites=tc.sites,
      masses=masses,
      tardy_tree=tt,
      potential_obj=None)
    compare_essence_and_fast_tardy_models(etm=etm)
    etm = scitbx.rigid_body.essence.tst_tardy.construct_tardy_model(
      labels=tc.labels,
      sites=tc.sites,
      masses=masses,
      tardy_tree=tt)
    assert etm.potential_obj is not None
    compare_essence_and_fast_tardy_models(etm=etm)
  print format_cpu_times()

if (__name__ == "__main__"):
  run(args=sys.argv[1:])
