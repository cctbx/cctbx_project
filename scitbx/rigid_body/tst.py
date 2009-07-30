import scitbx.rigid_body.essence.tardy
from scitbx.graph import tst_tardy_pdb
from scitbx.array_family import flex
from libtbx.test_utils import approx_equal, Exception_expected
from libtbx.utils import format_cpu_times
from libtbx import Auto
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
  ftm.flag_positions_as_changed()
  ftm.flag_velocities_as_changed()
  assert list(ftm.root_indices()) == etm.root_indices()
  fnosiet = ftm.number_of_sites_in_each_tree()
  assert list(fnosiet) == etm.number_of_sites_in_each_tree()
  assert approx_equal(
    ftm.sum_of_masses_in_each_tree(),
    etm.sum_of_masses_in_each_tree())
  #
  def xxx_check_spatial_inertia():
    e = [body.i_spatial for body in etm.bodies]
    #print "e:", [sum(ei) for ei in e]
    f = ftm.xxx_spatial_inertia()
    #print "f:", [sum(fi) for fi in f]
    for ei,fi in zip(e, f):
      assert approx_equal(fi, ei)
  xxx_check_spatial_inertia()
  #
  e = etm.mean_linear_velocity(number_of_sites_in_each_tree=Auto)
  assert approx_equal(
    ftm.mean_linear_velocity(number_of_sites_in_each_tree=None), e)
  assert approx_equal(
    ftm.mean_linear_velocity(number_of_sites_in_each_tree=fnosiet), e)
  ftm.subtract_from_linear_velocities(
    number_of_sites_in_each_tree=None,
    value=(0,0,0))
  ftm.subtract_from_linear_velocities(
    number_of_sites_in_each_tree=fnosiet,
    value=(0,0,0))
  e = etm.sites_moved()
  f = ftm.sites_moved()
  assert approx_equal(f, e)
  e = etm.e_pot()
  f = ftm.e_pot()
  assert approx_equal(f, e)
  e = etm.d_e_pot_d_sites()
  f = ftm.d_e_pot_d_sites()
  assert approx_equal(f, e)
  e = etm.e_kin()
  f = ftm.e_kin()
  assert approx_equal(f, e)
  e = etm.e_tot()
  f = ftm.e_tot()
  assert approx_equal(f, e)
  etm.reset_e_kin(e_kin_target=1)
  ftm.reset_e_kin(e_kin_target=1)
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
  e = etm.pack_q()
  f = ftm.pack_q()
  assert approx_equal(f, e)
  e = etm.pack_qd()
  f = ftm.pack_qd()
  assert approx_equal(f, e)
  assert len(e) == len(f)
  e = etm.featherstone_system_model().spatial_velocities()
  f = ftm.xxx_spatial_velocities()
  assert approx_equal(f, e)
  e = etm.e_kin()
  f = ftm.e_kin()
  xxx_check_spatial_inertia()
  assert approx_equal(f, e)

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
  print format_cpu_times()

if (__name__ == "__main__"):
  run(args=sys.argv[1:])
