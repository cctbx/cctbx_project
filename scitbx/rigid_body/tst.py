import scitbx.rigid_body.essence.tardy
from scitbx.graph import tst_tardy_pdb
from scitbx.array_family import flex
from libtbx.test_utils import approx_equal
from libtbx.utils import format_cpu_times
from libtbx import Auto
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
  emlv = etm.mean_linear_velocity(number_of_sites_in_each_tree=Auto)
  assert approx_equal(
    ftm.mean_linear_velocity(number_of_sites_in_each_tree=None), emlv)
  assert approx_equal(
    ftm.mean_linear_velocity(number_of_sites_in_each_tree=fnosiet), emlv)
  ftm.subtract_from_linear_velocities(
    number_of_sites_in_each_tree=None,
    value=(0,0,0))
  ftm.subtract_from_linear_velocities(
    number_of_sites_in_each_tree=fnosiet,
    value=(0,0,0))

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
