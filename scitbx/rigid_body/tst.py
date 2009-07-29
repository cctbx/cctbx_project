import scitbx.rigid_body.essence.tardy
from scitbx.graph import tst_tardy_pdb
from scitbx.array_family import flex
from libtbx.utils import format_cpu_times
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
