from __future__ import absolute_import, division, print_function
import cctbx.array_family.flex # import dependency
import boost_adaptbx.boost.python as bp
ext = bp.import_ext("mmtbx_ncs_ext")
from mmtbx_ncs_ext import *
from libtbx.test_utils import approx_equal

def esercise_00():
  p = pair(
    r = ([1,2,3,4,5,6,7,8,9]),
    t = ([10,11,12]),
    radius=13,
    radius_estimate=13.5,
    fracscat=15,
    rho_mn=([1,2,3,4,5]),
    id=0)
  assert approx_equal(p.r,[1,2,3,4,5,6,7,8,9])
  assert approx_equal(p.t,[10,11,12])
  assert approx_equal(p.radius, 13)
  assert approx_equal(p.radius_estimate, 13.5)
  assert approx_equal(p.fracscat, 15)
  assert approx_equal(p.rho_mn,[1,2,3,4,5])
  # exercise set
  p.set_radius(10)
  p.set_rhoMN([5,4,3,2,1])
  assert p.id==0
  p.set_id(9)
  assert p.id==9
  assert approx_equal(p.radius, 10)
  assert approx_equal(p.rho_mn, [5,4,3,2,1])

if (__name__ == "__main__"):
  esercise_00()
