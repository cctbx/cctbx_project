from __future__ import division
import cctbx.array_family.flex # import dependency
import boost.python
ext = boost.python.import_ext("mmtbx_ncs_ext")
from mmtbx_ncs_ext import *
from libtbx.test_utils import approx_equal

def esercise_00():
  p = pair(
    r = ([1,2,3,4,5,6,7,8,9]),
    t = ([10,11,12]),
    radius=13,
    weight=14,
    fracscat=15,
    rho_mn=([1,2,3,4,5,6,7,8,9,10,11,12,13,14,15]))
  assert approx_equal(p.r,[1,2,3,4,5,6,7,8,9])
  assert approx_equal(p.t,[10,11,12])
  assert approx_equal(p.radius, 13)
  assert approx_equal(p.weight, 14)
  assert approx_equal(p.fracscat, 15)
  assert approx_equal(p.rho_mn,[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15])

if (__name__ == "__main__"):
  esercise_00()
