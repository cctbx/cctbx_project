import math
from cctbx_boost.arraytbx import shared
from cctbx_boost import uctbx
from cctbx_boost import sgtbx
from cctbx_boost import miller
from cctbx import xutils
h0 = shared.miller_Index(((1,2,3), (-1,-2,-3), (2,3,4), (-2,-3,-4), (3,4,5)))
d0 = shared.double((1,2,3,4,5))
h1 = shared.miller_Index(((-1,-2,-3), (-2,-3,-4), (1,2,3), (2,3,4)))
d1 = shared.double((10,20,30,40))
js = miller.join_sets(h0, h1)
assert tuple(js.singles(0)) == (4,)
assert tuple(js.pairs()) == ((0,2), (1,0), (2,3), (3,1))
assert tuple(js.plus(d0, d1)) == (31, 12, 43, 24)
assert tuple(js.minus(d0, d1)) == (-29,-8,-37,-16)
assert tuple(js.multiplies(d0, d1)) == (30,20,120,80)
assert tuple(js.divides(d0, d1)) == (1/30.,2/10.,3/40.,4/20.)
assert list(js.additive_sigmas(d0, d1)) == [
  math.sqrt(x*x+y*y) for x,y in ((1,30), (2,10), (3,40), (4,20))]
jbm = miller.join_bijvoet_mates(h0)
assert tuple(jbm.singles()) == (4,)
assert tuple(jbm.pairs()) == ((0,1), (2,3))
assert tuple(jbm.minus(d0)) == (-1, -1)
assert tuple(jbm.average(d0)) == (3/2., 7/2.)
assert list(jbm.additive_sigmas(d0)) == [
  math.sqrt(x*x+y*y) for x,y in ((1,2), (3,4))]
ucell = uctbx.UnitCell()
sginfo = sgtbx.SpaceGroup().Info()
xtal = xutils.crystal_symmetry(ucell, sginfo)
miller_set = xutils.miller_set(xtal, h0)
data_set0 = xutils.reciprocal_space_array(miller_set, d0, d0)
anom_diffs = data_set0.anomalous_differences()
assert tuple(anom_diffs.H) == ((1,2,3), (2,3,4))
assert tuple(anom_diffs.F) == (-1, -1)
assert list(anom_diffs.sigmas) == [
  math.sqrt(x*x+y*y) for x,y in ((1,2), (3,4))]
miller_set = xutils.miller_set(xtal, h1)
data_set1 = xutils.reciprocal_space_array(miller_set, d1, d1)
sum = data_set0 + 3
assert sum.H == data_set0.H
assert tuple(sum.F) == ((4,5,6,7,8))
assert tuple(sum.sigmas) == ((4,5,6,7,8))
sum = data_set0 + data_set1
assert tuple(sum.H) == ((1,2,3), (-1,-2,-3), (2,3,4), (-2,-3,-4))
assert tuple(sum.F) == (31, 12, 43, 24)
assert list(sum.sigmas) == [
  math.sqrt(x*x+y*y) for x,y in ((1,30), (2,10), (3,40), (4,20))]
selected_data_set0 = data_set0.sigma_filter(cutoff_factor=2)
assert tuple(selected_data_set0.H) == ()
selected_data_set0 = data_set0.sigma_filter(cutoff_factor=2, negate=1)
assert tuple(selected_data_set0.H) == tuple(h0)
assert tuple(selected_data_set0.F) == tuple(d0)
data_set0.sigmas = shared.double((0.6, 0.4, 2., 0.5, 0.5))
selected_data_set0 = data_set0.sigma_filter(cutoff_factor=2)
assert tuple(selected_data_set0.H) == ((-1, -2, -3), (-2, -3, -4), (3, 4, 5))
assert tuple(selected_data_set0.F) == (2,4,5)
assert tuple(selected_data_set0.sigmas) == (0.4, 0.5, 0.5)
selected_data_set0 = data_set0.sigma_filter(cutoff_factor=2, negate=1)
assert tuple(selected_data_set0.H) == ((1, 2, 3), (2, 3, 4))
assert tuple(selected_data_set0.F) == (1,3)
assert tuple(selected_data_set0.sigmas) == (0.6, 2.)
data_set = xutils.reciprocal_space_array(
  xutils.miller_set(xtal, shared.miller_Index(((1,2,3), (2,3,4)))),
  shared.double((0,1)))
selected_data_set = data_set.rms_filter(cutoff_factor=0.5)
assert tuple(selected_data_set.H) == ((1,2,3),)
assert tuple(selected_data_set.F) == (0,)
data_set.friedel_flag = 0
selected_data_set = data_set.rms_filter(
  cutoff_factor=0.5, use_multiplicities=1)
assert tuple(selected_data_set.H) == ((1,2,3),)
assert tuple(selected_data_set.F) == (0,)
print "OK"
