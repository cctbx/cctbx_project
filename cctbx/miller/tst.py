import sys, os
import math
from cctbx_boost.arraytbx import shared
from cctbx_boost import uctbx
from cctbx_boost import sgtbx
from cctbx_boost import miller
from cctbx import xutils
from cctbx.development import debug_utils

def exercise_join_sets():
  h0 = shared.miller_Index(((1,2,3), (-1,-2,-3), (2,3,4), (-2,-3,-4), (3,4,5)))
  d0 = shared.double((1,2,3,4,5))
  h1 = shared.miller_Index(((-1,-2,-3), (-2,-3,-4), (1,2,3), (2,3,4)))
  d1 = shared.double((10,20,30,40))
  js = miller.join_sets(h0, h0)
  assert js.have_singles() == 0
  assert list(js.pairs()) == zip(range(5), range(5))
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

def show_binner_info(binner):
  for i_bin in binner.range_all():
    bin_d_range = binner.bin_d_range(i_bin)
    count = binner.count(i_bin)
    if (i_bin == binner.i_bin_d_too_large()):
      assert bin_d_range[0] == 0
      print "unused:              d > %8.4f: %5d" % (bin_d_range[1], count)
    elif (i_bin == binner.i_bin_d_too_small()):
      assert bin_d_range[1] == 0
      print "unused: %9.4f >  d           : %5d" % (bin_d_range[0], count)
    else:
      print "bin %2d: %9.4f >= d > %8.4f: %5d" % (
        (i_bin,) + bin_d_range + (count,))

def exercise_bins(SgInfo, n_bins=10, d_min=1):
  elements = ("N", "C", "C", "O", "N", "C", "C", "O")
  friedel_flag = 1
  xtal = debug_utils.random_structure(
    SgInfo, elements,
    volume_per_atom=50.,
    min_distance=1.5,
    general_positions_only=0)
  print "Unit cell:", xtal.UnitCell
  print "Space group:", xtal.SgInfo.BuildLookupSymbol()
  miller_set = xutils.build_miller_set(xtal, friedel_flag, d_min)
  binning1 = miller.binning(xtal.UnitCell, n_bins, miller_set.H)
  assert binning1.n_bins_used() == n_bins
  assert binning1.limits().size() == n_bins + 1
  assert binning1.n_bins_all() == n_bins + 2
  print "binning1.d_max():", binning1.d_max()
  print "binning1.d_min():", binning1.d_min()
  binner1 = miller.binner(binning1, miller_set.H)
  show_binner_info(binner1)
  assert binner1.count(binner1.i_bin_d_too_large()) == 0
  assert binner1.count(binner1.i_bin_d_too_small()) == 0
  counts = binner1.counts()
  for i_bin in binner1.range_all():
    assert binner1.count(i_bin) == counts[i_bin]
    assert binner1.bin_selection(i_bin).count(1) == counts[i_bin]
  assert list(binner1.range_all()) == range(binner1.n_bins_all())
  assert list(binner1.range_used()) == range(1, binner1.n_bins_used()+1)
  binning2 = miller.binning(xtal.UnitCell, n_bins - 2,
    binning1.bin_d_min(2),
    binning1.bin_d_min(n_bins))
  binner2 = miller.binner(binning2, miller_set.H)
  show_binner_info(binner2)
  assert tuple(binner1.counts())[1:-1] == tuple(binner2.counts())

def run():
  exercise_join_sets()
  Flags = debug_utils.command_line_options(sys.argv[1:], (
    "RandomSeed",
    "AllSpaceGroups",
  ))
  if (not Flags.RandomSeed): debug_utils.set_random_seed(0)
  symbols_to_stdout = 0
  auto_test = 0
  if (len(sys.argv) > 1 + Flags.n):
    symbols = sys.argv[1:]
  else:
    symbols = debug_utils.get_test_space_group_symbols(Flags.AllSpaceGroups)
    symbols_to_stdout = 1
  if (len(sys.argv) == 1 or (Flags.AllSpaceGroups and len(sys.argv) == 2)):
    auto_test = 1
  for RawSgSymbol in symbols:
    if (RawSgSymbol.startswith("--")): continue
    SgSymbols = sgtbx.SpaceGroupSymbols(RawSgSymbol)
    SgInfo = sgtbx.SpaceGroup(SgSymbols).Info()
    LookupSymbol = SgInfo.BuildLookupSymbol()
    sys.stdout.flush()
    print >> sys.stderr, LookupSymbol
    sys.stderr.flush()
    if (symbols_to_stdout):
      print LookupSymbol
      sys.stdout.flush()
    exercise_bins(SgInfo)
    sys.stdout.flush()

if (__name__ == "__main__"):
  run()
  t = os.times()
  print "u+s,u,s:", t[0] + t[1], t[0], t[1]
