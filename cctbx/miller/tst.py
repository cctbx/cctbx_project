import sys, os
import math
from cctbx_boost.arraytbx import shared
from cctbx_boost import uctbx
from cctbx_boost import sgtbx
from cctbx_boost import miller
from cctbx import xutils
from cctbx.development import debug_utils
random = debug_utils.random

def exercise_build_indices(SgInfo, index_abs_range = (6,6,6)):
  for friedel_flag in (1,0):
    miller_indices = miller.BuildIndices(SgInfo, friedel_flag, index_abs_range)
    miller_dict = {}
    for h in miller_indices: miller_dict[h] = 0
    sgops = SgInfo.SgOps()
    h = [0,0,0]
    for h[0] in range(-index_abs_range[0], index_abs_range[0]+1):
      for h[1] in range(-index_abs_range[1], index_abs_range[1]+1):
        for h[2] in range(-index_abs_range[2], index_abs_range[2]+1):
          if (sgops.isSysAbsent(h) or h == [0,0,0]): continue
          h_seq = miller.SymEquivIndices(sgops, h)
          found_h_asu = 0
          for i_eq in xrange(h_seq.M(friedel_flag)):
            h_eq = h_seq(i_eq).H()
            if (h_eq in miller_dict):
              assert found_h_asu == 0
              found_h_asu = 1
          assert found_h_asu != 0

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
  assert tuple(js.singles(1)) == ()
  assert tuple(js.pairs()) == ((0,2), (1,0), (2,3), (3,1))
  assert tuple(js.pair_selection(0)) == (1, 1, 1, 1, 0)
  assert tuple(js.single_selection(0)) == (0, 0, 0, 0, 1)
  assert tuple(js.pair_selection(1)) == (1, 1, 1, 1)
  assert tuple(js.single_selection(1)) == (0, 0, 0, 0)
  assert tuple(js.paired_miller_indices(0)) \
      == tuple(h0.select(js.pair_selection(0)))
  l1 = list(js.paired_miller_indices(1))
  l2 = list(h1.select(js.pair_selection(1)))
  l1.sort()
  l2.sort()
  assert l1 == l2
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
  miller_set.set_friedel_flag(0)
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

def get_random_structure(SgInfo, verbose=0):
  elements = ("N", "C", "C", "O", "N", "C", "C", "O")
  xtal = debug_utils.random_structure(
    SgInfo, elements,
    volume_per_atom=50.,
    min_distance=1.5,
    general_positions_only=0)
  if (0 or verbose):
    print "Unit cell:", xtal.UnitCell
    print "Space group:", xtal.SgInfo.BuildLookupSymbol()
  return xtal

def exercise_bins(SgInfo, n_bins=10, d_min=1, friedel_flag=1, verbose=0):
  xtal = get_random_structure(SgInfo, verbose)
  miller_set = xutils.build_miller_set(xtal, friedel_flag, d_min)
  fcalc_set = xutils.calculate_structure_factors_direct(
    miller_set, xtal, abs_F=1)
  binning1 = miller.binning(xtal.UnitCell, n_bins, miller_set.H)
  assert binning1.n_bins_used() == n_bins
  assert binning1.limits().size() == n_bins + 1
  assert binning1.n_bins_all() == n_bins + 2
  if (0 or verbose):
    print "binning1.d_max():", binning1.d_max()
    print "binning1.d_min():", binning1.d_min()
  binner1 = miller.binner(binning1, miller_set.H)
  if (0 or verbose): xutils.show_binner_summary(binner1)
  assert binner1.count(binner1.i_bin_d_too_large()) == 0
  assert binner1.count(binner1.i_bin_d_too_small()) == 0
  counts = binner1.counts()
  for i_bin in binner1.range_all():
    assert binner1.count(i_bin) == counts[i_bin]
    assert binner1(i_bin).count(1) == counts[i_bin]
  assert list(binner1.range_all()) == range(binner1.n_bins_all())
  assert list(binner1.range_used()) == range(1, binner1.n_bins_used()+1)
  binning2 = miller.binning(xtal.UnitCell, n_bins - 2,
    binning1.bin_d_min(2),
    binning1.bin_d_min(n_bins))
  binner2 = miller.binner(binning2, miller_set.H)
  if (0 or verbose): xutils.show_binner_summary(binner2)
  assert tuple(binner1.counts())[1:-1] == tuple(binner2.counts())
  array_indices = shared.size_t(tuple(xrange(miller_set.H.size())))
  perm_array_indices1 = shared.size_t()
  perm_array_indices2 = shared.size_t()
  for i_bin in binner1.range_all():
    perm_array_indices1.append(array_indices.select(binner1(i_bin)))
    perm_array_indices2.append(binner1.array_indices(i_bin))
  assert perm_array_indices1.size() == miller_set.H.size()
  assert perm_array_indices2.size() == miller_set.H.size()
  assert tuple(perm_array_indices1) == tuple(perm_array_indices2)
  assert tuple(perm_array_indices1.shuffle(perm_array_indices2)) \
      == tuple(array_indices)
  fcalc_set.setup_binner(reflections_per_bin=100)
  if (0 or verbose): xutils.show_binner_summary(fcalc_set.binner)
  for use_multiplicities in (0,1):
    filtered_set = []
    for negate in (0,1):
      filtered_set.append(fcalc_set.rms_filter(
        cutoff_factor=2,
        use_binning=1,
        use_multiplicities=use_multiplicities,
        negate=negate))
    assert filtered_set[0].H.size() > filtered_set[1].H.size()
    assert filtered_set[1].H.size() > 0
    js = miller.join_sets(filtered_set[0].H, filtered_set[1].H)
    assert js.pairs().size() == 0
    js0 = miller.join_sets(fcalc_set.H, filtered_set[0].H)
    assert js0.singles(1).size() == 0
    js1 = miller.join_sets(fcalc_set.H, filtered_set[1].H)
    assert js1.singles(1).size() == 0
    assert (js0.single_selection(0) | js1.single_selection(0)).count(0) == 0
    assert (js0.single_selection(0) & js1.single_selection(0)).count(1) == 0

def exercise_map_to_asu(SgInfo, d_min=2.5, friedel_flag=1, verbose=0):
  xtal = get_random_structure(SgInfo, verbose)
  miller_set = xutils.build_miller_set(xtal, friedel_flag, d_min)
  fcalc_set = xutils.calculate_structure_factors_direct(miller_set, xtal)
  fabs_set = xutils.reciprocal_space_array(fcalc_set, shared.abs(fcalc_set.F))
  phases = [shared.arg(fcalc_set.F, deg) for deg in (0,1)]
  hlc = shared.hendrickson_lattman()
  for i in fcalc_set.H.indices():
    hlc.append([random.random() for i in xrange(4)])
  h_random = shared.miller_Index()
  c_random = shared.complex_double()
  p_random = [shared.double(), shared.double()]
  hlc_random = shared.hendrickson_lattman()
  for i,h_asym in miller_set.H.items():
    h_seq = miller.SymEquivIndices(xtal.SgOps, h_asym)
    i_eq = random.randrange(h_seq.M(friedel_flag))
    h_eq = h_seq(i_eq)
    h_random.append(h_eq.H())
    c_random.append(h_eq.complex_eq(fcalc_set.F[i]))
    for deg in (0,1):
      p_random[deg].append(h_eq.phase_eq(phases[deg][i], deg))
    hlc_random.append(h_eq.hl_eq(hlc[i]))
  h_random_copy = h_random.deep_copy()
  miller.map_to_asu(SgInfo, friedel_flag, h_random_copy, c_random)
  for i,h_asym in miller_set.H.items():
    assert h_asym == h_random_copy[i]
  for i,f_asym in fcalc_set.F.items():
    assert abs(f_asym - c_random[i]) < 1.e-6
  h_random_copy = h_random.deep_copy()
  a_random = fabs_set.F.deep_copy()
  miller.map_to_asu(SgInfo, friedel_flag, h_random_copy, a_random)
  for i,h_asym in miller_set.H.items():
    assert h_asym == h_random_copy[i]
  for i,f_asym in fabs_set.F.items():
    assert f_asym == a_random[i]
  for deg in (0,1):
    h_random_copy = h_random.deep_copy()
    miller.map_to_asu(SgInfo, friedel_flag, h_random_copy, p_random[deg], deg)
    for i,h_asym in miller_set.H.items():
      assert h_asym == h_random_copy[i]
    for i,p_asym in phases[deg].items():
      assert debug_utils.phase_error(p_asym, p_random[deg][i], deg) < 1.e-5
  h_random_copy = h_random.deep_copy()
  miller.map_to_asu(SgInfo, friedel_flag, h_random_copy, hlc_random)
  for i,h_asym in miller_set.H.items():
    assert h_asym == h_random_copy[i]
  for i,hlc_asym in hlc.items():
    for j in xrange(4):
      assert abs(hlc_asym[j] - hlc_random[i][j]) < 1.e-5

def exercise_map_fft(SgInfo, d_min=2.5, verbose=0):
  xtal = get_random_structure(SgInfo, verbose)
  map_min = []
  map_max = []
  for friedel_flag in (0,1):
    miller_set = xutils.build_miller_set(xtal, friedel_flag, d_min)
    fcalc_set = xutils.calculate_structure_factors_direct(miller_set, xtal)
    map_set = xutils.fft_map(fcalc_set)
    s = shared.statistics(map_set.get_real_map())
    if (0 or verbose):
      print "friedel_flag=%d, map min: %f" % (map_set.friedel_flag, s.min())
      print "                map max: %f" % (s.max(),)
    map_min.append(s.min())
    map_max.append(s.max())
  assert abs(map_min[0] - map_min[1]) < 1.e-6
  assert abs(map_max[0] - map_max[1]) < 1.e-6

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
    exercise_build_indices(SgInfo)
    exercise_bins(SgInfo)
    exercise_map_to_asu(SgInfo)
    exercise_map_fft(SgInfo)
    sys.stdout.flush()

if (__name__ == "__main__"):
  run()
  t = os.times()
  print "u+s,u,s: %.2f %.2f %.2f" % (t[0] + t[1], t[0], t[1])
