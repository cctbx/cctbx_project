import sys, os
from cctbx_boost.arraytbx import shared
from cctbx_boost import sgtbx
from cctbx_boost import adptbx
from cctbx_boost import sftbx
from cctbx_boost import fftbx
from cctbx import xutils
from cctbx.development import debug_utils

def print_structure_factors(SgInfo, adp=0, d_min=1.):
  elements = ("N", "C", "C", "O", "N", "C", "C", "O")
  xtal = debug_utils.random_structure(
    SgInfo, elements,
    volume_per_atom=400.,
    min_distance=1.5,
    general_positions_only=0,
    anisotropic_displacement_parameters=adp)
  if (0):
    from cctbx_boost import uctbx
    from cctbx_boost.eltbx.caasf_wk1995 import CAASF_WK1995
    xtal.UnitCell = uctbx.UnitCell((10,10,10))
    xtal.Sites = shared.XrayScatterer()
    #s = sftbx.XrayScatterer("O", CAASF_WK1995("O"), 0j, (1./64,1./64,1./64), 1., 0)
    s = sftbx.XrayScatterer("O", CAASF_WK1995("O"), 0j, (1./32,1./16,1./6), 1., 0.5)
    s.ApplySymmetry(xtal.UnitCell, xtal.SgOps, 0.1, 0, 0)
    xtal.Sites.append(s)
    s = sftbx.XrayScatterer("O", CAASF_WK1995("O"), 0j, (8./32,5./16,3./6), 1., 0.5)
    s.ApplySymmetry(xtal.UnitCell, xtal.SgOps, 0.1, 0, 0)
    xtal.Sites.append(s)
  print xtal.UnitCell
  debug_utils.print_sites(xtal)
  MillerIndices = xutils.build_miller_indices(xtal, d_min)
  Fcalc = xutils.calculate_structure_factors(MillerIndices, xtal)
  max_q = 1. / (d_min**2)
  resolution_factor = 1./3
  max_prime = 5
  mandatory_grid_factors = xtal.SgOps.refine_gridding()
  grid_logical = sftbx.determine_grid(
    xtal.UnitCell,
    max_q, resolution_factor, max_prime, mandatory_grid_factors)
  fft = fftbx.real_to_complex_3d(grid_logical)
  sampled_density = sftbx.sampled_model_density(
    xtal.UnitCell, xtal.Sites,
    max_q, resolution_factor,
    fft.Nreal(), fft.Mreal())
  print "max_q:", sampled_density.max_q()
  print "resolution_factor:", sampled_density.resolution_factor()
  print "quality_factor:", sampled_density.quality_factor()
  print "wing_cutoff:", sampled_density.wing_cutoff()
  print "u_extra:", sampled_density.u_extra()
  map = sampled_density.map_as_shared()
  map_stats = shared.statistics(map)
  print "Electron density"
  print "max %.6g" % (map_stats.max())
  print "min %.6g" % (map_stats.min())
  print "mean %.6g" % (map_stats.mean())
  print "sigma %.6g" % (map_stats.sigma())
  if (0):
    map = sftbx.structure_factor_map(
      xtal.SgOps, Fcalc.H, Fcalc.F, fft.Ncomplex())
    fft.backward(map)
    map_stats = shared.statistics(map)
    print "True electron density"
    print "max %.6g" % (map_stats.max())
    print "min %.6g" % (map_stats.min())
    print "mean %.6g" % (map_stats.mean())
    print "sigma %.6g" % (map_stats.sigma())
  fft.forward(map)
  map_stats = shared.statistics(map)
  print "Structure factors"
  print "max %.6g" % (map_stats.max())
  print "min %.6g" % (map_stats.min())
  print "mean %.6g" % (map_stats.mean())
  print "sigma %.6g" % (map_stats.sigma())
  print "Ncomplex", fft.Ncomplex()
  if (0):
    map = sftbx.structure_factor_map(
      xtal.SgOps, Fcalc.H, Fcalc.F, fft.Ncomplex())
  miller_indices, fcal = sftbx.collect_structure_factors(
    xtal.UnitCell, xtal.SgInfo,
    max_q, map, fft.Ncomplex())
  js = shared.join_sets(MillerIndices.H, miller_indices)
  if (0):
    for i,j in js.pairs():
      print MillerIndices.H[i], miller_indices[j]
    print "singles 1:"
    for i in js.singles(0):
      print MillerIndices.H[i]
    print "singles 2:"
    for i in js.singles(1):
      print miller_indices[i]
  assert js.pairs().size() + js.singles(0).size() == MillerIndices.H.size()
  assert js.pairs().size() + js.singles(1).size() == miller_indices.size()
  assert js.pairs().size() == MillerIndices.H.size()
  for i in xrange(2):
    assert js.singles(i).size() == 0
  x = shared.double()
  y = shared.double()
  for i,j in js.pairs():
    assert MillerIndices.H[i] == miller_indices[j]
    print MillerIndices.H[i],
    print "(" + debug_utils.format_structure_factor(Fcalc.F[i]) + ")",
    print "(" + debug_utils.format_structure_factor(fcal[j]) + ")"
    x.append(abs(Fcalc.F[i]))
    y.append(abs(fcal[j]))
  xy_regr = shared.linear_regression(x, y)
  assert xy_regr.is_well_defined()
  print "cc:", xy_regr.cc()
  print "m:", xy_regr.m()
  if (0):
    x = list(x.as_tuple())
    x.sort()
    x = shared.double(tuple(x))
    y = list(y.as_tuple())
    y.sort()
    y = shared.double(tuple(y))
    xy_regr = shared.linear_regression(x, y)
    assert xy_regr.is_well_defined()
    print "cc:", xy_regr.cc()
    print "m:", xy_regr.m()

def run():
  Flags = debug_utils.command_line_options(sys.argv[1:], (
    "RandomSeed",
    "AllSpaceGroups",
    "Isotropic",
    "Anisotropic",
  ))
  if (not Flags.RandomSeed): debug_utils.set_random_seed(0)
  if (not (Flags.Isotropic or Flags.Anisotropic)):
    Flags.Isotropic = 1
    # XXX Flags.Anisotropic = 1
  symbols_to_stdout = 0
  if (len(sys.argv) > 1 + Flags.n):
    symbols = sys.argv[1:]
  else:
    symbols = debug_utils.get_test_space_group_symbols(Flags.AllSpaceGroups)
    symbols_to_stdout = 1
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
    if (Flags.Isotropic):
      print_structure_factors(SgInfo, adp=0)
    if (Flags.Anisotropic):
      print_structure_factors(SgInfo, adp=1)
    sys.stdout.flush()

if (__name__ == "__main__"):
  run()
  t = os.times()
  print "u+s,u,s:", t[0] + t[1], t[0], t[1]
