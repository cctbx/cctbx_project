import sys, os
from cctbx_boost.arraytbx import shared
from cctbx_boost import sgtbx
from cctbx_boost import adptbx
from cctbx_boost import fftbx
from cctbx_boost import sftbx
from cctbx import xutils
from cctbx.development import debug_utils

def exercise_map_collect(xtal, friedel_flag, fcalc, conjugate):
  max_q = xtal.UnitCell.max_Q(fcalc.H)
  grid_resolution_factor = 1/3.
  max_prime = 5
  mandatory_grid_factors = xtal.SgOps.refine_gridding()
  grid_logical = sftbx.determine_grid(
    xtal.UnitCell,
    max_q, grid_resolution_factor, max_prime, mandatory_grid_factors)
  rfft = fftbx.real_to_complex_3d(grid_logical)
  if (friedel_flag):
    n_complex = rfft.Ncomplex()
  else:
    cfft = fftbx.complex_to_complex_3d(grid_logical)
    n_complex = cfft.N()
  map = sftbx.structure_factor_map(
    xtal.SgOps, friedel_flag, fcalc.H, fcalc.F, n_complex, conjugate)
  if (friedel_flag):
    rfft.backward(map)
    rfft.forward(map)
  else:
    cfft.backward(map)
    cfft.forward(map)
  miller2, fcalc2 = sftbx.collect_structure_factors(
    xtal.UnitCell, xtal.SgInfo, friedel_flag,
    max_q, map, n_complex, conjugate)
  js = shared.join_sets(fcalc.H, miller2)
  assert js.pairs().size() == fcalc.H.size()
  debug_utils.show_structure_factor_correlation(
    "map/collect1", fcalc.H, js, fcalc.F, fcalc2,
    min_corr_ampl=0.9999, max_mean_w_phase_error=.01,
    verbose=0)
  fcalc2 = sftbx.collect_structure_factors(
    friedel_flag, fcalc.H, map, n_complex, conjugate)
  debug_utils.show_structure_factor_correlation(
    "map/collect2", fcalc.H, 0, fcalc.F, fcalc2,
    min_corr_ampl=0.9999, max_mean_w_phase_error=.01,
    verbose=0)

def print_structure_factors(SgInfo,
                            adp=0,
                            random_f_prime=0,
                            random_f_double_prime=0,
                            d_min=3.):
  elements = ("N", "C", "C", "O", "N", "C", "C", "O")
  grid_resolution_factor = 1./3
  if (random_f_prime): random_f_prime = grid_resolution_factor * d_min
  friedel_flag = not random_f_double_prime
  xtal = debug_utils.random_structure(
    SgInfo, elements,
    volume_per_atom=50.,
    min_distance=1.5,
    general_positions_only=0,
    random_f_prime_d_min=random_f_prime,
    random_f_double_prime=random_f_double_prime,
    anisotropic_displacement_parameters=adp)
  print xtal.UnitCell
  debug_utils.print_sites(xtal)
  MillerIndices = xutils.build_miller_indices(xtal, friedel_flag, d_min)
  Fcalc = xutils.calculate_structure_factors(MillerIndices, xtal)
  debug_utils.print_structure_factors(Fcalc)
  fcalc_fft = xutils.calculate_structure_factors_fft(MillerIndices, xtal)
  debug_utils.show_structure_factor_correlation(
    "direct/fft", Fcalc.H, 0, Fcalc.F, fcalc_fft.F,
    min_corr_ampl=0.99, max_mean_w_phase_error=3.)
  for i in xrange(2):
    exercise_map_collect(xtal, friedel_flag, Fcalc, i)

def run():
  Flags = debug_utils.command_line_options(sys.argv[1:], (
    "RandomSeed",
    "AllSpaceGroups",
    "Isotropic",
    "Anisotropic",
    "Dispersive",
    "Anomalous",
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
    if (not auto_test):
      print_structure_factors(SgInfo,
        random_f_prime=Flags.Dispersive,
        random_f_double_prime=Flags.Anomalous,
        adp=Flags.Anisotropic)
    else:
      for fdp in xrange(2):
        for adp in xrange(2):
          print_structure_factors(SgInfo,
            random_f_prime=1,
            random_f_double_prime=fdp,
            adp=adp)
    sys.stdout.flush()

if (__name__ == "__main__"):
  run()
  t = os.times()
  print "u+s,u,s:", t[0] + t[1], t[0], t[1]
