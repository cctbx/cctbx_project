import sys, os
from cctbx import xutils
from cctbx import statistics
from cctbx_boost import sgtbx
from cctbx_boost import adptbx
from cctbx.development import debug_utils
from cctbx.misc import python_utils

def exercise(SgInfo, d_min=1.0, reflections_per_bin=200, n_bins=10, verbose=0):
  elements = ("N", "C", "C", "O") * 5
  friedel_flag = 1
  xtal = debug_utils.random_structure(
    SgInfo, elements,
    volume_per_atom=50.,
    min_distance=1.5,
    general_positions_only=0,
    no_random_u=1,
    uiso=adptbx.B_as_U(10))
  if (0 or verbose):
    print "Unit cell:", xtal.UnitCell
    print "Space group:", xtal.SgInfo.BuildLookupSymbol()
    debug_utils.print_sites(xtal)
  miller_set = xutils.build_miller_set(xtal, friedel_flag, d_min)
  fcalc_set = xutils.calculate_structure_factors(miller_set, xtal, abs_F=1)
  fcalc_set.setup_binner(
    auto_binning=1, reflections_per_bin=reflections_per_bin, n_bins=n_bins)
  if (0 or verbose):
    xutils.show_binner_summary(fcalc_set.binner)
  asu_contents = python_utils.dict_with_default_value(0)
  for elem in elements: asu_contents[elem] += 1
  wp = statistics.wilson_plot(fcalc_set, asu_contents)
  print "wilson_k, wilson_b:", wp.wilson_k, wp.wilson_b
  assert 0.6 < wp.wilson_k < 1.4
  assert 9 < wp.wilson_b < 11
  if (0):
    f = open("tmp.xy", "w")
    sys.stdout = f
    wp.dump_raw_xy()
    sys.stdout = sys.__stdout__
    f.close()

def run():
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
    exercise(SgInfo)
    sys.stdout.flush()
  t = os.times()
  print "u+s,u,s: %.2f %.2f %.2f" % (t[0] + t[1], t[0], t[1])

if (__name__ == "__main__"):
  run()
