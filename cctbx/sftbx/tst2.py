import sys, os
from cctbx_boost import sgtbx
from cctbx_boost import adptbx
from cctbx import xutils
from cctbx.development import debug_utils

def print_structure_factors(SgInfo, adp=0, d_min=3.):
  Elements = ("N", "C", "C", "O", "N", "C", "C", "O")
  xtal = debug_utils.random_structure(
    SgInfo, Elements,
    volume_per_atom=50.,
    min_distance=1.5,
    general_positions_only=0,
    anisotropic_displacement_parameters=adp)
  xsym = xutils.crystal_symmetry(xtal.UnitCell, xtal.SgInfo)
  SymSites = xutils.symmetrized_sites(xsym, xtal.Sites)
  debug_utils.print_sites(SymSites)
  MillerIndices = xutils.miller_index_set(xsym, d_min)
  Fcalc = xutils.structure_factors(MillerIndices, SymSites)
  for i in xrange(len(Fcalc.H)):
    print Fcalc.H[i], "%.3g %.0f" % xutils.f_as_ampl_phase(Fcalc.F[i])

def run():
  Flags = debug_utils.command_line_options(sys.argv[1:], (
    "RandomSeed",
    "AllSpaceGroups",
    "Anisotropic",
  ))
  if (not Flags.RandomSeed): debug_utils.set_random_seed(0)
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
    print_structure_factors(SgInfo, Flags.Anisotropic)
    sys.stdout.flush()

if (__name__ == "__main__"):
  run()
  t = os.times()
  print "u+s,u,s:", t[0] + t[1], t[0], t[1]
