import sys, os, math
from cctbx_boost.arraytbx import shared
from cctbx_boost import sgtbx
from cctbx import xutils
from cctbx.development import debug_utils
from cctbx_boost import dmtbx

def inplace_divide(array, arg):
  for i in xrange(array.size()):
    array[i] = array[i] / arg

def erase_small(miller_indices, data, cutoff):
  assert miller_indices.size() == data.size()
  for i in xrange(miller_indices.size()-1, -1, -1):
    if (data[i] < cutoff):
      miller_indices.erase(i)
      data.erase(i)

def exercise(SgInfo,
             number_of_point_atoms = 10,
             d_min=1.,
             e_min=1.2):
  elements = ["const"] * number_of_point_atoms
  xtal = debug_utils.random_structure(
    SgInfo, elements,
    volume_per_atom=50.,
    min_distance=1.5,
    general_positions_only=0,
    no_random_u=1)
  print xtal.UnitCell
  debug_utils.print_sites(xtal)
  MillerIndices = xutils.build_miller_indices(xtal, friedel_flag=1,d_min=d_min)
  Fcalc = xutils.calculate_structure_factors(MillerIndices, xtal, abs_F=1)
  inplace_divide(Fcalc.F, math.sqrt(number_of_point_atoms))
  dmtbx.inplace_sort(MillerIndices.H, Fcalc.F, 1)
  s = shared.statistics(Fcalc.F)
  print "mean2:", s.mean2()
  print "number of structure factors:", Fcalc.F.size()
  erase_small(MillerIndices.H, Fcalc.F, e_min)
  print "number of structure factors:", Fcalc.F.size()
  debug_utils.print_structure_factors(Fcalc)
  t = dmtbx.triplet_invariants(xtal.SgInfo, MillerIndices.H, Fcalc.F)

def run():
  Flags = debug_utils.command_line_options(sys.argv[1:], (
    "RandomSeed",
    "AllSpaceGroups",
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
    exercise(SgInfo)
    sys.stdout.flush()

if (__name__ == "__main__"):
  run()
  t = os.times()
  print "u+s,u,s:", t[0] + t[1], t[0], t[1]
