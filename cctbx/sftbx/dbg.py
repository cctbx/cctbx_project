import sys, os
from cctbx_boost import sgtbx
from cctbx_boost import adptbx
from cctbx_boost import sftbx # XXX
from cctbx import xutils
from cctbx.development import debug_utils

def print_structure_factors(SgInfo, adp=0, d_min=3.):
  elements = ("N", "C", "C", "O", "N", "C", "C", "O")
  xtal = debug_utils.random_structure(
    SgInfo, elements,
    volume_per_atom=50.,
    min_distance=1.5,
    general_positions_only=0,
    anisotropic_displacement_parameters=adp)
  print xtal.UnitCell
  debug_utils.print_sites(xtal)
  d_min = 1.
  max_q = 1. / (d_min**2)
  resolution_factor = 1./3
  max_prime = 5
  mandatory_grid_factors = xtal.SgOps.refine_gridding()
  sampled_density = sftbx.sample_model_density(
    xtal.UnitCell, xtal.Sites,
    max_q, resolution_factor, max_prime, mandatory_grid_factors)
  print "max_q:", sampled_density.max_q()
  print "resolution_factor:", sampled_density.resolution_factor()
  print "quality_factor:", sampled_density.quality_factor()
  print "wing_cutoff:", sampled_density.wing_cutoff()
  print "u_extra:", sampled_density.u_extra()

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
