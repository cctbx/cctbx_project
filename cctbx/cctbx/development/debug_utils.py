import math
import random
from cctbx_boost.arraytbx import shared
from cctbx_boost import uctbx
from cctbx_boost import sgtbx
from cctbx_boost.eltbx.caasf_wk1995 import CAASF_WK1995
from cctbx_boost import adptbx
from cctbx_boost import sftbx
from cctbx import xutils
from cctbx.macro_mol import rotation_parameters

def set_random_seed(seed):
  random.seed(seed)

class command_line_options:

  def __init__(self, argv, keywords):
    self.n = 0
    for keyword in keywords:
      is_keyword = ("--" + keyword in argv)
      setattr(self, keyword, is_keyword)
      if (is_keyword):
        self.n += 1

def get_test_space_group_symbols(FlagAllSpaceGroups):
  if (FlagAllSpaceGroups):
    sgnumbers = xrange(1, 231)
  else:
    sgnumbers = (1,2,3,15,16,74,75,142,143,167,168,194,195,230)
  return [sgtbx.SpaceGroupSymbols(sgno).ExtendedHermann_Mauguin()
          for sgno in sgnumbers]

def get_compatible_unit_cell(SgInfo, Volume):
  if   (SgInfo.SgNumber() <   3):
    uc = uctbx.UnitCell((1., 1.3, 1.7, 83, 109, 129))
  elif (SgInfo.SgNumber() <  15):
    uc = uctbx.UnitCell((1., 1.3, 1.7, 90, 109, 90))
  elif (SgInfo.SgNumber() <  75):
    uc = uctbx.UnitCell((1., 1.3, 1.7, 90, 90, 90))
  elif (SgInfo.SgNumber() < 143):
    uc = uctbx.UnitCell((1., 1., 1.7, 90, 90, 90))
  elif (SgInfo.SgNumber() < 195):
    uc = uctbx.UnitCell((1., 1., 1.7, 90, 90, 120))
  else:
    uc = uctbx.UnitCell((1., 1., 1., 90, 90, 90))
  uc = uc.ChangeBasis(SgInfo.CBOp().M().as_tuple()[0])
  f = math.pow(Volume / uc.getVolume(), 1/3.)
  p = list(uc.getParameters())
  for i in xrange(3): p[i] *= f
  uc = uctbx.UnitCell(p)
  SgInfo.SgOps().CheckUnitCell(uc)
  return uc

def generate_positions(UnitCell, SgOps, N,
                       min_mate_distance = 1.5,
                       min_hetero_distance = 1.5,
                       general_positions_only = 0):
  SnapParameters = sgtbx.SpecialPositionSnapParameters(
    UnitCell, SgOps, 1, min_mate_distance)
  max_back_track = 100
  for back_track in xrange(max_back_track):
    Positions = []
    for i in xrange(N):
      max_placement_trials = 100
      have_position = 0
      for t in xrange(max_placement_trials):
        X = (random.random(), random.random(), random.random())
        SS = sgtbx.SiteSymmetry(SnapParameters, X)
        if (general_positions_only and SS.M() != SgOps.OrderZ()): continue
        SS.expand()
        sec = sgtbx.SymEquivCoordinates(SS)
        min_d2 = UnitCell.getLongestVector2()
        for p in Positions:
          min_d2 = min(min_d2, sec.getShortestDistance2(UnitCell, p))
        if (min_d2 >= min_hetero_distance**2):
          have_position = 1
          break
      if (not have_position): break
      Positions.append(SS.SnapPosition())
    if (len(Positions) == N): break
  assert len(Positions) == N, "Cannot find position matching all constraints."
  return Positions

def get_site_symmetry(UnitCell, SgOps, X,
                      min_mate_distance = 3.):
  SnapParameters = sgtbx.SpecialPositionSnapParameters(
    UnitCell, SgOps, 1, min_mate_distance)
  return sgtbx.SiteSymmetry(SnapParameters, X)

def perturb_coordinates(ucell, sgops, x, sigma, vary_z_only = 0):
  SSx = get_site_symmetry(ucell, sgops, x)
  xc = ucell.orthogonalize(x)
  max_placement_trials = 100
  have_position = 0
  for t in xrange(max_placement_trials):
    if (vary_z_only):
      xcp = [xc[0], xc[1], xc[2] + random.gauss(0, sigma)]
    else:
      xcp = [e + random.gauss(0, sigma) for e in xc]
    xp = ucell.fractionalize(xcp)
    xps = SSx.ApplySpecialOp(xp)
    SSxps = get_site_symmetry(ucell, sgops, xps)
    if (str(SSxps.SpecialOp()) == str(SSx.SpecialOp())):
      have_position = 1
      break
  assert have_position, "Cannot find position matching all constraints."
  return xps

def random_rotate_ellipsoid(Ucart):
  C = rotation_parameters.amore_alpha_beta_gamma_as_matrix(
    [random.uniform(0,360) for i in xrange(3)]).elems
  return adptbx.CondensedTensorTransformation(C, Ucart)

class random_structure:

  def __init__(self, SgInfo, Elements,
               volume_per_atom=1000.,
               min_distance=3.0,
               general_positions_only=0,
               anisotropic_displacement_parameters=0):
    self.SgInfo = SgInfo
    self.UnitCell = get_compatible_unit_cell(
      self.SgInfo,
      len(Elements) * volume_per_atom * self.SgInfo.SgOps().OrderZ())
    self.Sites = shared.XrayScatterer()
    Positions = generate_positions(
      self.UnitCell, self.SgInfo.SgOps(), len(Elements),
      min_mate_distance = min_distance,
      min_hetero_distance = min_distance,
      general_positions_only = general_positions_only)
    n = 0
    for Elem, Pos in zip(Elements, Positions):
      n += 1
      SF = CAASF_WK1995(Elem)
      U = 0.01 + abs(random.gauss(0, 0.05))
      if (anisotropic_displacement_parameters):
        Uiso = U
        run_away_counter = 0
        while 1:
          run_away_counter += 1
          assert run_away_counter < 100
          U = [(Uiso + abs(random.gauss(Uiso, 0.05))) / 2
            for i in xrange(3)] + [0.,0.,0.]
          U = random_rotate_ellipsoid(U)
          U = adptbx.Ucart_as_Ustar(self.UnitCell, U)
          Site = sftbx.XrayScatterer(Elem + str(n), SF, 0j, Pos, 1., U)
          Site.ApplySymmetry(
            self.UnitCell, self.SgInfo.SgOps(), min_distance, 0, 0)
          U = adptbx.Ustar_as_Ucart(self.UnitCell, Site.Uaniso())
          try:
            Ev = adptbx.Eigenvalues(U)
          except RuntimeError:
            continue
          else:
            if (min(Ev) > 0.001):
              break
        U = Site.Uaniso()
      Site = sftbx.XrayScatterer(Elem + str(n), SF, 0j, Pos, 1., U)
      Site.ApplySymmetry(
        self.UnitCell, self.SgInfo.SgOps(), min_distance, 0, 1)
      self.Sites.append(Site)

class shake_structure:

  def __init__(self, structure, sigma = 0.0001, vary_z_only = 0):
    self.SgInfo = structure.SgInfo
    self.UnitCell = structure.UnitCell
    self.Sites = shared.XrayScatterer()
    for Site in structure.Sites:
      Site.set_Coordinates(perturb_coordinates(
        self.UnitCell, self.SgInfo.SgOps(),
        Site.Coordinates(), sigma, vary_z_only))
      self.Sites.append(Site)

class modify_sites:

  def __init__(self, structure, value_type, sigma = 0.0001):
    self.SgInfo = structure.SgInfo
    self.UnitCell = structure.UnitCell
    self.Sites = shared.XrayScatterer()
    for Site in structure.Sites:
      if (value_type == "Occ"):
        Site.set_Occ(max(0.1,
          Site.Occ() - abs(random.gauss(0, sigma))), self.SgInfo.SgOps())
      elif (value_type == "Uiso"):
        Site.set_Uiso(max(0.1,
          Site.Uiso() - abs(random.gauss(0, sigma))))
      else:
        raise RuntimeError, \
          "value_type " + str(value_type) + " not recognized."
      self.Sites.append(Site)

def print_sites(SymSites):
  print "Label  M  Coordinates            Occ   Uiso or Uaniso"
  print "     fp     fdp"
  for Site in SymSites:
    print "%-4s" % (Site.Label(),),
    print "%3d" % (Site.M(),),
    print "%7.4f %7.4f %7.4f" % Site.Coordinates(),
    print "%4.2f" % (Site.Occ(),),
    if (not Site.isAnisotropic()):
      print "%6.4f" % (Site.Uiso(),),
    else:
      print ("%6.3f " * 5 + "%6.3f") % adptbx.Ustar_as_Ucart(
        SymSites.UnitCell, Site.Uaniso()),
    print
    if (abs(Site.fpfdp()) != 0):
      print "     %6.4f %6.4f" % (Site.fpfdp().real, Site.fpfdp().imag)

def write_kriber(file_name, structure, key="z"):
  f = open(file_name, "a")
  print >> f, "*" + key
  print >> f
  print >> f
  print >> f, structure.SgInfo.BuildLookupSymbol()
  print >> f, structure.UnitCell
  for site in structure.Sites:
    print >> f, site.Label(), "%.6g %.6g %.6g" % site.Coordinates()
  print >> f, "-" * 40
  f.close()

def write_pdb(file_name, structure):
  f = open(file_name, "w")
  i = 0
  for site in structure.Sites:
    i += 1
    print >> f, "ATOM  %5d %-4s %-3s  %4d%3s" % (
      i, site.Label().upper(), site.Label().upper(), i, ""),
    print >> f, "%8.3f%8.3f%8.3f%6.2f%6.2f" % (
      structure.UnitCell.orthogonalize(site.Coordinates()) +
      (site.Occ(), adptbx.U_as_B(site.Uiso())))
  print >> f, "CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f %-11s" % (
    structure.UnitCell.getParameters()
    + (structure.SgInfo.BuildLookupSymbol(),))
  print >> f, "END"
  f.close()

def round_scaled(x, scale):
  y = round(x * scale)
  if (-0.5 < y <= 0): return 0.
  return y / scale

def print_structure_factors(F, precision_ampl=3, precision_phase=0):
  for i in xrange(len(F.H)):
    a, p = xutils.f_as_ampl_phase(F.F[i])
    a = round_scaled(a, 10**precision_ampl)
    p = round_scaled(p, 10**precision_phase) % 360
    if (p <= round_scaled(-180., 10**precision_phase)): p += 360
    print F.H[i], (
      "%%.%dg %%.%df" % (precision_ampl, precision_phase)) % (a, p)
