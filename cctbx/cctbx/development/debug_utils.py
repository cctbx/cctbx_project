import math
import random
from cctbx_boost.arraytbx import shared
from cctbx_boost import uctbx
from cctbx_boost import sgtbx
from cctbx_boost.eltbx.caasf_wk1995 import CAASF_WK1995
from cctbx_boost import adptbx
from cctbx_boost import sftbx
from cctbx.misc import python_utils
from cctbx import xutils
from cctbx.macro_mol import rotation_parameters

def set_random_seed(seed):
  random.seed(seed)

class command_line_options:

  def __init__(self, argv, keywords):
    self.keywords = keywords
    self.n = 0
    for keyword in keywords:
      setattr(self, keyword, 0)
    self.regular_args = []
    for arg in argv:
      if (not arg.startswith("--")):
        self.regular_args.append(arg)
        continue
      if (not arg[2:] in keywords):
        raise AssertionError, "Unknown option: " + arg
      setattr(self, arg[2:], 1)
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
  return xutils.crystal_symmetry(uc, SgInfo)

def generate_positions(N,
                       xsym, special_position_snap_parameters,
                       min_hetero_distance = 1.5,
                       general_positions_only = 0,
                       grid = None,
                       N_on_grid = 0,
                       existing_positions = []):
  existing_positions = list(existing_positions)
  max_back_track = 100
  max_placement_trials = 100
  for back_track in xrange(max_back_track):
    new_positions = []
    for i in xrange(N):
      have_position = 0
      for t in xrange(max_placement_trials):
        if (not grid or i < N - N_on_grid):
          X = (random.random(), random.random(), random.random())
        else:
          X = [random.randrange(g) / float(g) for g in grid]
        SS = sgtbx.SiteSymmetry(special_position_snap_parameters, X)
        if (general_positions_only and SS.M() != xsym.SgOps.OrderZ()):
          continue
        SS.expand()
        sec = sgtbx.SymEquivCoordinates(SS)
        min_d2 = xsym.UnitCell.getLongestVector2()
        for p in new_positions + existing_positions:
          min_d2 = min(min_d2, sec.getShortestDistance2(xsym.UnitCell, p))
        if (min_d2 >= min_hetero_distance**2):
          have_position = 1
          break
      if (not have_position): break
      new_positions.append(SS.SnapPosition())
    if (len(new_positions) == N): break
  assert len(new_positions) == N, \
    "Cannot find position matching all constraints."
  return new_positions

def random_adp():
  return random.random()

def random_rotate_ellipsoid(Ucart):
  C = rotation_parameters.amore_alpha_beta_gamma_as_matrix(
    [random.uniform(0,360) for i in xrange(3)]).elems
  return adptbx.CondensedTensorTransformation(C, Ucart)

class random_structure(xutils.symmetrized_sites):

  def __init__(self, SgInfo, elements,
               volume_per_atom = 1000.,
               min_distance = 3.0,
               general_positions_only = 0,
               random_f_prime_d_min = 0,
               random_f_double_prime = 0,
               no_random_u = 0,
               uiso = 0,
               anisotropic_displacement_parameters = 0,
               defer_build = 0):
    python_utils.adopt_init_args(self, locals())
    xutils.symmetrized_sites.__init__(self,
      get_compatible_unit_cell(
        SgInfo,
        len(self.elements) * self.volume_per_atom * SgInfo.SgOps().OrderZ()),
      MinMateDistance = min_distance,
      Ustar_tolerance = 0,
      TestPositiveDefiniteness = 1)
    if (not defer_build):
      self.build_sites()

  def build_sites(self, grid = None, N_on_grid = 0):
    assert self.Sites.size() == 0
    positions = generate_positions(len(self.elements),
      self, self.SnapParameters,
      min_hetero_distance = self.MinMateDistance,
      general_positions_only = self.general_positions_only,
      grid = grid,
      N_on_grid = N_on_grid)
    fp = 0
    fdp = 0
    n = 0
    for Elem, Pos in zip(self.elements, positions):
      n += 1
      SF = CAASF_WK1995(Elem)
      if (self.random_f_prime_d_min):
        f0 = SF.stol(1./(2*self.random_f_prime_d_min))
        assert f0 > 0
        fp = -min(f0*0.9, random.gauss(f0, f0/3))
      if (self.random_f_double_prime):
        fdp = random.gauss(10, 3)
      fpfdp = complex(fp, fdp)
      U = self.uiso
      if (not U and not self.no_random_u):
        U = random_adp()
      if (self.anisotropic_displacement_parameters):
        run_away_counter = 0
        while 1:
          run_away_counter += 1
          assert run_away_counter < 100
          U = [math.fmod(i * random_adp(), 1) for i in xrange(3)] + [0.,0.,0.]
          U = random_rotate_ellipsoid(U)
          U = adptbx.Ucart_as_Ustar(self.UnitCell, U)
          Site = sftbx.XrayScatterer(Elem + str(n), SF, fpfdp, Pos, 1., U)
          Site.ApplySymmetry(
            self.UnitCell, self.SgOps, self.MinMateDistance, 0, 0)
          U = adptbx.Ustar_as_Ucart(self.UnitCell, Site.Uaniso())
          try:
            Ev = adptbx.Eigenvalues(U)
          except RuntimeError:
            continue
          else:
            if (min(Ev) > 0.001):
              break
        U = Site.Uaniso()
      Site = sftbx.XrayScatterer(Elem + str(n), SF, fpfdp, Pos, 1., U)
      self.add_site(Site)

def shake_position(xsym, special_position_snap_parameters,
                   x, sigma, max_diff = 0, vary_z_only = 0):
  SSx = sgtbx.SiteSymmetry(special_position_snap_parameters, x)
  xc = xsym.UnitCell.orthogonalize(x)
  max_placement_trials = 100
  have_position = 0
  for t in xrange(max_placement_trials):
    if (vary_z_only):
      xcp = [xc[0], xc[1], xc[2] + random.gauss(0, sigma)]
    else:
      xcp = [e + random.gauss(0, sigma) for e in xc]
    xp = xsym.UnitCell.fractionalize(xcp)
    xps = SSx.ApplySpecialOp(xp)
    SSxps = sgtbx.SiteSymmetry(special_position_snap_parameters, xps)
    if (str(SSxps.SpecialOp()) != str(SSx.SpecialOp())): continue
    if (max_diff):
      diff = xsym.UnitCell.Distance2(x, xps)
      if (diff > max_diff): continue
    have_position = 1
    break
  assert have_position, "Cannot find position matching all constraints."
  return xps

def shake_structure(xtal, sigma = 0.0001, vary_z_only = 0):
  xtal_shake = xtal.copy_attributes()
  for site in xtal:
    site.set_Coordinates(shake_position(
      xtal, xtal.SnapParameters,
      site.Coordinates(), sigma, vary_z_only))
    xtal_shake.add_site(site)
  return xtal_shake

def random_modify_atomic_parmeters(xtal, parameter_type, sigma = 0.0001):
  xtal_mod = xtal.copy_attributes()
  for site in xtal:
    if (parameter_type == "Occ"):
      site.set_Occ(max(0.1,
        site.Occ() - abs(random.gauss(0, sigma))), xtal.SgOps)
    elif (parameter_type == "Uiso"):
      site.set_Uiso(max(0.1,
        site.Uiso() - abs(random.gauss(0, sigma))))
    else:
      raise RuntimeError, \
        "parameter_type " + str(parameter_type) + " not recognized."
    xtal_mod.add_site(site)
  return xtal_mod

def print_sites(xtal):
  print "Label  M  Coordinates            Occ   Uiso or Uaniso"
  print "     fp     fdp"
  for Site in xtal:
    print "%-4s" % (Site.Label(),),
    print "%3d" % (Site.M(),),
    print "%7.4f %7.4f %7.4f" % Site.Coordinates(),
    print "%4.2f" % (Site.Occ(),),
    if (not Site.isAnisotropic()):
      print "%6.4f" % (Site.Uiso(),),
    else:
      print ("%6.3f " * 5 + "%6.3f") % adptbx.Ustar_as_Ucart(
        xtal.UnitCell, Site.Uaniso()),
    print
    if (abs(Site.fpfdp()) != 0):
      print "     %6.4f %6.4f" % (Site.fpfdp().real, Site.fpfdp().imag)

def write_kriber(file_name, xtal, key="z"):
  f = open(file_name, "a")
  print >> f, "*" + key
  print >> f
  print >> f
  print >> f, xtal.SgInfo.BuildLookupSymbol()
  print >> f, xtal.UnitCell
  for site in xtal.Sites:
    print >> f, site.Label(), "%.6g %.6g %.6g" % site.Coordinates()
  print >> f, "-" * 40
  f.close()

def write_pdb(file_name, xtal):
  f = open(file_name, "w")
  i = 0
  for site in xtal.Sites:
    i += 1
    print >> f, "ATOM  %5d %-4s %-3s  %4d%3s" % (
      i, site.Label().upper(), site.Label().upper(), i, ""),
    print >> f, "%8.3f%8.3f%8.3f%6.2f%6.2f" % (
      xtal.UnitCell.orthogonalize(site.Coordinates()) +
      (site.Occ(), adptbx.U_as_B(site.Uiso())))
  print >> f, "CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f %-11s" % (
    xtal.UnitCell.getParameters()
    + (xtal.SgInfo.BuildLookupSymbol(),))
  print >> f, "END"
  f.close()

def round_scaled(x, scale):
  y = round(x * scale)
  if (-0.5 < y <= 0): return 0.
  return y / scale

def round_scaled_list(l, scale):
  return [round_scaled(x, scale) for x in l]

def format_round_scaled_list(l, precision=4, fmt="f"):
  s = ""
  if (len(l) == 0): return s
  for x in round_scaled_list(l, 10**precision):
    s += (" %%.%d%s" % (precision, fmt)) % (x,)
  return s[1:]

def format_structure_factor(f, precision_ampl=3, precision_phase=0):
  a, p = xutils.f_as_ampl_phase(f)
  a = round_scaled(a, 10**precision_ampl)
  p = round_scaled(p, 10**precision_phase) % 360
  if (p <= round_scaled(-180., 10**precision_phase)): p += 360
  return ("%%.%dg %%.%df" % (precision_ampl, precision_phase)) % (a, p)

def print_structure_factors(F, precision_ampl=3, precision_phase=0):
  if (type(F) == type(0j)):
    for i in xrange(len(F.H)):
      print F.H[i], format_structure_factor(
        F.F[i], precision_ampl, precision_phase)
  else:
    for i in xrange(len(F.H)):
      print F.H[i], round_scaled(F.F[i], 10**precision_ampl)

def show_regression(x, y, label, min_correlation = 0):
  xy_regr = shared.linear_regression(x, y)
  assert xy_regr.is_well_defined()
  print label, "cc: %.4f m: %.3f" % (xy_regr.cc(), xy_regr.m())
  assert min_correlation == 0 or xy_regr.cc() >= min_correlation

def phase_error(p1, p2):
  d_as_r = math.pi / 180
  return math.acos(math.cos((p1 - p2) * d_as_r)) / d_as_r

class structure_factor_comparison:

  def __init__(self, label, min_corr_ampl=0, max_mean_w_phase_error=0,
               verbose=0):
    python_utils.adopt_init_args(self, locals())
    self.amp1 = shared.double()
    self.amp2 = shared.double()
    self.sum_amp1_minus_amp2_sq = 0
    self.sum_amp1_sq = 0
    self.sum_w_phase_error = 0
    self.sum_w = 0

  def add(self, h, f1, f2):
    a1, p1 = xutils.f_as_ampl_phase(f1)
    a2, p2 = xutils.f_as_ampl_phase(f2)
    if (self.verbose):
      print h
      print " ", a1, p1
      print " ", a2, p2
      print " " * 20, phase_error(p1, p2)
    self.amp1.append(a1)
    self.amp2.append(a2)
    self.sum_amp1_minus_amp2_sq += (a1 - a2)**2
    self.sum_amp1_sq += a1**2
    self.sum_w_phase_error += (a1 + a2) * phase_error(p1, p2)
    self.sum_w += (a1 + a2)

  def report(self):
    if (self.sum_amp1_sq):
      r = self.sum_amp1_minus_amp2_sq / self.sum_amp1_sq
      print self.label, "R-factor: %.3f" % (r,)
    if (self.sum_w):
      self.mean_w_phase_error = self.sum_w_phase_error / self.sum_w
    show_regression(
      self.amp1, self.amp2, self.label + " ampl", self.min_corr_ampl)
    print self.label + (" mean weighted phase error: %.2f" % (
      self.mean_w_phase_error,))
    if (self.max_mean_w_phase_error):
      assert self.mean_w_phase_error <= self.max_mean_w_phase_error

def show_structure_factor_correlation(label, h1, joined_sets, f1, f2,
                                      min_corr_ampl=0,
                                      max_mean_w_phase_error=0,
                                      verbose=0):
  sf_cmp = structure_factor_comparison(
    label, min_corr_ampl, max_mean_w_phase_error, verbose)
  if (joined_sets == 0):
    assert f1.size() == f2.size()
    for i in xrange(f1.size()):
      sf_cmp.add(h1[i], f1[i], f2[i])
  else:
    for i,j in joined_sets.pairs():
      sf_cmp.add(h1[i], f1[i], f2[j])
  sf_cmp.report()
