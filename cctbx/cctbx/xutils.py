from cctbx_boost.arraytbx import shared
from cctbx_boost import sgtbx
from cctbx_boost import sftbx

class crystal_symmetry:

  def __init__(self, UnitCell, SgInfo):
    self.UnitCell = UnitCell
    self.SgInfo = SgInfo
    self.SgOps = self.SgInfo.SgOps() # FUTURE: avoid the copy?

  def cell_equivalent_p1(self):
    return crystal_symmetry(self.UnitCell, sgtbx.SpaceGroupInfo())

class miller_index_set(crystal_symmetry):

  def __init__(self, xsym, H):
    crystal_symmetry.__init__(self, xsym.UnitCell, xsym.SgInfo)
    self.H = H

  def expand_indices_to_p1(self, friedel_flag = 1):
    set_p1_H = shared.Miller_Index()
    sgtbx.expand_to_p1(self.SgOps, self.H, set_p1_H, friedel_flag)
    return miller_index_set(self.cell_equivalent_p1(), set_p1_H)

class reciprocal_space_array(miller_index_set):

  def __init__(self, miller_indices, F):
    miller_index_set.__init__(self, miller_indices, miller_indices.H)
    self.F = F

class symmetrized_sites(crystal_symmetry):

  def __init__(self, xsym, sites = None,
               MinMateDistance = 0.5,
               Ustar_tolerance = 0.1,
               TestPositiveDefiniteness = 1):
    crystal_symmetry.__init__(self, xsym.UnitCell, xsym.SgInfo)
    self.MinMateDistance = MinMateDistance
    self.Ustar_tolerance = Ustar_tolerance
    self.TestPositiveDefiniteness = TestPositiveDefiniteness
    self.Sites = shared.XrayScatterer()
    self.SpecialPositionOps = shared.RTMx()
    self.SnapParameters = sgtbx.SpecialPositionSnapParameters(
      self.UnitCell, self.SgOps, 1, self.MinMateDistance)
    if (sites):
      self.add_sites(sites)

  def copy_attributes(self):
    return symmetrized_sites(self, sites = None,
      MinMateDistance = self.MinMateDistance,
      Ustar_tolerance = self.Ustar_tolerance,
      TestPositiveDefiniteness = self.TestPositiveDefiniteness)

  def n_special_positions(self):
    n = 0
    for special_op in self.SpecialPositionOps:
      if (not special_op.isUnit()):
        n += 1
    return n

  def get_site_symmetry(self, X):
    return sgtbx.SiteSymmetry(self.SnapParameters, X)

  def add_site(self, site):
    site_symmetry = site.ApplySymmetry(
      self.UnitCell, self.SgOps,
      self.MinMateDistance,
      self.Ustar_tolerance,
      self.TestPositiveDefiniteness)
    special_op = site_symmetry.SpecialOp()
    self.Sites.append(site)
    self.SpecialPositionOps.append(special_op)

  def add_sites(self, sites):
    for site in sites:
      self.add_site(site)

  def __len__(self): return self.size()
  def size(self): return self.Sites.size()

  def __getitem__(self, key):
    return self.Sites[key]

def build_miller_indices(xsym, d_min): # XXX friedel_flag
  return miller_index_set(
    xsym, sftbx.BuildMillerIndices(xsym.UnitCell, xsym.SgInfo, d_min))

def calculate_structure_factors(miller_indices, xtal, abs_F = 0):
  F = sftbx.StructureFactorArray(
    miller_indices.UnitCell, miller_indices.SgOps, miller_indices.H,
    xtal.Sites)
  if (abs_F): F = shared.abs(F)
  result = reciprocal_space_array(miller_indices, F)
  result.symmetrized_sites = xtal
  return result

def f_as_ampl_phase(f, deg=1):
  from math import hypot, atan2, pi
  if (abs(f) == 0): return (0., 0.)
  x, y = f.real, f.imag
  if deg:
    return hypot(x, y), atan2(y, x) * 180. / pi
  else:
    return hypot(x, y), atan2(y, x)

def ampl_phase_as_f(ampl_phase, deg=1):
  from math import cos, sin, pi
  a, p = ampl_phase
  if deg:
    p *= pi / 180.0
  return complex(a * cos(p), a * sin(p))
