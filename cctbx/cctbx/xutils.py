from cctbx_boost.arraytbx import shared
from cctbx_boost import sftbx

class crystal_symmetry:

  def __init__(self, UnitCell, SgInfo):
    self.UnitCell = UnitCell
    self.SgInfo = SgInfo

  def SgOps(self):
    return self.SgInfo.SgOps()

class miller_index_set(crystal_symmetry):

  def __init__(self, xsym, d_min):
    self.__dict__.update(xsym.__dict__)
    self.d_min = d_min
    self.H = sftbx.BuildMillerIndices(self.UnitCell, self.SgInfo, d_min)

class symmetrized_sites(crystal_symmetry):

  def __init__(self, xsym, sites,
               MinMateDistance = 0.5,  
               Ustar_tolerance = 0.1,  
               TestPositiveDefiniteness = 1):
    self.__dict__.update(xsym.__dict__)
    self.Sites = sites
    self.MinMateDistance = MinMateDistance
    self.Ustar_tolerance = Ustar_tolerance
    self.TestPositiveDefiniteness = TestPositiveDefiniteness
    self.SymSites = shared.XrayScatterer()
    self.SpecialPositionOps = shared.RTMx()
    self.nSpecialPositions = 0
    for site in sites:
      site_symmetry = site.ApplySymmetry(
        xsym.UnitCell, xsym.SgOps(),
        MinMateDistance, Ustar_tolerance, TestPositiveDefiniteness)
      special_op = site_symmetry.SpecialOp()
      self.SymSites.append(site)
      self.SpecialPositionOps.append(special_op)
      if (not special_op.isUnit()):
        self.nSpecialPositions += 1

  def __getitem__(self, key):
    return self.SymSites[key]

class structure_factors(miller_index_set, symmetrized_sites):

  def __init__(self, MillerIndices, SymSites, abs_F = 0):
    self.__dict__.update(SymSites.__dict__)
    self.__dict__.update(MillerIndices.__dict__)
    self.F = sftbx.StructureFactorArray(
      self.UnitCell, self.SgOps(), self.H, self.SymSites)
    if (abs_F):
      abs_F = shared.double()
      for F in self.F:
        abs_F.append(abs(F))
      self.F = abs_F

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
