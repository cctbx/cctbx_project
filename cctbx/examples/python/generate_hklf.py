"""
# This example script reads in the parameters for a structure
# (unit cell, space group, fractional coordinates, etc.) from a file,
# builds a list of Miller indices and computes structure factors
# for an asymmetric unit of reflections.
# The resolution limit is specified on the command line. E.g.:
#
#       python generate_hklf.py 2. vni.structure
#
# This is, the high resolution limit is 2 Angstrom, and the
# structure parameters are read from the file vni.structure.
"""

# Revision history:
#   Created April 2001 (Ralf W. Grosse-Kunstleve)

import sys, exceptions, string, fileinput, math

import uctbx
import sgtbx
from eltbx.caasf_wk1995 import CAASF_WK1995

class Vector:

  def __init__(self, Coordinates = (0., 0., 0.)):
    self.Coordinates = list(Coordinates)
  def __str__(self):
    return "%.6g %.6g %.6g" % tuple(self.Coordinates)
  def __add__(self, other):
    Sum = list(self.Coordinates)
    for i in xrange(3): Sum[i] = Sum[i] + other.Coordinates[i]
    return Vector(Sum)
  def __sub__(self, other):
    Diff = list(self.Coordinates)
    for i in xrange(3): Diff[i] = Diff[i] - other.Coordinates[i]
    return Vector(Diff)
  def __rmul__(self, other):
    if (hasattr(other, "as_tuple")):
      S = other.as_tuple()
      R = Matrix(S[0]) # rotation part
      T = Vector(S[1]) # translation part
      return R * self + T
    assert(len(other) == len(self.Coordinates))
    result = 0.
    for i in xrange(3): result = result + other[i] * self.Coordinates[i]
    return result
  def ModShort(self):
    Short = list(self.Coordinates)
    for i in xrange(3):
      Short[i] = math.fmod(Short[i], 1.)
      if (Short[i] < -.5):
        Short[i] = Short[i] + 1.
      elif (Short[i] > .5):
        Short[i] = Short[i] - 1.
    return Vector(Short)
  def Length2(self):
    L2 = 0.
    for i in xrange(3):
      L2 = L2 + float(self.Coordinates[i]) ** 2
    return L2

class Matrix:

  def __init__(self, M = (1, 0, 0, 0, 1, 0, 0, 0, 1)):
    self.M = list(M)
  def __str__(self):
    return "((%.6g,%.6g,%.6g), (%.6g,%.6g,%.6g), (%.6g,%.6g,%.6g))" % tuple(self.M)
  def __mul__(self, other):
    assert (type(other) == type(Vector()))
    result = Vector()
    for i in xrange(3):
      s = 0.
      for j in xrange(3):
        s = s + self.M[i * 3 + j] * other.Coordinates[j]
      result.Coordinates[i] = s
    return result

class SiteInfo:

  def __init__(self, flds, default_Uiso = 0.035):
    # Label [ScatFact] x y z [Occ [Uiso]]
    try:
      self.Label = flds[0]
      try:
        string.atof(flds[1])
      except:
        offs = 2
        self.Sf = CAASF_WK1995(flds[1], 1)
      else:
        offs = 1
        try: self.Sf = CAASF_WK1995(flds[0], 0)
        except: # special case for T (my zeolite past shows through)
          if (    flds[0][0] == "T"
              and (len(flds[0]) == 1 or flds[0][1] in string.digits)):
            self.Sf = CAASF_WK1995("Si", 1)
          else:
            raise
      coordinates = flds[offs : offs + 3]
      for i in xrange(3):
        coordinates[i] = string.atof(coordinates[i])
      self.Coordinates = Vector(coordinates)
      self.Occ = 1.
      self.Uiso = default_Uiso
      if (len(flds) >= offs + 4):
        self.Occ = flds[offs + 3]
        if (len(flds) == offs + 5):
          self.Uiso = flds[offs + 4]
        else:
          raise FormatError
      self.Multiplicity = 0
    except:
      raise FormatError

  def __repr__(self):
    return "%s %s (%d) (%s) %.6g %.6g" % (
      self.Label, self.Sf.Label(), self.Multiplicity,
      self.Coordinates, self.Occ, self.Uiso)

class FormatError(exceptions.Exception):
  pass

def strip_comment(line):
  i = string.find(line, "#")
  if (i >= 0): line = line[:i]
  return line

def strip_keyword(line, keyword):
  return string.strip(string.lstrip(line)[len(keyword) + 1:])

class StructureInfo:

  def __init__(self, UnitCell = uctbx.UnitCell(()),
                     SpaceGroupSymbols = sgtbx.SpaceGroupSymbols(1),
                     Resolution_d_min = 1.):
    self.Titles = []
    self.UnitCell = UnitCell
    self.SpaceGroupSymbols = SpaceGroupSymbols
    self.Resolution_d_min = Resolution_d_min
    self.Sites = []

  def read(self, files):
    iline = 0
    try:
      input = fileinput.input(files)
      for line in input:
        iline = iline + 1
        line = line[:-1]
        flds = string.split(string.strip(strip_comment(line)))
        if (len(flds) == 0): continue
        keyword = string.capwords(flds[0])

        if (keyword == "Title"):
          self.Titles.append(strip_keyword(line, keyword))

        elif (keyword == "Unitcell"):
          uc_params = flds[1:]
          if (not (0 < len(uc_params) <= 6)):
            raise FormatError
          for i in xrange(len(uc_params)):
            uc_params[i] = string.atof(uc_params[i])
          self.UnitCell = uctbx.UnitCell(uc_params)

        elif (keyword == "Spacegroup"):
          sym = string.join(flds[1:])
          self.SpaceGroupSymbols = sgtbx.SpaceGroupSymbols(sym)

        elif (keyword == "Resolution"):
          if (   string.lower(flds[1]) != "d_min"
              or len(flds) != 3): raise FormatError
          self.Resolution_d_min = string.atof(flds[2])

        elif (keyword == "End"):
          input.nextfile()

        else:
          if (keyword == "Site"):
            site = SiteInfo(flds[1:])
          else:
            site = SiteInfo(flds)
          self.Sites.append(site)

    except FormatError:
      print "Input:", line
      raise

class MillerIndexSet:

  def __init__(self, SpaceGroupSymbols, UnitCell):
    self.SpaceGroupSymbols = SpaceGroupSymbols
    self.UnitCell = UnitCell
    self.IndexDict = {}

  def BuildIndices(self, Resolution_d_min, FriedelSym = 1):
    SgOps = sgtbx.SgOps(self.SpaceGroupSymbols.Hall())
    SgOps.CheckUnitCell(self.UnitCell)
    CutP = SgOps.getCutParameters(FriedelSym)
    Hmax = self.UnitCell.MaxMillerIndices(Resolution_d_min)
    Qlow = 0.
    Qhigh = 1. / (Resolution_d_min * Resolution_d_min)
    Hmin = [0] * 3
    for i in xrange(3): Hmin[i] = CutP[i] * Hmax[i]

    # short-cuts for (slightly) better performance
    isSysAbsent = SgOps.isSysAbsent
    UnitCell_Q = self.UnitCell.Q

    # loop over all possible indices
    H = [0] * 3
    for H[0] in xrange(Hmin[0], Hmax[0] + 1):
      for H[1] in xrange(Hmin[1], Hmax[1] + 1):
        for H[2] in xrange(Hmin[2], Hmax[2] + 1):
          Q = UnitCell_Q(H)
          if (Q != 0 and Qlow <= Q <= Qhigh): # resolution filter
            if (not isSysAbsent(H)):
              Master = SgOps.getMasterIndex(H, CutP, 1)
              key = tuple(H)
              if (key == Master.H()):
                self.IndexDict[key] = Q

def getMinDelta2(UnitCell, SymEquivCoordinates, SX):
  MinDelta2 = None
  for C in SymEquivCoordinates:
    Delta = (SX - C).ModShort()
    CartDelta = Vector(UnitCell.orthogonalize(Delta.Coordinates))
    CartDelta2 = CartDelta.Length2()
    if (MinDelta2 == None or MinDelta2 > CartDelta2):
      MinDelta2 = CartDelta2
  return MinDelta2

def getEquivCoordinates(UnitCell, SgOps, X, MinimumDistance = 0.05):
  SymEquivCoordinates = []
  for i in xrange(SgOps.OrderZ()):
    SX = SgOps(i) * X
    MinDelta2 = getMinDelta2(UnitCell, SymEquivCoordinates, SX)
    if (MinDelta2 == None or MinDelta2 >= MinimumDistance ** 2):
      SymEquivCoordinates.append(SX)
  assert(SgOps.OrderZ() % len(SymEquivCoordinates) == 0)
  return SymEquivCoordinates

def ComputeStructureFactors(Sites, IndexSet):
  EightPiSquared = 8. * math.pi * math.pi
  TwoPi = 2. * math.pi
  SgOps = sgtbx.SgOps(IndexSet.SpaceGroupSymbols.Hall())
  FcalcDict = {}
  for H in IndexSet.IndexDict.keys():
    FcalcDict[H] = 0j
  for Site in Sites:
    SymEquivCoordinates = getEquivCoordinates(
      IndexSet.UnitCell, SgOps, Site.Coordinates)
    Site.Multiplicity = len(SymEquivCoordinates)
    for H in IndexSet.IndexDict.keys():
      Q = IndexSet.IndexDict[H]
      stol2 = Q / 4.
      f0 = Site.Sf(stol2)
      B = EightPiSquared * Site.Uiso
      fs = f0 * math.exp(-B * stol2) * Site.Occ
      F = FcalcDict[H]
      for X in SymEquivCoordinates:
        phase = TwoPi * (H * X)
        F = F + complex(fs * math.cos(phase), fs * math.sin(phase))
      FcalcDict[H] = F
  return FcalcDict

# http://www.python.org/topics/scicomp/recipes_in_python.html
"""
Convert from polar (r,w) to rectangular (x,y)
    x = r cos(w)
    y = r sin(w)
"""
def rect(r, w, deg=0):          # radian if deg=0; degree if deg=1
    from math import cos, sin, pi
    if deg:
        w = pi * w / 180.0
    return r * cos(w), r * sin(w)

"""
Convert from rectangular (x,y) to polar (r,w)
    r = sqrt(x^2 + y^2)
    w = arctan(y/x) = [-\pi,\pi] = [-180,180]
"""
def polar(x, y, deg=0):         # radian if deg=0; degree if deg=1
    from math import hypot, atan2, pi
    if deg:
        return hypot(x, y), 180.0 * atan2(y, x) / pi
    else:
        return hypot(x, y), atan2(y, x)

def cpolar(c, deg = 0):
  return polar(c.real, c.imag, deg)

if (__name__ == "__main__"):

  if (len(sys.argv) < 2):
    print __doc__
    sys.exit(1)

  Resolution_d_min = string.atof(sys.argv[1])

  Structure = StructureInfo()
  Structure.read(sys.argv[2:])

  for T in Structure.Titles: print "Title:", T
  print "Unit cell:", Structure.UnitCell
  print "Space group symbol:", Structure.SpaceGroupSymbols.Hermann_Mauguin()
  print "Resolution d_min:", Resolution_d_min
  print "Number of sites:", len(Structure.Sites)
  print

  IndexSet = MillerIndexSet(Structure.SpaceGroupSymbols, Structure.UnitCell)
  IndexSet.BuildIndices(Resolution_d_min)

  FcalcDict = ComputeStructureFactors(Structure.Sites, IndexSet)

  print "Label ScatterLabel Multiplicity Coordinates Occupancy Uiso"
  for S in Structure.Sites: print S
  print

  print "H K L Fcalc"
  for H in FcalcDict.keys(): print H, cpolar(FcalcDict[H], 1)
  print
