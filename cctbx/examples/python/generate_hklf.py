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
#   2001-04-20 Using new sgtbx::SymEquivCoordinates (Ralf W. Grosse-Kunstleve)
#   Created April 2001 (Ralf W. Grosse-Kunstleve)

import sys, exceptions, string, fileinput, math

import uctbx
import sgtbx
from eltbx.caasf_wk1995 import CAASF_WK1995

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
        self.Sf = CAASF_WK1995(flds[0], 0)
      coordinates = flds[offs : offs + 3]
      for i in xrange(3):
        coordinates[i] = string.atof(coordinates[i])
      self.Coordinates = coordinates
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

  def __str__(self):
    return "%s %s (%d) (%.6g %.6g %.6g) %.6g %.6g" % (
      (self.Label, self.Sf.Label(), self.Multiplicity)
      + tuple(self.Coordinates) + (self.Occ, self.Uiso))

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

  def read(self, files, SpaceGroupSymbolConvention = ""):
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
          self.SpaceGroupSymbols = sgtbx.SpaceGroupSymbols(sym,
            SpaceGroupSymbolConvention)

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

def ComputeStructureFactors(Sites, IndexSet):
  EightPiSquared = 8. * math.pi * math.pi
  SgOps = sgtbx.SgOps(IndexSet.SpaceGroupSymbols.Hall())
  FcalcDict = {}
  for H in IndexSet.IndexDict.keys():
    FcalcDict[H] = 0j
  for Site in Sites:
    SymEquivCoordinates = sgtbx.SymEquivCoordinates(
      IndexSet.UnitCell, SgOps, Site.Coordinates)
    Site.Multiplicity = SymEquivCoordinates.M()
    for H in IndexSet.IndexDict.keys():
      Q = IndexSet.IndexDict[H]
      stol2 = Q / 4.
      f0 = Site.Sf(stol2)
      B = EightPiSquared * Site.Uiso
      f = f0 * math.exp(-B * stol2) * Site.Occ
      FcalcDict[H] = FcalcDict[H] + f * SymEquivCoordinates.StructureFactor(H)
  return FcalcDict

def polar(c):
  from math import hypot, atan2, pi
  return hypot(c.real, c.imag), 180.0 * atan2(c.imag, c.real) / pi

if (__name__ == "__main__"):

  if (len(sys.argv) < 2):
    print __doc__
    sys.exit(1)

  Resolution_d_min = string.atof(sys.argv[1])

  Structure = StructureInfo()
  Structure.read(sys.argv[2:])

  for T in Structure.Titles: print "Title:", T
  print "Unit cell:", Structure.UnitCell
  print "Space group symbol:", \
        Structure.SpaceGroupSymbols.ExtendedHermann_Mauguin()
  print "Resolution d_min:", Resolution_d_min
  print "Number of sites:", len(Structure.Sites)
  print

  IndexSet = MillerIndexSet(Structure.SpaceGroupSymbols, Structure.UnitCell)
  IndexSet.BuildIndices(Resolution_d_min)

  FcalcDict = ComputeStructureFactors(Structure.Sites, IndexSet)

  print "Label ScatterLabel Multiplicity Coordinates Occupancy Uiso"
  for S in Structure.Sites: print S
  print

  print "Number of Miller indicies:", len(FcalcDict)
  print "H K L Fcalc"
  for H in FcalcDict.keys(): print H, "(%.6g, %.3f)" % polar(FcalcDict[H])
  print
