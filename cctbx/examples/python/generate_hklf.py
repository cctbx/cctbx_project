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
#   2001 Sep 13: SpaceGroupType -> SpaceGroupInfo (R.W. Grosse-Kunstleve)
#   2001 Jul 02: Merged from CVS branch sgtbx_special_pos (rwgk)
#   2001 Jun 21: use WyckoffTable (R.W. Grosse-Kunstleve)
#   2001 May 31: merged from CVS branch sgtbx_type (R.W. Grosse-Kunstleve)
#   2001-04-20 Using new sgtbx::SymEquivCoordinates (Ralf W. Grosse-Kunstleve)
#   Created April 2001 (Ralf W. Grosse-Kunstleve)

import sys, exceptions, string, fileinput, math
class FormatError(exceptions.Exception): pass

import uctbx
import sgtbx
from eltbx.caasf_wk1995 import CAASF_WK1995

class SiteInfo:

  def __init__(self, flds, default_Biso = 3.0):
    # Label [ScatFact] x y z [Occ [Biso]]
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
      self.Biso = default_Biso
      if (len(flds) >= offs + 4):
        self.Occ = string.atof(flds[offs + 3])
        if (len(flds) == offs + 5):
          self.Biso = string.atof(flds[offs + 4])
        elif (len(flds) > offs + 5):
          raise FormatError, flds
      self.WyckoffMapping = None
      self.SiteSymmetry = None
    except:
      raise FormatError, flds

  def __str__(self):
    return "%s %s (%d %s %s) (%.6g %.6g %.6g) %.6g %.6g" % (
      (self.Label, self.Sf.Label(),
       self.WyckoffMapping.WP().M(), self.WyckoffMapping.WP().Letter(),
       self.SiteSymmetry)
      + tuple(self.Coordinates) + (self.Occ, self.Biso))

def strip_comment(line):
  i = string.find(line, "#")
  if (i >= 0): line = line[:i]
  return line

def strip_keyword(line, keyword):
  return string.strip(string.lstrip(line)[len(keyword) + 1:])

class StructureInfo:

  def __init__(self, UnitCell = uctbx.UnitCell(()),
                     SgOps = sgtbx.SpaceGroup(),
                     Resolution_d_min = 1.):
    self.Titles = []
    self.UnitCell = UnitCell
    self.SgOps = SgOps
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
          SpaceGroupSymbols = sgtbx.SpaceGroupSymbols(sym)
          self.SgOps = sgtbx.SpaceGroup(SpaceGroupSymbols.Hall())
          self.SgInfo = self.SgOps.Info()
          self.WyckoffTable = sgtbx.WyckoffTable(self.SgInfo, 1)

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

    self.SgOps.CheckUnitCell(self.UnitCell)

    SnapParameters = \
      sgtbx.SpecialPositionSnapParameters(self.UnitCell, self.SgOps)
    for Site in self.Sites:
      SS = sgtbx.SiteSymmetry(SnapParameters, Site.Coordinates, 0)
      Site.WyckoffMapping = self.WyckoffTable.getWyckoffMapping(SS)
      Site.SiteSymmetry = SS.PointGroupType()

class MillerIndexSet:

  def __init__(self, SgOps, UnitCell):
    self.SgOps = SgOps
    self.UnitCell = UnitCell
    self.IndexDict = {}

  def BuildIndices(self, Resolution_d_min):
    SgOps = self.SgOps
    CutP = SgOps.getCutParameters()
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
                self.IndexDict[key] = Q / 4.

def ComputeStructureFactors(Sites, IndexSet):
  FcalcDict = {}
  for H in IndexSet.IndexDict.keys(): FcalcDict[H] = 0j
  for Site in Sites:
    SEC = sgtbx.SymEquivCoordinates(Site.WyckoffMapping,
                                    Site.Coordinates)
    for H in IndexSet.IndexDict.keys():
      stol2 = IndexSet.IndexDict[H]
      f0 = Site.Sf.stol2(stol2)
      f = f0 * math.exp(-Site.Biso * stol2) * Site.Occ
      FcalcDict[H] = FcalcDict[H] + f * SEC.StructureFactor(H)
  return FcalcDict

def polar(c):
  from math import hypot, atan2, pi
  if (c == 0j): return 0., 0.
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
  print "Space group symbol:", Structure.SgInfo.BuildLookupSymbol()
  print "Resolution d_min:", Resolution_d_min
  print "Number of sites:", len(Structure.Sites)
  print

  print "Label ScatterLabel WyckoffPosition SiteSymmetry",
  print "Coordinates Occupancy Biso"
  for S in Structure.Sites: print S
  print

  IndexSet = MillerIndexSet(Structure.SgOps, Structure.UnitCell)
  IndexSet.BuildIndices(Resolution_d_min)

  FcalcDict = ComputeStructureFactors(Structure.Sites, IndexSet)

  print "Number of Miller indices:", len(FcalcDict)
  print "H K L Fcalc"
  for H in FcalcDict.keys(): print H, "(%.6g, %.3f)" % polar(FcalcDict[H])
  print
