#! /usr/local/Python-2.1/bin/python

# Revision history:
#   2001 Sep 13: SpaceGroupType -> SpaceGroupInfo (R.W. Grosse-Kunstleve)
#   2001 Aug 09: Derived from generate_hklf.py (rwgk)

PATH_cctbx_lib_python = "/net/boa/srv/html/cci/cctbx/lib/python"

import sys
sys.stderr = sys.stdout

print "Content-type: text/html"
print

import traceback
import exceptions
class FormatError(exceptions.Exception): pass

import math, string, cgi

sys.path.insert(0, PATH_cctbx_lib_python)
import sgtbx
import uctbx
from eltbx.caasf_wk1995 import CAASF_WK1995

print "sgtbx version:", sgtbx.__version__
print "<br>"
print "uctbx version:", uctbx.__version__
print "<p>"
print "<pre>"
InTable = 0

class Empty: pass

def GetFormData():
  form = cgi.FieldStorage()
  inp = Empty()
  for key in (("ucparams", "1 1 1 90 90 90"),
              ("sgsymbol", "P1"),
              ("convention", ""),
              ("d_min", "1"),
              ("MinMateDistance", "0.5"),
              ("coor_type", None)):
    if (form.has_key(key[0])):
      inp.__dict__[key[0]] = string.strip(form[key[0]].value)
    else:
      inp.__dict__[key[0]] = key[1]
  inp.coordinates = []
  if (form.has_key("coordinates")):
    lines = string.split(form["coordinates"].value, "\015\012")
    for l in lines:
      s = string.strip(l)
      if (len(s) != 0): inp.coordinates.append(s)
  return inp

def ShowInputSymbol(sgsymbol, convention, label):
  if (sgsymbol != ""):
    print label, "space group symbol:", sgsymbol
    print "Convention:",
    if   (convention == "A1983"):
      print "International Tables for Crystallography, Volume A 1983"
    elif (convention == "I1952"):
      print "International Tables for Crystallography, Volume I 1952"
    elif (convention == "Hall"):
      print "Hall symbol"
    else:
      print "Default"

def HallSymbol_to_SgOps(HallSymbol):
  try:
    ps = sgtbx.parse_string(HallSymbol)
    SgOps = sgtbx.SpaceGroup(ps)
  except RuntimeError, e:
    print "-->" + ps.string() + "<--"
    print ("-" * (ps.where() + 3)) + "^"
    raise
  return SgOps

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

def BuildMillerIndices(UnitCell, SgInfo, Resolution_d_min):
  MIG = sgtbx.MillerIndexGenerator(UnitCell, SgInfo, Resolution_d_min)
  MillerIndices = []
  for H in MIG: MillerIndices.append(H)
  return MillerIndices

def ComputeStructureFactors(Sites, MillerIndices):
  FcalcDict = {}
  for H in MillerIndices: FcalcDict[H] = 0j
  for Site in Sites:
    SEC = sgtbx.SymEquivCoordinates(Site.WyckoffMapping,
                                    Site.Coordinates)
    for H in MillerIndices:
      stol2 = UnitCell.Q(H) / 4.
      f0 = Site.Sf.stol2(stol2)
      f = f0 * math.exp(-Site.Biso * stol2) * Site.Occ
      FcalcDict[H] += f * SEC.StructureFactor(H)
  return FcalcDict

def polar(c):
  from math import hypot, atan2, pi
  if (c == 0j): return 0., 0.
  return hypot(c.real, c.imag), 180.0 * atan2(c.imag, c.real) / pi

if (__name__ == "__main__"):

  inp = GetFormData()

  try:
    u = string.split(inp.ucparams)
    for i in xrange(len(u)): u[i] = string.atof(u[i])
    UnitCell = uctbx.UnitCell(u)
    print "Unit cell parameters:", UnitCell
    print
    ShowInputSymbol(inp.sgsymbol, inp.convention, "Input ")
    Symbols_Inp = sgtbx.SpaceGroupSymbols(inp.sgsymbol, inp.convention)
    SgOps = HallSymbol_to_SgOps(Symbols_Inp.Hall())
    SgInfo = SgOps.Info()
    print "Space group: (%d) %s" % (
      SgInfo.SgNumber(), SgInfo.BuildLookupSymbol())
    print

    SgOps.CheckUnitCell(UnitCell)

    d_min = string.atof(inp.d_min)
    print "Minimum d-spacing:", d_min
    if (d_min <= 0.):
      raise ValueError, "d-spacing must be greater than zero."
    print

    MinMateDistance = string.atof(inp.MinMateDistance)
    print "Minimum distance between symmetry mates: ", MinMateDistance

    SnapParameters = \
      sgtbx.SpecialPositionSnapParameters(UnitCell, SgOps, 1, MinMateDistance)
    WyckoffTable = sgtbx.WyckoffTable(SgInfo, 1)

    print "</pre><table border=2 cellpadding=2>"
    InTable = 1
    print "<tr>"
    print "<th>Label"
    print "<th>Scattering<br>factor<br>label"
    print "<th>Multiplicty"
    print "<th>Wyckoff<br>position"
    print "<th>Site<br>symmetry"
    print "<th colspan=3>Fractional coordinates"
    print "<th>Occupancy<br>factor"
    print "<th>Biso"
    print "<tr>"
    Sites = []
    print
    for line in inp.coordinates:
      flds = string.split(line)
      site = SiteInfo(flds)
      if (inp.coor_type != "Fractional"):
        site.Coordinates = UnitCell.fractionalize(site.Coordinates)
      SS = sgtbx.SiteSymmetry(SnapParameters, site.Coordinates, 0)
      site.WyckoffMapping = WyckoffTable.getWyckoffMapping(SS)
      site.SiteSymmetry = SS.PointGroupType()
      Sites.append(site)
      print "<tr>"
      print (  "<td>%s<td>%s"
             + "<td align=center>%d<td align=center>%s<td align=center>%s"
             + "<td><tt>%.6g</tt><td><tt>%.6g</tt><td><tt>%.6g</tt>"
             + "<td align=center><tt>%.6g</tt>"
             + "<td align=center><tt>%.6g</tt>") % (
      (site.Label, site.Sf.Label(),
       site.WyckoffMapping.WP().M(), site.WyckoffMapping.WP().Letter(),
       site.SiteSymmetry)
      + tuple(site.Coordinates) + (site.Occ, site.Biso))
    print "</table><pre>"
    InTable = 0
    print

    MillerIndices = BuildMillerIndices(UnitCell, SgInfo, d_min)
    FcalcDict = ComputeStructureFactors(Sites, MillerIndices)

    print "Number of Miller indices:", len(FcalcDict)
    print
    print "</pre><table border=2 cellpadding=2>"
    InTable = 1
    print "<tr>"
    print "<th>hkl<th>Amplitude<th>Phase"
    for H in FcalcDict.keys():
      print "<tr>"
      print "<td>%3d %3d %3d<td>%.6g<td align=right>%.3f" % (
        tuple(H) + polar(FcalcDict[H]))
    print "</table><pre>"
    InTable = 0
    print

  except RuntimeError, e:
    if (InTable): print "</table><pre>"
    print e
  except:
    if (InTable): print "</table><pre>"
    ei = sys.exc_info()
    print traceback.format_exception_only(ei[0], ei[1])[0]

  print "</pre>"
