#! /usr/local_cci/Python-2.2.1/bin/python

# Revision history:
#   2001 Nov: Use sftbx (rwgk).
#   2001 Sep: SpaceGroupType -> SpaceGroupInfo (R.W. Grosse-Kunstleve)
#   2001 Aug: Derived from generate_hklf.py (rwgk)

PATH_cctbx_lib_python = "/net/boa/srv/html/cci/cctbx/lib_python"

import sys
sys.stderr = sys.stdout

print "Content-type: text/html"
print

import traceback
import exceptions
class FormatError(exceptions.Exception): pass

import math, string, cgi

sys.path.insert(0, PATH_cctbx_lib_python)
from cctbx_boost.arraytbx import flex
from cctbx_boost import sgtbx
from cctbx_boost import uctbx
from cctbx_boost.eltbx.caasf_wk1995 import CAASF_WK1995
from cctbx_boost import adptbx
from cctbx_boost import miller
from cctbx_boost import sftbx

print "sgtbx version:", sgtbx.__version__
print "<br>"
print "uctbx version:", uctbx.__version__
print "<br>"
print "adptbx version:", adptbx.__version__
print "<br>"
print "sftbx version:", sftbx.__version__
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
        self.SF = CAASF_WK1995(flds[1], 1)
      else:
        offs = 1
        self.SF = CAASF_WK1995(flds[0], 0)
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
    except:
      raise FormatError, flds

  def as_XrayScatterer(self, UnitCell, SpaceGroup, MinMateDistance):
    Scatterer = sftbx.XrayScatterer(
      self.Label, self.SF, 0j, self.Coordinates,
      self.Occ, adptbx.B_as_U(self.Biso))
    Scatterer.ApplySymmetry(UnitCell, SpaceGroup, MinMateDistance, 0.1, 1)
    return Scatterer

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
    Sites = flex.XrayScatterer()
    print
    for line in inp.coordinates:
      flds = string.split(line)
      Site = SiteInfo(flds)
      if (inp.coor_type != "Fractional"):
        Site.Coordinates = UnitCell.fractionalize(Site.Coordinates)
      Site = Site.as_XrayScatterer(UnitCell, SgOps, MinMateDistance)
      Sites.append(Site)
      SS = sgtbx.SiteSymmetry(SnapParameters, Site.Coordinates())
      WyckoffPosition = WyckoffTable.getWyckoffMapping(SS).WP()
      SiteSymmetryPointGroupLabel = SS.PointGroupType()
      print "<tr>"
      print (  "<td>%s<td>%s"
             + "<td align=center>%d<td align=center>%s<td align=center>%s"
             + "<td><tt>%.6g</tt><td><tt>%.6g</tt><td><tt>%.6g</tt>"
             + "<td align=center><tt>%.6g</tt>"
             + "<td align=center><tt>%.6g</tt>") % (
        (Site.Label(), Site.CAASF().Label(),
         WyckoffPosition.M(), WyckoffPosition.Letter(),
         SiteSymmetryPointGroupLabel)
       + Site.Coordinates()
       + (Site.Occ(), adptbx.U_as_B(Site.Uiso())))
    print "</table><pre>"
    InTable = 0
    print

    MillerIndices = miller.BuildIndices(UnitCell, SgInfo, 1, d_min)
    Fcalc = sftbx.StructureFactorArray(UnitCell, SgOps, MillerIndices, Sites)

    print "Number of Miller indices:", len(Fcalc)
    print
    print "</pre><table border=2 cellpadding=2>"
    InTable = 1
    print "<tr>"
    print "<th>hkl<th>Amplitude<th>Phase"
    for i in xrange(len(MillerIndices)):
      print "<tr>"
      print "<td>%3d %3d %3d<td>%.6g<td align=right>%.3f" % (
        MillerIndices[i] + polar(Fcalc[i]))
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
