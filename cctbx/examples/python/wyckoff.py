#! /usr/local/Python-2.1/bin/python

# When studying a crystal structure it can be helpful to know the Wyckoff
# positions of the atoms in the structure. This script can be used to
# assign Wyckoff letters.

PATH_cctbx_lib_python = "/net/boa/srv/html/cci/cctbx/lib_python"

import sys
sys.stderr = sys.stdout

print "Content-type: text/html"
print

import traceback
import exceptions
class FormatError(exceptions.Exception): pass

import string, cgi

sys.path.insert(0, PATH_cctbx_lib_python)
from cctbx import sgtbx
from cctbx import uctbx

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
              ("MinMateDistance", "0.5"),
              ("coor_type", None),
              ("skip_columns", "0")):
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

def InterpretCoordinateLine(line, skip_columns):
  flds = string.split(line)
  if (len(flds) < skip_columns + 3): raise FormatError, line
  coordinates = [0,0,0]
  for i in xrange(3):
    try: coordinates[i] = string.atof(flds[skip_columns + i])
    except: raise FormatError, line
  return string.join(flds[:skip_columns]), coordinates

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

  MinMateDistance = string.atof(inp.MinMateDistance)
  SnapParameters = \
    sgtbx.SpecialPositionSnapParameters(UnitCell, SgOps, 1, MinMateDistance)
  WyckoffTable = sgtbx.WyckoffTable(SgInfo, 1)

  skip_columns = string.atoi(inp.skip_columns)
  if (skip_columns < 0):
    raise FormatError, "Negative number for columns to skip."

  print "</pre><table border=2 cellpadding=2>"
  InTable = 1
  print "<tr>"
  if (skip_columns): print "<th>"
  print "<th colspan=3>" + inp.coor_type + " coordinates"
  print "<th>Multiplicity"
  print "<th>Wyckoff letter"
  print "<th>Site symmetry<br>point group type"
  print "<th>Special position operator"
  print "</tr>"
  for line in inp.coordinates:
    skipped, coordinates = InterpretCoordinateLine(line, skip_columns)
    if (inp.coor_type != "Fractional"):
      coordinates = UnitCell.fractionalize(coordinates)
    SS = sgtbx.SiteSymmetry(SnapParameters, coordinates, 0)
    SnapPosition = SS.SnapPosition()
    WyckoffMapping = WyckoffTable.getWyckoffMapping(SS)
    if (inp.coor_type != "Fractional"):
      SnapPosition = UnitCell.orthogonalize(SnapPosition)
    print "<tr>"
    if (skip_columns): print "<td>", skipped
    for elem in SnapPosition: print "<td><tt>%.6g</tt>" % (elem,)
    print "<td align=center>", WyckoffMapping.WP().M()
    print "<td align=center>", WyckoffMapping.WP().Letter()
    print "<td align=center>", SS.PointGroupType()
    print "<td><tt>" + str(SS.SpecialOp()) + "</tt>"
    print "</tr>"
  print "</table><pre>"
  InTable = 0

except RuntimeError, e:
  if (InTable): print "</table><pre>"
  print e
except:
  if (InTable): print "</table><pre>"
  ei = sys.exc_info()
  print traceback.format_exception_only(ei[0], ei[1])[0]

print "</pre>"
