#! /usr/local/Python-2.1/bin/python

PATH_cctbx_lib_python = "/net/boa/srv/html/cci/cctbx"

import sys
sys.stderr = sys.stdout

print "Content-type: text/plain"
print

import traceback
import exceptions
class FormatError(exceptions.Exception): pass

import string, cgi

sys.path.insert(0, PATH_cctbx_lib_python)
import sgtbx
import uctbx

print "sgtbx version:", sgtbx.__version__
print "uctbx version:", uctbx.__version__
print

class Empty: pass

def GetFormData():
  form = cgi.FieldStorage()
  inp = Empty()
  for key in (("ucparams_old", "1 1 1 90 90 90"),
              ("sgsymbol_old", "P1"),
              ("convention_old", ""),
              ("sgsymbol_new", ""),
              ("convention_new", ""),
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

def Symbol_to_SgOps(sgsymbol, convention):
  if (convention == "Hall"):
    HallSymbol = sgsymbol
  else:
    Symbols_Inp = sgtbx.SpaceGroupSymbols(sgsymbol, convention)
    HallSymbol = Symbols_Inp.Hall()
  try:
    ps = sgtbx.parse_string(HallSymbol)
    SgOps = sgtbx.SgOps(ps)
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
  u = string.split(inp.ucparams_old)
  for i in xrange(len(u)): u[i] = string.atof(u[i])
  UnitCell_old = uctbx.UnitCell(u)
  print "Old unit cell parameters:", UnitCell_old
  print
  ShowInputSymbol(inp.sgsymbol_old, inp.convention_old, "Old")
  SgOps_old = Symbol_to_SgOps(inp.sgsymbol_old, inp.convention_old)
  SgType_old = SgOps_old.getSpaceGroupType()
  print

  if (len(inp.sgsymbol_new) == 0):
    inp.sgsymbol_new = sgtbx.SpaceGroupSymbols(
      SgType_old.SgNumber()).ExtendedHermann_Mauguin()
    inp.convention_new = ""

  ShowInputSymbol(inp.sgsymbol_new, inp.convention_new, "New")
  SgOps_new = Symbol_to_SgOps(inp.sgsymbol_new, inp.convention_new)
  SgType_new = SgOps_new.getSpaceGroupType()
  print

  print "Old space group: (%d) %s" % (
    SgType_old.SgNumber(), SgOps_old.BuildLookupSymbol(SgType_old))
  print "New space group: (%d) %s" % (
    SgType_new.SgNumber(), SgOps_new.BuildLookupSymbol(SgType_new))
  print

  if (SgType_old.SgNumber() != SgType_new.SgNumber()):
    print "Space group numbers are not equal!"
  else:
    M = SgType_new.CBOp().InvM() * SgType_old.CBOp().M()
    CBOp = sgtbx.ChOfBasisOp(M)
    print "Change-of-basis matrix:", CBOp.M()
    print "               Inverse:", CBOp.InvM()
    print

    SgOps_old.CheckUnitCell(UnitCell_old)
    UnitCell_new = UnitCell_old.ChangeBasis(CBOp.InvM().as_tuple()[0])
    print "New unit cell parameters: ", UnitCell_new
    SgOps_new.CheckUnitCell(UnitCell_new)
    print

    print inp.coor_type, "coordinates:"
    print
    skip_columns = string.atoi(inp.skip_columns)
    if (skip_columns < 0):
      raise FormatError, "Negative number for columns to skip."
    for line in inp.coordinates:
      skipped, coordinates = InterpretCoordinateLine(line, skip_columns)
      if (inp.coor_type == "Fractional"):
        c = CBOp(coordinates)
      else:
        c = UnitCell_old.fractionalize(coordinates)
        c = CBOp(c)
        c = UnitCell_new.orthogonalize(c)
      print skipped, "%.6g %.6g %.6g" % tuple(c)

except RuntimeError, e:
  print e
except:
  ei = sys.exc_info()
  print traceback.format_exception_only(ei[0], ei[1])[0]
