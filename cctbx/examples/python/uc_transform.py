#! /usr/local/Python-2.1/bin/python

import sys
sys.stderr = sys.stdout

print "Content-type: text/plain"
print

import exceptions
class InternalError(exceptions.Exception):
  pass

import string, re, cgi

sys.path.insert(0, "/net/boa/srv/html/sgtbx") # for sgtbx, uctbx
import sgtbx
import uctbx

print "sgtbx version:", sgtbx.__version__
print "uctbx version:", uctbx.__version__
print

class Empty: pass

def GetFormData():
  form = cgi.FieldStorage()
  inp = Empty()
  for key in (("ucparams", "1 1 1 90 90 90"),
              ("sgsymbol_1", "P1"),
              ("convention_1", ""),
              ("sgsymbol_2", ""),
              ("convention_2", "")):
    if (form.has_key(key[0])):
      inp.__dict__[key[0]] = string.strip(form[key[0]].value)
    else:
      inp.__dict__[key[0]] = key[1]
  return inp

def ShowInputSymbol(sgsymbol, convention):
  if (sgsymbol != ""):
    print "Input space group symbol:", sgsymbol
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

inp = GetFormData()

try:
  u = string.split(inp.ucparams)
  for i in xrange(len(u)): u[i] = string.atof(u[i])
  UnitCell = uctbx.UnitCell(u)
  print "Unit cell parameters:", UnitCell
  print
  ShowInputSymbol(inp.sgsymbol_1, inp.convention_1)
  SgOps_1 = Symbol_to_SgOps(inp.sgsymbol_1, inp.convention_1)
  SgType_1 = SgOps_1.getSpaceGroupType()
  print

  if (len(inp.sgsymbol_2) == 0):
    inp.sgsymbol_2 = sgtbx.SpaceGroupSymbols(
      SgType_1.SgNumber()).ExtendedHermann_Mauguin()
    inp.convention_2 = ""

  ShowInputSymbol(inp.sgsymbol_2, inp.convention_2)
  SgOps_2 = Symbol_to_SgOps(inp.sgsymbol_2, inp.convention_2)
  SgType_2 = SgOps_2.getSpaceGroupType()
  print

  print "Space group 1: (%d) %s" % (
    SgType_1.SgNumber(), SgOps_1.BuildLookupSymbol(SgType_1))
  print "Space group 2: (%d) %s" % (
    SgType_2.SgNumber(), SgOps_2.BuildLookupSymbol(SgType_2))
  print

  if (SgType_1.SgNumber() != SgType_2.SgNumber()):
    print "Space group numbers are not equal!"
  else:
    M = SgType_2.CBOp().InvM() * SgType_1.CBOp().M()
    CBOp = sgtbx.ChOfBasisOp(M)
    print "Change-of-basis matrix:", CBOp.M()
    print "               Inverse:", CBOp.InvM()
    print
    SgOps_1.CheckUnitCell(UnitCell)
    NewUnitCell = UnitCell.ChangeBasis(CBOp.InvM().as_tuple()[0])
    print "Transformed unit cell parameters: ", NewUnitCell
    SgOps_2.CheckUnitCell(NewUnitCell)

except RuntimeError, e:
  print e
