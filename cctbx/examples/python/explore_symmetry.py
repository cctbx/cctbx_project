#! /usr/local/Python-2.1/bin/python

PATH_cctbx_lib_python = "/net/boa/srv/html/cci/cctbx/lib/python"

import sys
sys.stderr = sys.stdout

print "Content-type: text/html"
print

import exceptions
class InternalError(exceptions.Exception):
  pass

import string, re, cgi, math

sys.path.insert(0, PATH_cctbx_lib_python)
import sgtbx
import uctbx

print "sgtbx version:", sgtbx.__version__
print "<p>"
print "<pre>"
InTable = 0

class Empty: pass

def GetFormData():
  form = cgi.FieldStorage()
  inp = Empty()
  for key in (("sgsymbol", ""),
              ("convention", "")):
    if (form.has_key(key[0])):
      inp.__dict__[key[0]] = string.strip(form[key[0]].value)
    else:
      inp.__dict__[key[0]] = key[1]
  inp.symxyz = []
  if (form.has_key("symxyz")):
    lines = string.split(form["symxyz"].value, "\015\012")
    for l in lines:
      for s in string.split(l, ";"):
        s = string.strip(s)
        if (s != ""): inp.symxyz.append(string.strip(s))
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
  print

def StrEV(EV):
  return "[%d,%d,%d]" % EV

def RTMxAnalysisHeader():
  print "<tr>"
  print "<th>Matrix"
  print "<th>Rotation-part type"
  print "<th>Axis direction"
  print "<th>Screw/glide component"
  print "<th>Origin shift"
  print "</tr>"

def RTMxAnalysis(M):
  print "<tr>"
  print "<td><tt>" + str(M) + "</tt>"
  RI = M.getRotMxInfo()
  TI = M.analyzeTpart()
  if (RI.Rtype() == 1):
    print "<td>1<td>-<td>-<td>-"
  elif (RI.Rtype() == -1):
    print "<td>%d<td>-<td>-<td>(%s)" % (
      RI.Rtype(), TI.OriginShift())
  elif (abs(RI.Rtype()) == 2):
    print "<td>%d<td>%s<td>(%s)<td>(%s)" % (
      RI.Rtype(), StrEV(RI.EV()),
      TI.IntrinsicPart(), TI.OriginShift())
  else:
    print "<td>%d^%d<td>%s<td>(%s)<td>(%s)" % (
       RI.Rtype(), RI.SenseOfRotation(), StrEV(RI.EV()),
       TI.IntrinsicPart(), TI.OriginShift())
  print "</tr>"

def ShowSgOpsGeneric(SgOps):
  print "Number of lattice translations:", SgOps.nLTr()
  if (SgOps.isCentric()):
    print "Space group is centric."
  else:
    print "Space group is acentric."
  if (SgOps.isChiral()):
    print "Space group is chiral."
  if (SgOps.isEnantiomorphic()):
    print "Space group is enantiomorphic."
  print "Number of representative symmetry operations:", SgOps.nSMx()
  print "Total number of symmetry operations:", SgOps.OrderZ()
  print
  print "Parallelepiped containing an asymmetric unit:"
  try: Brick = SgOps.getBrick()
  except RuntimeError, e:
    print " ", e
  else:
    print " ", Brick
  print
  print "List of symmetry operations:"
  print "</pre><table border=2 cellpadding=2>"
  global InTable
  InTable = 1
  RTMxAnalysisHeader()
  for M in SgOps: RTMxAnalysis(M)
  print "</table><pre>"
  InTable = 0
  print


def ShowSymbols(Symbols):
  print "  Space group number:", Symbols.SgNumber()
  print "  Schoenflies symbol:", Symbols.Schoenflies()
  print "  Hermann-Mauguin symbol:", Symbols.Hermann_Mauguin()
  E = Symbols.Extension()
  if (E != ""):
    if (E in "12"):
      print "  Origin choice:", E
    elif (E == "H"):
      print "  Trigonal using hexagonal axes"
    elif (E == "R"):
      print "  Trigonal using rhombohedral axes"
    else:
      raise InternalError
  Q = Symbols.Qualifier()
  if (Q != ""):
    if (Symbols.SgNumber() < 16):
      if (Q[-1] in "123"):
        UniqueAxis = Q[:-1]
        CellChoice = Q[-1]
      else:
        UniqueAxis = Q
        CellChoice = ""
      print "  Unique axis:", UniqueAxis
      if (CellChoice != ""):
        print "  Cell choice:", CellChoice
    else:
      print "  Relation to standard setting:", Q
  print "  Hall symbol:", string.strip(Symbols.Hall())

def get_unitcell(SgType):
  if (143 <= SgType.SgNumber() < 195):
    RefUnitCell = uctbx.UnitCell((10, 10, 10, 90, 90, 120))
  else:
    RefUnitCell = uctbx.UnitCell((10, 10, 10, 90, 90, 90))
  return RefUnitCell.ChangeBasis(SgType.CBOp().M().as_tuple()[0])

inp = GetFormData()

print "<pre>"
ShowInputSymbol(inp.sgsymbol, inp.convention)

try:
  Symbols_Inp = None
  lookup_symbol = inp.sgsymbol
  if (lookup_symbol == ""): lookup_symbol = "P 1"
  if (inp.convention == "Hall"):
    HallSymbol = lookup_symbol
  else:
    Symbols_Inp = sgtbx.SpaceGroupSymbols(lookup_symbol, inp.convention)
    HallSymbol = Symbols_Inp.Hall()
    if (Symbols_Inp.SgNumber() == 0):
      Symbols_Inp = None
      inp.convention = "Hall"
    else:
      print "Result of symbol lookup:"
      ShowSymbols(Symbols_Inp)
      print

  try:
    ps = sgtbx.parse_string(HallSymbol)
    SgOps = sgtbx.SgOps(ps)
  except RuntimeError, e:
    print "--&gt;" + ps.string() + "&lt;--"
    print ("-" * (ps.where() + 3)) + "^"
    raise

  if (len(inp.symxyz) != 0):
    print "Addition of symmetry operations:"
    print "</pre><table border=2 cellpadding=2>"
    InTable = 1
    RTMxAnalysisHeader()
    for s in inp.symxyz:
      ps = sgtbx.parse_string(s)
      try:
        M = sgtbx.RTMx(ps)
      except RuntimeError, e:
        print "</table><pre>"
        InTable = 0
        print "--&gt;" + ps.string() + "&lt;--"
        print ("-" * (ps.where() + 3)) + "^"
        raise
      RTMxAnalysis(M)
      SgOps.expandSMx(M)
    print "</table><pre>"
    InTable = 0
    print

  ShowSgOpsGeneric(SgOps)
  SgType = SgOps.getSpaceGroupType(1)

  if (inp.convention == "Hall" or len(inp.symxyz) != 0):
    Symbols_Match = SgOps.MatchTabulatedSettings()
    if (Symbols_Match.SgNumber() != 0):
      if (   Symbols_Inp == None
          or    Symbols_Inp.ExtendedHermann_Mauguin()
             != Symbols_Match.ExtendedHermann_Mauguin()):
        print "Symmetry operations match:"
        ShowSymbols(Symbols_Match)
        print
      else:
        print "Additional symmetry operations are redundant."
        print
    else:
      print "Space group number:", SgType.SgNumber()
      print "Conventional Hermann-Mauguin symbol:", \
        sgtbx.SpaceGroupSymbols(SgType.SgNumber()).ExtendedHermann_Mauguin()
      print "Hall symbol:", SgOps.BuildHallSymbol(SgType, 1)
      print "Change-of-basis matrix:", SgType.CBOp().M()
      print "               Inverse:", SgType.CBOp().InvM()
      print

  UnitCell = get_unitcell(SgType)
  SnapParameters = sgtbx.SpecialPositionSnapParameters(
    UnitCell, SgOps, 1, 1.e-6)
  WyckoffTable = sgtbx.WyckoffTable(SgType)
  print "List of Wyckoff positions:"
  print "</pre><table border=2 cellpadding=2>"
  InTable = 1
  print "<tr>"
  print "<th>Wyckoff letter"
  print "<th>Multiplicity"
  print "<th>Site symmetry<br>point group type"
  print "<th>Representative special position operator"
  print "</tr>"
  for WP in WyckoffTable:
    # The site symmetry point group type is not tabulated.
    # Generate dummy coordinates to obtain it the indirectly using
    # the SpecialPosition class.
    X = WP.SpecialOp().multiply((math.sqrt(2),math.sqrt(3),math.sqrt(5)))
    SP = sgtbx.SpecialPosition(SnapParameters, X, 0, 1)
    assert SP.DistanceMoved2() < 1.e-5
    WMap = WyckoffTable.getWyckoffMapping(SP)
    assert WMap.WP().Letter() == WP.Letter()
    print "<tr>"
    print "<td>%s<td>%d<td>%s<td><tt>%s</tt>" % (
      WP.Letter(), WP.M(), SP.getPointGroupType(), str(WP.SpecialOp()))
    print "</tr>"
  print "</table><pre>"
  InTable = 0
  print

  print "Additional generators of Euclidean normalizer:"
  ss = sgtbx.StructureSeminvariant(SgOps)
  print "  Number of structure-seminvariant vectors and moduli:", len(ss)
  if (len(ss)):
    print "    Vector    Modulus"
    for i in xrange(len(ss)): print "   ", ss.V(i), ss.M(i)
  K2L = SgType.getAddlGeneratorsOfEuclideanNormalizer(1, 0)
  L2N = SgType.getAddlGeneratorsOfEuclideanNormalizer(0, 1)
  if (len(K2L)):
    print "  Inversion through a centre at:",
    print K2L[0].analyzeTpart().OriginShift()
  if (len(L2N)):
    print "  Further generators:"
    print "</pre><table border=2 cellpadding=2>"
    InTable = 1
    RTMxAnalysisHeader()
    for M in L2N: RTMxAnalysis(M)
    print "</table><pre>"
    InTable = 0
  print

except RuntimeError, e:
  if (InTable): print "</table><pre>"
  print e

print "</pre>"
