#! /usr/local/Python-2.1/bin/python

PATH_cctbx_lib_python = "/net/boa/srv/html/cci/cctbx/lib/python"

import sys
sys.stderr = sys.stdout

print "Content-type: text/html"
print

import exceptions
class InternalError(exceptions.Exception):
  pass

import string, re, cgi

sys.path.insert(0, PATH_cctbx_lib_python)
import sgtbx

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
      RI.Rtype(), TI.OriginShift().as_xyz(),)
  elif (abs(RI.Rtype()) == 2):
    print "<td>%d<td>%s<td>(%s)<td>(%s)" % (
      RI.Rtype(), StrEV(RI.EV()),
      TI.IntrinsicPart().as_xyz(), TI.OriginShift().as_xyz())
  else:
    print "<td>%d^%d<td>%s<td>(%s)<td>(%s)" % (
       RI.Rtype(), RI.SenseOfRotation(), StrEV(RI.EV()),
       TI.IntrinsicPart().as_xyz(), TI.OriginShift().as_xyz())
  print "</tr>"

def ShowSgOpsGeneric(SgOps):
  print "Number of lattice translations:", SgOps.nLTr()
  if (SgOps.isCentric()):
    print "Space group is acentric."
  else:
    print "Space group is centric."
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
      SgType = SgOps.getSpaceGroupType(1)
      print "Space group number:", SgType.SgNumber()
      print "Conventional Hermann-Mauguin symbol:", \
        sgtbx.SpaceGroupSymbols(SgType.SgNumber()).ExtendedHermann_Mauguin()
      print "Hall symbol:", SgOps.BuildHallSymbol(SgType, 1)
      print "Change-of-basis matrix:", SgType.CBOp().M().as_xyz()
      print "               Inverse:", SgType.CBOp().InvM().as_xyz()
      print

# Not implemented.
#  print "Additional generators of Euclidean normalizer:"
#  ssVM = SgOps.get_ss()
#  print "  Number of structure-seminvariant vectors and moduli:", ssVM["N"]
#  print "    vector        modulus"
#  for i in xrange(ssVM["N"]):
#    v = ssVM["VM"][i][0]
#    m = ssVM["VM"][i][1]
#    print "    (%2d, %2d, %2d)  %d" % (v[0], v[1], v[2], m)
#  AddlG = SgOps.get_AddlGenEuclNorm(K2L=1, L2N=1)
#  print "  Number of rotational generators:", AddlG["N"]
#  for i in xrange(AddlG["N"]):
#    print "   ", sglite.RTMx2XYZ(RTMx=AddlG["SMx"][i],
#                                 RBF=sglite.SRBF, TBF=sglite.STBF)
#  print

except RuntimeError, e:
  if (InTable): print "</table><pre>"
  print e

print "</pre>"
