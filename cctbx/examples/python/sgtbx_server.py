#! /usr/local/Python-2.1/bin/python

import sys
sys.stderr = sys.stdout

print "Content-type: text/plain"
print

import exceptions
class InternalError(exceptions.Exception):
  pass

import string, re, cgi

sys.path.insert(0, "/net/boa/srv/html/sgtbx") # for sgtbx
import sgtbx

print "sgtbx version:", sgtbx.__version__
print

def GetFormData():
  form = cgi.FieldStorage()

  sgsymbol = ""
  if (form.has_key("sgsymbol")):
    sgsymbol = string.strip(form["sgsymbol"].value)

  convention = ""
  if (form.has_key("convention")):
    convention = form["convention"].value

  symxyz = []
  if (form.has_key("symxyz")):
    lines = string.split(form["symxyz"].value, "\015\012")
    for l in lines:
      for s in string.split(l, ";"):
        s = string.strip(s)
        if (s != ""): symxyz.append(string.strip(s))

  return (sgsymbol, convention, symxyz)


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

def RTMxAnalysis(M):
  RI = M.getRotMxInfo()
  TI = M.analyzeTpart()
  if (RI.Rtype() == 1):
    return "rotation=1"
  elif (RI.Rtype() == -1):
    return "rotation=%d shift=(%s)" % (
      RI.Rtype(), TI.OriginShift().as_xyz(),)
  elif (abs(RI.Rtype()) == 2):
    return "rotation=%d%s intrinsic=(%s) shift=(%s)" % (
      RI.Rtype(), StrEV(RI.EV()),
      TI.IntrinsicPart().as_xyz(), TI.OriginShift().as_xyz())
  else:
    return "rotation=%d^%d%s intrinsic=(%s) shift=(%s)" % (
       RI.Rtype(), RI.SenseOfRotation(), StrEV(RI.EV()),
       TI.IntrinsicPart().as_xyz(), TI.OriginShift().as_xyz())

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
  print "List of symmetry operations:"
  for M in SgOps:
    print " ", M.as_xyz(), "", RTMxAnalysis(M)
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


(sgsymbol, convention, symxyz) = GetFormData()
#sgsymbol = "P 4"
#symxyz = ["x,y,z", "x,y,z"]
ShowInputSymbol(sgsymbol, convention)

try:
  Symbols_Inp = None
  lookup_symbol = sgsymbol
  if (lookup_symbol == ""): lookup_symbol = "P 1"
  if (convention == "Hall"):
    HallSymbol = lookup_symbol
  else:
    Symbols_Inp = sgtbx.SpaceGroupSymbols(lookup_symbol, convention)
    HallSymbol = Symbols_Inp.Hall()
    if (Symbols_Inp.SgNumber() == 0):
      Symbols_Inp = None
      convention = "Hall"
    else:
      print "Result of symbol lookup:"
      ShowSymbols(Symbols_Inp)
      print

  try:
    ps = sgtbx.parse_string(HallSymbol)
    SgOps = sgtbx.SgOps(ps)
  except RuntimeError, e:
    print "-->" + ps.string() + "<--"
    print ("-" * (ps.where() + 3)) + "^"
    raise

  if (len(symxyz) != 0):
    print "Addition of symmetry operations:"
    for s in symxyz:
      ps = sgtbx.parse_string(s)
      try:
        M = sgtbx.RTMx(ps)
      except RuntimeError, e:
        print "-->" + ps.string() + "<--"
        print ("-" * (ps.where() + 3)) + "^"
        raise
      print " ", M.as_xyz(), "", RTMxAnalysis(M)
      SgOps.expandSMx(M)
    print

  ShowSgOpsGeneric(SgOps)

  if (convention == "Hall" or len(symxyz) != 0):
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
  print e
