#! /usr/local/Python-2.1/bin/python

PATH_cctbx_lib_python = "/net/boa/srv/html/cci/cctbx/lib/python"
URL_explore_symmetry = "http://cci.lbl.gov/cctbx/explore_symmetry.py"

import sys
sys.stderr = sys.stdout

print "Content-type: text/html"
print

import traceback
import exceptions
class FormatError(exceptions.Exception): pass

import string, cgi, urllib

sys.path.insert(0, PATH_cctbx_lib_python)
import sgtbx

print "sgtbx version:", sgtbx.__version__
print "<p>"

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
  return inp

def HallSymbol_to_SgOps(HallSymbol):
  try:
    ps = sgtbx.parse_string(HallSymbol)
    SgOps = sgtbx.SpaceGroup(ps)
  except RuntimeError, e:
    print "-->" + ps.string() + "<--"
    print ("-" * (ps.where() + 3)) + "^"
    raise
  return SgOps

inp = GetFormData()

try:
  SgNumber = 0
  if (len(string.strip(inp.sgsymbol)) != 0):
    Symbols_Inp = sgtbx.SpaceGroupSymbols(inp.sgsymbol, inp.convention)
    SgOps = HallSymbol_to_SgOps(Symbols_Inp.Hall())
    SgInfo = sgtbx.SpaceGroupInfo(SgOps)
    SgNumber = SgInfo.SgNumber()
  nSettings = 0
  print "<table border=2 cellpadding=2>"
  print "<tr>"
  print "<th>Space group<br>No."
  print "<th>Schoenflies<br>symbol"
  print "<th>Hermann-Mauguin<br>symbol"
  print "<th>Hall<br>symbol"
  for SgSymbols in sgtbx.SpaceGroupSymbolIterator():
    if (SgNumber == 0 or SgSymbols.SgNumber() == SgNumber):
      print "<tr>"
      print "<td>(%d)<td>%s" % (
        SgSymbols.SgNumber(), SgSymbols.Schoenflies())
      print "<td><a href=\"%s?sgsymbol=%s\">%s</a>" % (
        URL_explore_symmetry,
        urllib.quote_plus(SgSymbols.ExtendedHermann_Mauguin()),
        SgSymbols.ExtendedHermann_Mauguin())
      print "<td>%s" % (SgSymbols.Hall(),)
      nSettings = nSettings + 1
  print "</table>"
  if (SgNumber == 0):
    print "<p>"
    print "Number of settings listed:", nSettings

except RuntimeError, e:
  print e
except:
  ei = sys.exc_info()
  print traceback.format_exception_only(ei[0], ei[1])[0]
