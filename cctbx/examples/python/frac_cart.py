#! /usr/local_cci/Python-2.1.1/bin/python

# Convert coordinates from/to Fractional or Orthogonal.

PATH_cctbx_lib_python = "/net/boa/srv/html/cci/cctbx/lib_python"

import sys
sys.stderr = sys.stdout

print "Content-type: text/plain"
print

import traceback
import exceptions
class FormatError(exceptions.Exception): pass

import string, cgi

sys.path.insert(0, PATH_cctbx_lib_python)
from cctbx import uctbx

print "uctbx version:", uctbx.__version__
print

class Empty: pass

def GetFormData():
  form = cgi.FieldStorage()
  inp = Empty()
  for key in (("ucparams", "1 1 1 90 90 90"),
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

  if (inp.coor_type == "Fractional"):
    print "Cartesian coordinates:"
  else:
    print "Fractional coordinates:"
  print

  skip_columns = string.atoi(inp.skip_columns)
  if (skip_columns < 0):
    raise FormatError, "Negative number for columns to skip."
  for line in inp.coordinates:
    skipped, coordinates = InterpretCoordinateLine(line, skip_columns)
    if (inp.coor_type == "Fractional"):
      c = UnitCell.orthogonalize(coordinates)
    else:
      c = UnitCell.fractionalize(coordinates)
    print skipped, "%.6g %.6g %.6g" % tuple(c)

except RuntimeError, e:
  print e
except:
  ei = sys.exc_info()
  print traceback.format_exception_only(ei[0], ei[1])[0]
