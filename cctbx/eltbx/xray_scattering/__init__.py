from __future__ import division
import scitbx.math.gaussian # base class for gaussian

import boost.python
ext = boost.python.import_ext("cctbx_eltbx_xray_scattering_ext")
from cctbx_eltbx_xray_scattering_ext import *

import sys

# grep 'DATA KQ' shelxl.f | cut -d"'" -f2 | grep -v NCSY
shelxl_97_2_980324_tabulated_chemical_elements = """\
H HE LI BE B C N O F NE NA MG AL SI P S CL AR K CA SC TI V CR MN FE
CO NI CU ZN GA GE AS SE BR KR RB SR Y ZR NB MO TC RU RH PD AG CD IN
SN SB TE I XE CS BA LA CE PR ND PM SM EU GD TB DY HO ER TM YB LU HF
TA W RE OS IR PT AU HG TL PB BI PO AT RN FR RA AC TH PA U NP PU""".split()

def get_element_and_charge_symbols(scattering_type, exact=True):
  sl = get_standard_label(label=scattering_type, exact=exact, optional=True)
  if (sl is None): return "", ""
  if (sl == "Hiso"): return "H", ""
  if (sl == "Cval"): return "C", ""
  if (sl == "Sival"): return "Si", ""
  if (sl[-1] in ["+", "-"]):
    return sl[:-2], sl[-2:]
  return sl, ""

class _(boost.python.injector, ext.gaussian):

  def show(self, f=None, format=None):
    if (f is None): f = sys.stdout
    if (format is None): format = "%.8g"
    for l,v in (("a:", self.array_of_a()), ("b:", self.array_of_b())):
      print >> f, l, " ".join([format % x for x in v])
    print >> f, "c:", format % self.c()
    return self

  def electron_density(self, r, b_iso):
    from math import pi, exp
    result = 0
    def ft(b):
      # Agarwal (1978). Acta Cryst. A34, 791-809.
      # Page 796 before equation (42).
      return (4*pi/(b+b_iso))**(3/2) * exp(-4*pi**2*r**2/(b+b_iso))
    for a,b in zip(self.array_of_a(), self.array_of_b()):
      result += a * ft(b)
    if (self.use_c()):
      result += self.c() * ft(0)
    return result

def best_approximation(scattering_type):
  if (scattering_type == "const"):
    return gaussian(1)
  return wk1995(scattering_type, True).fetch()

class two_gaussian_agarwal_isaacs(object):
  source="ccp4/lib/data/atomsf.lib Revision 1.4, Thu Feb 13 14:10:58 1997 UTC"
  source_short = "CCP4 atomsf.lib Rev. 1.4"
  table = {
    "H": gaussian([0.7932, 0.1949], [24.2157, 2.1089]),
    "C": gaussian([2.9972, 2.9791], [30.016701, 2.8886]),
    "N": gaussian([2.9924, 3.9986], [25.3766, 3.5004]),
    "O": gaussian([2.4485, 5.5589], [24.756199, 4.1372]),
    "S": gaussian([5.5480, 10.4241], [33.7108, 1.9034]),
  }

class two_gaussian_agarwal_1978(object):
  source = "Agarwal, R.C. (1978). Acta Cryst. A34, 791-809, Table 1."
  source_short = "Agarwal (1978)"
  table = {
    "H": gaussian([0.4866, 0.5098], [34.1284, 8.8996]),
    "C": gaussian([3.0102, 2.9705], [29.9132, 2.8724]),
    "N": gaussian([3.0492, 3.9432], [25.0383, 3.4059]),
    "O": gaussian([3.2942, 4.6968], [20.0401, 3.1184]),
    "S": gaussian([5.6604, 10.3140], [33.0400, 1.8160]),
    "Fe3+": gaussian([10.3568, 12.6329], [8.1324, 0.8137]),
    "Fe2+": gaussian([11.6635, 12.3057], [9.0361, 0.5749]),
    "Zn2+": gaussian([5.7826, 22.2163], [11.7082, 1.8234]),
    "Ba2+": gaussian([12.1432, 41.8442], [21.7090, 1.4090]),
  }

class one_gaussian_agarwal_1978(object):
  source = "Agarwal, R.C. (1978). Acta Cryst. A34, 791-809, Table 3."
  source_short = "Agarwal (1978)"
  table = {
    "C": gaussian([5.9074], [1.2913]),
    "N": gaussian([7.0411], [0.2065]),
    "O": gaussian([8.1561], [-0.8941]),
    "S": gaussian([15.8448], [-2.1392]),
  }

class fitted_gaussian(gaussian):

  def __init__(self, stol, gaussian_sum, max_error=None):
    gaussian.__init__(self, gaussian_sum)
    self.stol = stol
    self.max_error = max_error

  def __getinitargs__(self):
    return (self.stol, gaussian(self), self.max_error)

  def sort(self):
    return fitted_gaussian(self.stol, gaussian.sort(self), self.max_error)

  def show(self, f=None, format=None):
    if (f is None): f = sys.stdout
    if (self.max_error is None):
      e = ""
    else:
      e = ", max_error: %.4f" % self.max_error
    print >> f, "stol: %.2f # d_min: %.2f%s" % (self.stol, 1/(2*self.stol), e)
    return gaussian.show(self, f, format)
