import cctbx.array_family.flex # for tuple mappings

from scitbx.python_utils import misc
ext = misc.import_ext("cctbx_boost.eltbx.xray_scattering_ext")
misc.import_regular_symbols(globals(), ext.__dict__)
del misc

from scitbx.boost_python_utils import injector
import sys

class _gaussian(injector, ext.gaussian):

  def show(self, f=None, format=None):
    if (f is None): f = sys.stdout
    if (format is None): format = "%.8g"
    for l,v in (("a:", self.a()), ("b:", self.b())):
      print >> f, l, " ".join([format % x for x in v])
    print >> f, "c:", format % self.c()
    return self

def best_approximation(scattering_type):
  if (scattering_type == "const"):
    return gaussian(1)
  return wk1995(scattering_type, 1).fetch()

class two_gaussian_agarwal_isaacs:
  source="ccp4/lib/data/atomsf.lib Revision 1.4, Thu Feb 13 14:10:58 1997 UTC"
  source_short = "CCP4 atomsf.lib Rev. 1.4"
  table = {
    "H": gaussian([0.7932, 0.1949], [24.2157, 2.1089]),
    "C": gaussian([2.9972, 2.9791], [30.016701, 2.8886]),
    "N": gaussian([2.9924, 3.9986], [25.3766, 3.5004]),
    "O": gaussian([2.4485, 5.5589], [24.756199, 4.1372]),
    "S": gaussian([5.5480, 10.4241], [33.7108, 1.9034]),
  }

class two_gaussian_agarwal_1978:
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

class one_gaussian_agarwal_1978:
  source = "Agarwal, R.C. (1978). Acta Cryst. A34, 791-809, Table 3."
  source_short = "Agarwal (1978)"
  table = {
    "C": gaussian([5.9074], [1.2913]),
    "N": gaussian([7.0411], [0.2065]),
    "O": gaussian([8.1561], [-0.8941]),
    "S": gaussian([15.8448], [-2.1392]),
  }

class fitted_gaussian(gaussian):

  def __init__(self, stol, a, b, c=0):
    gaussian.__init__(self, a, b, c)
    self.stol = stol

  def __getinitargs__(self):
    return (self.stol, self.a(), self.b(), self.c())

  def show(self, f=None, format=None):
    if (f is None): f = sys.stdout
    print "stol: %.2f # d_min: %.2f" % (self.stol, 1/(2*self.stol))
    gaussian.show(self, f, format)
