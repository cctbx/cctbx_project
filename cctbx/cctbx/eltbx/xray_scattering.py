import cctbx.array_family.flex # for tuple mappings

from scitbx.python_utils import misc
ext = misc.import_ext("cctbx_boost.eltbx.xray_scattering_ext")
misc.import_regular_symbols(globals(), ext.__dict__)
del misc

def best_approximation(scattering_type):
  if (scattering_type == "const"):
    return gaussian(1)
  return wk1995(scattering_type, 1).fetch()

class two_gaussian_agarwal_isaacs:
  source="ccp4/lib/data/atomsf.lib Revision 1.4, Thu Feb 13 14:10:58 1997 UTC"
  table = {
    "H": gaussian([0.7932, 0.1949], [24.2157, 2.1089], 0),
    "C": gaussian([2.9972, 2.9791], [30.016701, 2.8886], 0),
    "N": gaussian([2.9924, 3.9986], [25.3766, 3.5004], 0),
    "O": gaussian([2.4485, 5.5589], [24.756199, 4.1372], 0),
    "S": gaussian([5.5480, 10.4241], [33.7108, 1.9034], 0),
  }

class two_gaussian_agarwal_1978:
  source = "Agarwal, R.C. (1978). Acta Cryst. A34, 791-809, Table 1."
  table = {
    "H": gaussian([0.4866, 0.5098], [34.1284, 8.8996], 0),
    "C": gaussian([3.0102, 2.9705], [29.9132, 2.8724], 0),
    "N": gaussian([3.0492, 3.9432], [25.0383, 3.4059], 0),
    "O": gaussian([3.2942, 4.6968], [20.0401, 3.1184], 0),
    "S": gaussian([5.6604, 10.3140], [33.0400, 1.8160], 0),
    "Fe3+": gaussian([10.3568, 12.6329], [8.1324, 0.8137], 0),
    "Fe2+": gaussian([11.6635, 12.3057], [9.0361, 0.5749], 0),
    "Zn2+": gaussian([5.7826, 22.2163], [11.7082, 1.8234], 0),
    "Ba2+": gaussian([12.1432, 41.8442], [21.7090, 1.4090], 0),
  }

class one_gaussian_agarwal_1978:
  source = "Agarwal, R.C. (1978). Acta Cryst. A34, 791-809, Table 3."
  table = {
    "C": gaussian([5.9074], [1.2913], 0),
    "N": gaussian([7.0411], [0.2065], 0),
    "O": gaussian([8.1561], [-0.8941], 0),
    "S": gaussian([15.8448], [-2.1392], 0),
  }
