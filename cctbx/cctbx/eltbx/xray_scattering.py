import cctbx.array_family.flex # for tuple mappings

from scitbx.python_utils import misc
ext = misc.import_ext("cctbx_boost.eltbx.xray_scattering_ext")
misc.import_regular_symbols(globals(), ext.__dict__)
del misc

def best_approximation(scattering_type):
  if (scattering_type == "const"):
    return gaussian(1)
  return wk1995(scattering_type, 1).fetch()

two_gaussian_agarwal_isaacs = {
  "source":
      "ccp4/lib/data/atomsf.lib Revision 1.4, Thu Feb 13 14:10:58 1997 UTC ",
  "H":
      gaussian((0.7932, 0.1949), (24.2157, 2.1089), 0),
  "C":
      gaussian((2.9972, 2.9791), (30.016701, 2.8886), 0),
  "N":
      gaussian((2.9924, 3.9986), (25.3766, 3.5004), 0),
  "O":
      gaussian((2.4485, 5.5589), (24.756199, 4.1372), 0),
  "S":
      gaussian((5.5480, 10.4241), (33.7108, 1.9034), 0),
}
