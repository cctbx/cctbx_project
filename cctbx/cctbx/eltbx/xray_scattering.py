import cctbx.array_family.flex # for tuple mappings

from scitbx.python_utils import misc
ext = misc.import_ext("cctbx_boost.eltbx.xray_scattering_ext")
misc.import_regular_symbols(globals(), ext.__dict__)
del misc

def best_approximation(scattering_type):
  if (scattering_type == "const"):
    return gaussian(1)
  return wk1995(scattering_type, 1).fetch()
