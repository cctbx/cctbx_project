import cctbx.array_family.flex # for tuple mappings

from scitbx.python_utils import misc
ext = misc.import_ext("cctbx_boost.adptbx_ext")
misc.import_regular_symbols(globals(), ext.__dict__)
del misc

from cctbx.macro_mol import rotation_parameters
import random

def random_rotate_ellipsoid(u_cart):
  c = rotation_parameters.amore_alpha_beta_gamma_as_matrix(
    [random.uniform(0,360) for i in xrange(3)]).elems
  return c_u_c_transpose(c, u_cart)
