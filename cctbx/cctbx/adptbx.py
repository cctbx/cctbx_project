from scitbx.python_utils.misc import import_regular_symbols
from cctbx_boost import adptbx_ext as ext
import_regular_symbols(globals(), ext.__dict__)
del import_regular_symbols

from cctbx.macro_mol import rotation_parameters
import random

def random_rotate_ellipsoid(u_cart):
  c = rotation_parameters.amore_alpha_beta_gamma_as_matrix(
    [random.uniform(0,360) for i in xrange(3)]).elems
  return c_u_c_transpose(c, u_cart)
