import cctbx.array_family.flex # for tuple mappings

import boost.python
ext = boost.python.import_ext("cctbx_adptbx_ext")
from cctbx_adptbx_ext import *

from scitbx.math import euler_angles_as_matrix
import random

def random_rotate_ellipsoid(u_cart):
  c = euler_angles_as_matrix(
    [random.uniform(0,360) for i in xrange(3)]).elems
  return c_u_c_transpose(c, u_cart)
