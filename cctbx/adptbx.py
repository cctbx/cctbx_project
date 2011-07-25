from cctbx.array_family import flex # for tuple mappings

import boost.python
ext = boost.python.import_ext("cctbx_adptbx_ext")
from cctbx_adptbx_ext import *

import scitbx.math
import random
import math

u_as_b_factor = u_as_b(1)
b_as_u_factor = b_as_u(1)

mtps = -2 * math.pi**2
mtpss = mtps**2

def random_rotate_ellipsoid(u_cart, r_min = 0, r_max = 360):
  c = scitbx.math.euler_angles_as_matrix(
    [random.uniform(r_min,r_max) for i in xrange(3)], deg=True).elems
  return c_u_c_transpose(c, u_cart)

def random_u_cart(u_scale=1, u_min=0):
  return random_rotate_ellipsoid(u_cart=[random.random()*u_scale+u_min
    for i in xrange(3)] + [0,0,0])

def debye_waller_factor_u_star_gradients(h, u_star):
  return flex.double(debye_waller_factor_u_star_gradient_coefficients(h)) \
       * (mtps * debye_waller_factor_u_star(h, u_star))

def debye_waller_factor_u_star_curvatures(h, u_star):
  return debye_waller_factor_u_star_curvature_coefficients(h=h) \
       * (mtpss * debye_waller_factor_u_star(h, u_star))
