from cctbx import sgtbx
from cctbx import matrix
from scitbx.python_utils.math_utils import iround
from scitbx.python_utils.misc import store
from scitbx.array_family import flex
import math

# Y. Le Page
# The derivation of the axes of the conventional unit cell from the
# dimensions of the Buerger-reduced cell
# J. Appl. Cryst. (1982). 15, 255-259

def two_fold_matrix_from_axis_direction(ev_cart):
  f = 2. / ev_cart.norm()
  x,y,z = ev_cart.elems
  return matrix.sqr((f*x*x-1., f*x*y,    f*x*z,
                     f*y*x,    f*y*y-1., f*y*z,
                     f*z*x,    f*z*y,    f*z*z-1.))

def as_integer_matrix(m):
  result = []
  for e in m.elems:
    result.append(iround(e))
  return result

class group_search:

  def __init__(self, modulus=2):
    self.modulus = modulus
    self._potential_axes = None

  def compute_potential_axes(self):
    result = []
    for u in xrange(0,self.modulus+1):
      for v in xrange(-self.modulus,self.modulus+1):
        for w in xrange(-self.modulus,self.modulus+1):
          for h in xrange(0,self.modulus+1):
            for k in xrange(-self.modulus,self.modulus+1):
              for l in xrange(-self.modulus,self.modulus+1):
                abs_uh = abs(u*h+v*k+w*l)
                if (abs_uh in (1,2)):
                  result.append(store(
                    u=matrix.col((u,v,w)),
                    h=matrix.row((h,k,l)),
                    abs_uh=abs_uh))
    return result

  def potential_axes(self):
    if (self._potential_axes == None):
      self._potential_axes = self.compute_potential_axes()
    return self._potential_axes

  def __call__(self, niggli_cell, max_delta=3.):
    max_delta *= math.pi/180
    frac = matrix.sqr(niggli_cell.fractionalization_matrix())
    orth = matrix.sqr(niggli_cell.orthogonalization_matrix())
    deltas = flex.double()
    ts = []
    for axis in self.potential_axes():
      t = orth * axis.u
      tau = matrix.col((axis.h * frac).elems)
      delta = abs(math.atan2(abs(t.cross(tau)), axis.abs_uh))
      if (delta < max_delta):
        deltas.append(delta)
        ts.append(t)
    perm = flex.sort_permutation(deltas)
    group = sgtbx.space_group()
    for i in perm:
      w_cart = two_fold_matrix_from_axis_direction(ts[i])
      w_frac = as_integer_matrix(frac*w_cart*orth)
      s = sgtbx.rt_mx(w_frac, (0,0,0))
      expanded_group = sgtbx.space_group(group)
      try:
        expanded_group.expand_smx(s)
      except RuntimeError, e:
        ee = "cctbx Error: Non-crystallographic rotation matrix encountered."
        if (str(e) != ee): raise
        break
      group = expanded_group
    return group

group = group_search()
