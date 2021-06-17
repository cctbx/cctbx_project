from __future__ import absolute_import, division, print_function
from cctbx import crystal
from cctbx import sgtbx
from boost_adaptbx.boost import rational
import math
import sys
from six.moves import range

def dot3(a, b):
  return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]

def mod_positive3(x_frac):
  result = [math.fmod(x, 1.) for x in x_frac]
  for i in range(3):
    while (result[i] <  0.): result[i] += 1.
    while (result[i] >= 1.): result[i] -= 1.
  return tuple(result)

class plane_fractional(object):

  def __init__(self, s, n, p, c):
    self.s = s
    self.n = n
    self.p = p
    self.c = c

  def operation(self):
    return sgtbx.rt_mx(self.s.r().minus_unit_mx(), self.s.t())

  def algebraic(self):
    op = self.operation()
    # eliminate minus signs to make the expression look more like
    # what people are used to see
    r = op.r().num()
    d = op.r().den()
    signs = [None,None,None]
    for j in range(3):
      for i in range(3):
        rij = r[i*3+j]
        if (signs[i] is None and rij != 0):
          if (rij < 0): signs[i] = -1
          else: signs[i] = 1
    m = list(sgtbx.rot_mx(d, d).num())
    for i in range(3):
      if (signs[i] == -1): m[i*4] *= -1
    return str(sgtbx.rt_mx(op.r().multiply(sgtbx.rot_mx(m,d)), op.t()))

class planes_fractional(object):

  def __init__(self, space_group):
    self.space_group = space_group
    nc_dict = {}
    self.list = []
    for i_smx in range(space_group.order_p()):
      s = space_group(i_smx)
      r_info = s.r().info()
      if (r_info.type() < 2): continue
      if (r_info.sense() < 0): continue
      n = r_info.ev()
      p = s.t().mod_positive()
      assert p.den() == space_group.t_den()
      c = -dot3(n, p.num())
      if ((n,c) not in nc_dict):
        nc_dict[(n,c)] = 0
        self.list.append(
          plane_fractional(s=s, n=n, p=p, c=rational.int(c,p.den())))

class plane_cartesian(object):

  def __init__(self, n, c):
    self.n = n
    self.c = c

  def distance(self, x):
    return abs(dot3(self.n, x) + self.c)

class planes_cartesian(crystal.symmetry):

  def __init__(self, crystal_symmetry):
    crystal.symmetry._copy_constructor(self, crystal_symmetry)
    planes_f = planes_fractional(self.space_group())
    self.list = []
    for pf in planes_f.list:
      n_c = self.unit_cell().orthogonalize(pf.n)
      p_c = self.unit_cell().orthogonalize(pf.p.as_double())
      len_n_c = math.sqrt(dot3(n_c, n_c))
      n_c = tuple([x/len_n_c for x in n_c])
      pl = plane_cartesian(n_c, -dot3(n_c, p_c))
      assert pl.distance(p_c) < 1.e-6
      self.list.append(pl)
    self.patterson_group = self.space_group().build_derived_patterson_group()

  def count(self):
    return len(self.list)

  def min_distance(self, x_frac):
    result = None
    for s in self.patterson_group:
      x = mod_positive3(s * x_frac)
      for u0 in (-1,0,1):
        for u1 in (-1,0,1):
          for u2 in (-1,0,1):
            xu = [x[0]+u0, x[1]+u1, x[2]+u2]
            xuc = self.unit_cell().orthogonalize(xu)
            for pl in self.list:
              if (result is None): result = pl.distance(xuc)
              else: result = min(result, pl.distance(xuc))
    return result

def test_call_back(flags, space_group_info):
  planes = planes_fractional(space_group_info.group())
  for plane in planes.list:
    print(plane.s, plane.algebraic(), plane.n, plane.c, plane.p)

def test():
  from cctbx.development import debug_utils
  debug_utils.parse_options_loop_space_groups(sys.argv[1:], test_call_back)

if (__name__ == "__main__"):
  test()
