from cctbx import crystal
import math

def dot3(a, b):
  return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]

def mod_positive3(x_frac):
  result = [math.fmod(x, 1.) for x in x_frac]
  for i in xrange(3):
    while (result[i] <  0.): result[i] += 1.
    while (result[i] >= 1.): result[i] -= 1.
  return tuple(result)

class plane_fractional:

  def __init__(self, n, p):
    self.n = n
    self.p = p

class planes_fractional:

  def __init__(self, space_group):
    self.space_group = space_group
    nc_dict = {}
    self.list = []
    for i_smx in xrange(space_group.order_p()):
      s = space_group(i_smx)
      r_info = s.r().info()
      if (r_info.type() < 2): continue
      if (r_info.sense() < 0): continue
      n = r_info.ev()
      p = s.t().mod_positive()
      assert p.den() == space_group.t_den()
      c = -dot3(n, p.num())
      if (not nc_dict.has_key((n,c))):
        nc_dict[(n,c)] = 0
        self.list.append(plane_fractional(n=n, p=p))

class plane_cartesian:

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
      p_c = self.unit_cell().orthogonalize(float(pf.p))
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
              if (result == None): result = pl.distance(xuc)
              else: result = min(result, pl.distance(xuc))
    return result
