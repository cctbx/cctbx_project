from scitbx.math import eigensystem
from scitbx.math import superpose_kearsley_rotation
from scitbx import matrix
from stdlib import math
from scitbx.array_family import flex
from stdlib import math as smath

def kearsley_rotation(reference_sites, other_sites):
  """
  Kearsley, S.K. (1989). Acta Cryst. A45, 208-210.
  On the orthogonal transformation used for structural comparison

  Added by Peter H. Zwart, Nov 3rd, 2006.
  Converted to C++ by Gabor Bunkoczi, Apr 2008.
  """
  return matrix.sqr(superpose_kearsley_rotation(
    reference_sites=reference_sites,
    other_sites=other_sites))

def kabsch_rotation(reference_sites, other_sites):
  """
Kabsch, W. (1976). Acta Cryst. A32, 922-923.
A solution for the best rotation to relate two sets of vectors

Based on a prototype by Erik McKee and Reetal K. Pai.

This implementation does not handle degenerate situations correctly
(e.g. if all atoms are on a line or plane) and should therefore not
be used in applications. It is retained here for development purposes
only.
  """
  assert reference_sites.size() == other_sites.size()
  sts = matrix.sqr(other_sites.transpose_multiply(reference_sites))
  eigs = eigensystem.real_symmetric((sts * sts.transpose()).as_sym_mat3())
  vals = list(eigs.values())
  vecs = list(eigs.vectors())
  a3 = list(matrix.col(vecs[:3]).cross(matrix.col(vecs[3:6])))
  a = matrix.sqr(list(vecs[:6])+a3)
  b = list(a * sts)
  for i in xrange(3):
    d = math.sqrt(math.fabs(vals[i]))
    if (d > 0):
      for j in xrange(3):
        b[i*3+j] /= d
  b3 = list(matrix.col(b[:3]).cross(matrix.col(b[3:6])))
  b = matrix.sqr(b[:6]+b3)
  return b.transpose() * a

class least_squares_fit(object):

  def __init__(self, reference_sites, other_sites, method="kearsley"):
    assert method in [None, "kearsley", "kabsch"]
    if (method is None): method = "kearsley"
    self.reference_sites = reference_sites
    self.other_sites = other_sites
    self.reference_shift = reference_sites.mean()
    self.other_shift = other_sites.mean()
    self.r = None
    if method == "kearsley":
      self.r = kearsley_rotation(
        reference_sites-self.reference_shift,
        other_sites-self.other_shift)
    if method == "kabsch":
      self.r = kabsch_rotation(
        reference_sites-self.reference_shift,
        other_sites-self.other_shift)
    self.t = matrix.col(self.reference_shift) \
           - self.r * matrix.col(self.other_shift)

  def other_sites_best_fit(self, additional_sites=None):
    if additional_sites:
      return self.r.elems * self.other_sites.concatenate(additional_sites) \
             + self.t.elems
    else:
      return self.r.elems * self.other_sites + self.t.elems

  def rt(self):
    return matrix.rt(tuple_r_t=(self.r, self.t))



"""
The NSD engine is a simple implementation of the normalized spatial discrepancy
as discussed by Kochin & Svergun J. Appl. Cryst 2001, 34, 33-41.

It can be used as a target for structure super positioning when no point-point correspondence is known.
"""

class nsd_engine(object):
  def __init__(self, fixed, d_fixed=None, d_moving=None):
    self.fixed = fixed
    self.d_fixed = d_fixed
    self.d_moving = d_moving
    if self.d_fixed is None:
      self.d_fixed = self.get_mean_distance(self.fixed)


  def get_mean_distance(self,xyz):
    N = xyz.size()
    d = 0
    count=0
    for ii in range(N):
      for jj in range(ii+1,N):
        dd = flex.double(xyz[ii])-flex.double(xyz[jj])
        d += dd.norm()
        count += 1
    d = d / count
    return d*d


  def nsd(self,moving,d_moving=None):
    if d_moving is None:
      self.d_moving = self.get_mean_distance(moving)
    else:
      self.d_moving = d_moving

    # loop over all sites in fixed, find the minimum for each site
    tot_rho_mf = 0
    tot_rho_fm = 0
    for site in moving:
      dd = self.fixed-site
      dd = flex.min( dd.norms() )
      tot_rho_mf+=dd*dd

    for site in self.fixed:
      dd = moving-site
      dd = flex.min( dd.norms() )
      tot_rho_fm+=dd
    tot_rho_fm = tot_rho_fm / (self.fixed.size()*self.d_fixed )
    tot_rho_mf = tot_rho_mf / (moving.size()*self.d_moving )
    result = smath.sqrt((tot_rho_fm+tot_rho_mf)/2.0)
    return result



def tst_nsd():
  moving1 = flex.vec3_double()
  moving2 = flex.vec3_double()
  fixed  = flex.vec3_double()
  for ii in range(20):
    noise = flex.random_double(3)*2-1.0
    xyz = flex.random_double(3)*5
    fixed.append( list(xyz) )
    moving1.append(  list(xyz + noise/10) )
    moving2.append(  list(xyz + noise/2) )

  ne = nsd_engine(fixed)
  a = ne.nsd(fixed)
  b = ne.nsd(moving1)
  c = ne.nsd(moving2)
  assert abs(a)<1e-6
  assert(b<=c)



if __name__ == "__main__":
  tst_nsd()
  print "OK"
