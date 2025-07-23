from __future__ import absolute_import, division, print_function
from scitbx.linalg import eigensystem
from scitbx.math import superpose_kearsley_rotation
from scitbx import matrix
from scitbx.stdlib import math, random
from scitbx.array_family import flex
from scitbx import differential_evolution as de
from scitbx.math import euler_angles as euler
from six.moves import range

def kearsley_rotation(reference_sites, other_sites):
  """
  Kearsley, S.K. (1989). Acta Cryst. A45, 208-210.
  On the orthogonal transformation used for structural comparison

  Added by Peter H. Zwart, Nov 3rd, 2006.
  Converted to C++ by Gabor Bunkoczi, Apr 2008.
  """
  assert reference_sites.size() == other_sites.size(), "%d != %d" % (
      reference_sites.size(), other_sites.size())
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
  for i in range(3):
    d = math.sqrt(math.fabs(vals[i]))
    if (d > 0):
      for j in range(3):
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

# TODO test
def sieve_fit(sites_fixed,
               sites_moving,
               selection=None,
               frac_discard=0.5):
  """
  Reference: Chothia & Lesk???
  """
  assert (sites_fixed.size() == sites_moving.size() > 0)
  if (selection is None):
    selection = flex.bool(sites_fixed.size(), True)
  # step 1: superpose using originally selected atoms
  sites_fixed_aln = sites_fixed.select(selection)
  sites_moving_aln = sites_moving.select(selection)
  lsq_fit_obj = least_squares_fit(
    reference_sites=sites_fixed_aln,
    other_sites=sites_moving_aln)
  sites_moving_new = lsq_fit_obj.other_sites_best_fit()
  # step 2: discard 50% of sites that deviate the most, and superpose again
  deltas = (sites_fixed_aln - sites_moving_new).norms()
  deltas_sorted = flex.sorted(deltas)
  cutoff = deltas_sorted[int((1-frac_discard)*deltas.size())]
  selection = (deltas > cutoff)
  if (selection.count(True) == 0):
    return lsq_fit_obj
  sites_fixed_aln = sites_fixed_aln.select(selection)
  sites_moving_aln = sites_moving_aln.select(selection)
  lsq_fit_obj = least_squares_fit(
    reference_sites=sites_fixed_aln,
    other_sites=sites_moving_aln)
  return lsq_fit_obj

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
    if self.d_moving is None:
      self.d_moving = self.get_mean_distance(moving)
    if d_moving is not None:
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
    result = math.sqrt((tot_rho_fm+tot_rho_mf)/2.0)
    return result


class nsd_rigid_body_fitter(object):
  def __init__(self, fixed, moving):
    self.fixed = fixed
    self.moving = moving

    self.nsde = nsd_engine(self.fixed)

    self.d_fixed  = math.sqrt(self.nsde.d_fixed)
    self.d_moving = math.sqrt(self.nsde.get_mean_distance( self.moving ) )

    self.m_com = self.moving.mean()
    self.f_com = self.fixed.mean()
    self.n_mov = self.moving-self.m_com

    self.d = (self.d_fixed+self.d_moving)/12

    self.n = 6
    pi = math.pi
    self.domain = [ (-pi,pi),(-pi,pi), (-pi,pi), (-self.d,self.d),(-self.d,self.d), (-self.d,self.d) ]
    self.x = None
    self.optimizer = de.differential_evolution_optimizer(self, population_size=12,f=0.85,cr=0.95, n_cross=2,eps=1e-2,
                         show_progress=False,show_progress_nth_cycle=20)




  def move_points(self,vector):
    shift = (vector[0], vector[1], vector[2])
    abg   =  ( vector[3], vector[4], vector[5] )
    matrix = euler.zyz_matrix( vector[3], vector[4], vector[5] )
    new_xyz = matrix*self.n_mov+self.f_com+shift
    return new_xyz


  def target(self,vector):
    nxyz = self.move_points(vector)
    result = self.nsde.nsd(nxyz)
    return result

  def print_status(self, min_target, mean_target, best_vector, txt):
    print(min_target,mean_target,list(best_vector), txt)

  def best_shifted(self):
    nxyz = self.move_points(self.x)
    return nxyz


def tst_nsd():
  moving1 = flex.vec3_double()
  moving2 = flex.vec3_double()
  fixed  = flex.vec3_double()
  max_noise = 0
  for ii in range(10):
    noise = flex.random_double(3)*2-1.0
    if noise.norm() > max_noise:
      max_noise = noise.norm()
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

  matrix = euler.zyz_matrix(0.7,1.3,2.1)
  fixed_r = matrix*moving1+(8,18,28)
  fitter = nsd_rigid_body_fitter( fixed,fixed_r)
  nxyz = fitter.best_shifted()
  dd = nxyz[0:fixed.size()]-fixed
  dd = dd.norms()
  dd = flex.max(dd)
  assert (dd<2.00*max_noise/10)


if __name__ == "__main__":
  random.seed(0)
  flex.set_random_seed(0)
  for ii in range(10):
    tst_nsd()
  print("OK")
