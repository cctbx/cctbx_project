from scitbx.math import eigensystem
from scitbx import matrix
from stdlib import math
from scitbx.array_family import flex

def kearsley_rotation(reference_sites, other_sites):
  """
  Kearsley, S.K. (1989). Acta Cryst. A45, 208-210.
  On the orthogonal transformation used for structural comparison

  Added by Peter H. Zwart, Nov 3rd, 2006.
  """
  assert reference_sites.size() == other_sites.size()
  m = reference_sites - other_sites
  p = reference_sites + other_sites
  #extract the x y and z coords using dot products,
  #I did not see any other easy way of obtaining this info.
  dx=flex.vec3_double( [(1.0,0.0,0.0)]*reference_sites.size() )
  dy=flex.vec3_double( [(0.0,1.0,0.0)]*reference_sites.size() )
  dz=flex.vec3_double( [(0.0,0.0,1.0)]*reference_sites.size() )
  xm = m.dot( dx ); ym = m.dot( dy ); zm = m.dot( dz )
  xp = p.dot( dx ); yp = p.dot( dy ); zp = p.dot( dz )
  # make matrix elements
  a11 = flex.sum( xm*xm + ym*ym + zm*zm )
  a12 = flex.sum( yp*zm - ym*zp )
  a13 = flex.sum( xm*zp - xp*zm )
  a14 = flex.sum( xp*ym - xm*yp )
  a22 = flex.sum( yp*yp+zp*zp+xm*xm )
  a23 = flex.sum( xm*ym-xp*yp )
  a24 = flex.sum( xm*zm-xp*zp )
  a33 = flex.sum( xp*xp + zp*zp + ym*ym )
  a34 = flex.sum( ym*zm - yp*zp )
  a44 = flex.sum( xp*xp + yp*yp + zm*zm )
  #setup a 4*4 matrix
  grid = flex.grid(4,4)
  mtrx = flex.double( grid )
  #fil the elements
  mtrx[(0,0)] = a11; mtrx[(0,1)] = a12; mtrx[(0,2)] = a13; mtrx[(0,3)] = a14
  mtrx[(1,0)] = a12; mtrx[(1,1)] = a22; mtrx[(1,2)] = a23; mtrx[(1,3)] = a24
  mtrx[(2,0)] = a13; mtrx[(2,1)] = a23; mtrx[(2,2)] = a33; mtrx[(2,3)] = a34
  mtrx[(3,0)] = a14; mtrx[(3,1)] = a24; mtrx[(3,2)] = a34; mtrx[(3,3)] = a44
  # get the eigenvectors of this matrix please
  eigs = eigensystem.real_symmetric( mtrx )
  # Note that rms = sqrt(eigs.values()[3]/n_sites)!
  #
  # The values of the eigenvector with the smallest eigenvalue
  # form a unit quaternion
  q1 = eigs.vectors()[12]
  q2 = eigs.vectors()[13]
  q3 = eigs.vectors()[14]
  q4 = eigs.vectors()[15]
  # make a rotation matrix out of it for our convenience.
  a11 = q1*q1 + q2*q2 - q3*q3 - q4*q4
  a22 = q1*q1 + q3*q3 - q2*q2 - q4*q4
  a33 = q1*q1 + q4*q4 - q2*q2 - q3*q3
  a12 = 2.0*(q2*q3+q1*q4); a21 = 2.0*(q2*q3-q1*q4)
  a13 = 2.0*(q2*q4-q1*q3); a31 = 2.0*(q2*q4+q1*q3)
  a23 = 2.0*(q3*q4+q1*q2); a32 = 2.0*(q3*q4-q1*q2)
  r = matrix.sqr( [a11,a12,a13,a21,a22,a23,a31,a32,a33] )
  # done
  return r

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
