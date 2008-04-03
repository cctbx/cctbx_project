from scitbx.math import eigensystem
from scitbx.math import superpose_kearsley_rotation
from scitbx import matrix
from stdlib import math
from scitbx.array_family import flex

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
