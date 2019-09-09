# -*- coding: utf-8 -*-

from __future__ import absolute_import, division, print_function

import itertools

import scitbx.matrix
from rstbx.dps_core.cell_assessment import unit_cell_too_small
from rstbx.indexing_api import cpp_absence_test

def _is_collinear(x,y): # X x Y cross product is zero
  return x[0]*y[1]-x[1]*y[0]==x[1]*y[2]-x[2]*y[1]==x[2]*y[0]-x[0]*y[2]==0

def _is_coplanar(x,y,z):
  #triple product is zero; (X x Y).Z
  x = scitbx.matrix.row(x)
  y = scitbx.matrix.row(y)
  z = scitbx.matrix.row(z)
  return x.cross(y).dot(z)==0

def _generate_reindex_transformations():
    '''This implementation is based on the algorithm described in §2.5
    steps 1-3 of Sauter et al. (2004). J. Appl. Cryst. 37, 399-409.
    https://doi.org/10.1107/S0021889804005874

    The reindex transformations are specific for a particular
    presence condition, such as H + 2K + 3L = 5n.  The transformation
    is applied in reciprocal space, and is intended to change the
    original incorrect basis set a*',b*',c*' into the correct basis
    set a*,b*,c*. The meaning of the correction matrix A is as follows:

           a* = A00(a*') + A01(b*') + A02(c*')
           b* = A10(a*') + A11(b*') + A12(c*')
           c* = A20(a*') + A21(b*') + A22(c*')

    The choice of A is not unique, we use an algorithm to select a
    particular solution.  Briefly, for the first row of A choose the row
    vector HKL which satisfies the presence condition, and is shortest in
    length.  For the second row choose the next shortest allowed row vector
    that is not collinear with the first.  The third allowed row vector is
    the next shortest not coplanar with the first two.  We check to see
    that the determinant is positive (or switch first two rows) and of
    magnitude equal to the mod factor; this assures that the unit cell will
    be reduced in volume by the appropriate factor.

    Our approach sometimes backfires:  an already too-small unit cell can
    produce a positive absence test; the cell will then be reduced in volume
    indefinitely.  Therefore the application always uses a cell volume filter
    after making the correction.
    '''
    # modularities 2,3,5 were sufficient for every two-image case
    # need up to 11 for Fig 4 in the single-image indexing
    modularities = [2,3,5]

    mod_range = range(max(modularities), -max(modularities)-1,-1)
    points = itertools.product(mod_range, mod_range, mod_range)
    # sort by increasing distance and descending size
    spiral_order = list(sorted(points, key=lambda v: (sum(c*c for c in v), -sum(v))))
    spiral_order.remove((0,0,0))  # G0 in the paper (step 1)

    representatives = []  # G1 in the paper (step 2)
    # The vector representations connote systematic absence conditions.
    # For example, the vector v = (1,2,3) means H + 2K + 3L = ?n,
    # where the ? represents the modularity (2,3,5,...) specified elsewhere
    for vector in spiral_order:
      if sum(c*c for c in vector) > 6: break
      if any(_is_collinear(vector, item) for item in representatives): continue
      representatives.append(vector)

    # Now generate the matrices for every reflection condition (step 3)
    reindex = []
    for vec in representatives:
      for mod in modularities:
        candidate_points = (pt for pt in spiral_order if sum(v*p for v,p in zip(vec,pt))%mod == 0)
        # find three points that are not coplanar
        first = next(candidate_points)
        while True:
          second = next(candidate_points)
          if not _is_collinear(first, second):
            break
        while True:
          third = next(candidate_points)
          if not _is_coplanar(first, second, third):
            break
        A = scitbx.matrix.sqr(first + second + third)
        if A.determinant() < 0:
          A = scitbx.matrix.sqr(second + first + third)
        assert A.determinant() == mod
        reindex.append({'mod':mod, 'vec':vec, 'trans':A})
    return reindex

R = _generate_reindex_transformations()

class AbsenceHandler:
  def __init__(self):
    self.recursion_limit=8

  def absence_detected(self,hkllist):
    self.hkl = hkllist
    self.N   = self.hkl.size()
    self.flag = None

    for test in R:
        cum = cpp_absence_test(self.hkl,test['mod'],test['vec'])
        for counter in range(test['mod']):
          #print test['vec'],test['mod'],float(cum[counter])/self.N
          if float(cum[counter])/self.N > 0.8 and counter==0:
            # (if counter != 0 there is no obvious way to correct this)
            #print "Detected exclusive presence of %dH %dK %dL = %dn, remainder %d"%(
            #         test['vec'][0],test['vec'][1],test['vec'][2],test['mod'],counter)
            self.flag = {'vec':test['vec'],'mod':test['mod'],
                         'remainder':counter, 'trans':test['trans'].elems}
            return 1
    return 0

  def correct(self,orientation):
    if self.flag is None:
      raise RuntimeError("no correction necessary")
    M1 = scitbx.matrix.sqr(self.flag['trans'])
    corrected = orientation.change_basis(M1.transpose().elems)
    unit_cell_too_small(corrected.unit_cell(),cutoff = 100.)
    return corrected

  def list(self,hkllist):
      self.hkl = hkllist
      count = 0
      for m in self.hkl:
        print(count,m)
        print("                ",(m[0]%2,m[1]%2,m[2]%2), end=' ')
        print((m[0]%3,m[1]%3,m[2]%3), end=' ')
        print(((m[1]-m[2])%2,(m[1]+m[2])%2), end=' ')
        print(((m[2]-m[0])%2,(m[2]+m[0])%2), end=' ')
        print(((m[0]-m[1])%2,(m[0]+m[1])%2))
        count+=1

if __name__=='__main__':
  def pelem(arg):
    return arg.elems.__repr__()
  scitbx.matrix.sqr.__repr__ = pelem
  import pprint
  pprint.pprint( R)
  print(len(R))
