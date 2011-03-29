import cctbx.array_family.flex # import dependency
from cctbx import uctbx # import dependency
import boost.python
ext = boost.python.import_ext("cctbx_orientation_ext")
from cctbx_orientation_ext import *

class basis_type:
  direct = False
  reciprocal = True

class _(boost.python.injector,ext.crystal_orientation):

  def __getattr__(self,tag):
    mm = self.unit_cell().metrical_matrix()
    if tag=='A':
      return mm[0]
    elif tag=='B':
      return mm[1]
    elif tag=='C':
      return mm[2]
    elif tag=='D':
      return mm[5]
    elif tag=='E':
      return mm[4]
    elif tag=='F':
      return mm[3]
    mm = self.unit_cell().reciprocal().metrical_matrix()
    if tag=='As':
      return mm[0]
    elif tag=='Bs':
      return mm[1]
    elif tag=='Cs':
      return mm[2]
    elif tag=='Ds':
      return mm[5]
    elif tag=='Es':
      return mm[4]
    elif tag=='Fs':
      return mm[3]
    else:
      return

  def make_positive(self):
    from scitbx import matrix
    signed_volume = matrix.sqr(self.direct_matrix()).determinant()
    if signed_volume<0:
      return self.change_basis((-1.,0.,0.,0.,-1.,0.,0.,0.,-1.))
    return self

  def __copy__(self):
    return crystal_orientation(self.reciprocal_matrix(),basis_type.reciprocal)

  def __getinitargs__(self):
    return (self.reciprocal_matrix(),basis_type.reciprocal)

  def __str__(self):
    return "A-star:%10.6f%10.6f%10.6f%10.6f%10.6f%10.6f%10.6f%10.6f%10.6f\ncell:"%self.reciprocal_matrix()+\
       str(self.unit_cell())+"%.0f"%self.unit_cell().volume()

  def __eq__(self,other):
    S = self.reciprocal_matrix()
    O = other.reciprocal_matrix()
    return S[0]==O[0] and S[1]==O[1] and S[2]==O[2] and S[3]==O[3] and \
           S[4]==O[4] and S[5]==O[5] and S[6]==O[6] and S[7]==O[7] and S[8]==O[8]
