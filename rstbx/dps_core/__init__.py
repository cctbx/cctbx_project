from cctbx.array_family import flex

import boost.python
boost.python.import_ext("rstbx_ext")
from rstbx_ext import *
import rstbx_ext as ext

import types,math
from rstbx.dps_core.constrainment import s_minimizer

from cctbx.crystal_orientation import basis_type
from cctbx.crystal_orientation import ext as coext

class Orientation(coext.crystal_orientation):

  def __init__(self,either_matrix, basis_type_flag=basis_type.reciprocal):
    if isinstance(either_matrix,Orientation) or \
       isinstance(either_matrix,coext.crystal_orientation):
      coext.crystal_orientation.__init__(self,either_matrix)
    else:
      coext.crystal_orientation.__init__(self,either_matrix,basis_type_flag)

  def constrain(self,crystal_system):
    S = s_minimizer(self,constraint=crystal_system)
    newOrient = S.newOrientation()
    return newOrient

def combocmp(a,b):
  #gives -1,0,1 depending on closeness of combo to (0,0,0)
  a_measure = a[0]+a[1]+a[2]
  b_measure = b[0]+b[1]+b[2]
  if a_measure<b_measure: return -1
  if a_measure==b_measure: return 0
  return 1

def directional_show(direction,message):
  print message,"%.4f %8.2f %8.2f kmax=%2d kval=%5.1f kval2=%5.1f kval3=%5.1f"%(
    direction.real,180*direction.psi/math.pi, 180.*direction.phi/math.pi,
    direction.kmax, direction.kval,direction.kval2,direction.kval3)

class dps_core(ext.dps_core):
  def __init__(self):
    ext.dps_core.__init__(self)

  def combos(self,basis=10):
    """All interesting combinations of the directional candidates.
    Parameter MAXINDEX is the maximum id number considered when choosing
    combos.  With combos sorted before return, as below, I feel more
    comfortable increasing this parameter if needed to index a difficult case"""
    nc = self.n_candidates()
    MAXINDEX=basis
    bases = range(min(MAXINDEX,nc))
    comb=[]
    for pp in xrange(len(bases)-2):
      for qq in xrange(pp+1,len(bases)-1):
        for rr in xrange(qq+1, len(bases)):
          comb.append((bases[pp],bases[qq],bases[rr]))
    comb.sort(combocmp)
    return comb

  def setA(self,argument):
         from scitbx import matrix as vector # to clarify role of column vector
         assert type(argument) == types.ListType
         assert not 0 in [isinstance(x,Direction) for x in argument]
         self.combo_state = argument
         #case of a list of dptbx.Directions.
         realaxis=[]
         for i in xrange(3):
           realaxis.append(  vector.col(argument[i].dvec) * argument[i].real )
         matA = [  realaxis[0].elems[0],realaxis[0].elems[1],realaxis[0].elems[2],
                   realaxis[1].elems[0],realaxis[1].elems[1],realaxis[1].elems[2],
                   realaxis[2].elems[0],realaxis[2].elems[1],realaxis[2].elems[2]  ]
         ext.dps_core.set_orientation_direct_matrix(self,matA)

  def getOrientation(self):
    return Orientation(ext.dps_core.getOrientation(self))

  def niggli(self, cutoff = 25.):
    from rstbx.dps_core.zuoreduction import rwgk_niggli as support_niggli
    ni = support_niggli(self.getOrientation(),cutoff=cutoff) # orientation of reduced cell
    self.setOrientation(ni)

