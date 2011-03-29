from cctbx.array_family import flex # import dependency

import boost.python
boost.python.import_ext("rstbx_ext")
from rstbx_ext import *
import rstbx_ext as ext

import math

from cctbx.crystal_orientation import basis_type
from cctbx.crystal_orientation import ext as coext

class _(boost.python.injector, coext.crystal_orientation):

  def constrain(self,constraints):

    #algorithm 1.  Use pre-defined crystal_systems to give hard-coded restraints.
    # dps_core.constrainment.s_minimizer uses LBFGS minimizer to adapt
    # all 9 components of the orientation matrix.   This gives the best-fit
    # to the starting matrix (better than algorithm #2), but the disadvantage
    # is that it is keyed to the crystal_system descriptors.  It is therefore
    # not adapted to all small-molecule space groups (monoclinics),
    # and will not take into account non-standard settings.

    if constraints in ["triclinic","monoclinic",'orthorhombic','tetragonal',
                       "cubic","rhombohedral",'hexagonal']:

      from rstbx.dps_core.constrainment import s_minimizer
      S = s_minimizer(self,constraint=constraints)
      return S.newOrientation()

    #algorithm 2:  Tensor_rank_2 symmetrization
    # Advantages:  constraints are calculated directly from the space
    # group, so will account for non-standard settings.
    # Disadvantages:  drift away from starting orientation is greater than
    # for algorithm #1.

    from cctbx.sgtbx import space_group
    if isinstance(constraints,space_group):

      from labelit.symmetry.metricsym import a_g_conversion
      converter = a_g_conversion.AG()
      converter.forward(self)
      average_cell = constraints.average_unit_cell(self.unit_cell())
      converter.validate_and_setG( average_cell.reciprocal().metrical_matrix() )
      return Orientation(converter.back(),basis_type.reciprocal)

    #Future plans:  a hybrid approach.  Use algorithm 2 to do the
    # symmetrization, as it is clearly the best approach for 1) conciseness of
    # code, 2) supporting non-standard settings.  Then use an LBFGS
    # minimizer to minimize the psi, phi and theta offsets to the
    # original orientation.

class Orientation(coext.crystal_orientation):

  def __init__(self,either_matrix, basis_type_flag=basis_type.reciprocal):
    if isinstance(either_matrix,Orientation) or \
       isinstance(either_matrix,coext.crystal_orientation):
      coext.crystal_orientation.__init__(self,either_matrix)
    else:
      coext.crystal_orientation.__init__(self,either_matrix,basis_type_flag)

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

  def getOrientation(self):
    return Orientation(ext.dps_core.getOrientation(self))

  def niggli(self, cutoff = 25.):
    from rstbx.dps_core.zuoreduction import rwgk_niggli as support_niggli
    ni = support_niggli(self.getOrientation(),cutoff=cutoff) # orientation of reduced cell
    self.setOrientation(ni)
