from cctbx.uctbx import unit_cell,fast_minimum_reduction
from rstbx.dps_core import Orientation as orient3d
from rstbx.dps_core.cell_assessment import unit_cell_too_small
from scitbx import matrix

def rwgk_niggli(UC,epsilon=None,cutoff=100.):

  '''Reference:
  R.W. Grosse-Kunstleve, N.K. Sauter and P.D. Adams.
  Numerically stable algorithms for the computation of reduced unit cells.
  Acta Cryst. A60, 1-6 (2004)
  '''

  if isinstance(UC,unit_cell):
    #return UC.niggli_cell(epsilon)
    return fast_minimum_reduction(UC).as_unit_cell()
  elif isinstance(UC,orient3d):
    uc = UC.unit_cell()
    unit_cell_too_small(uc,cutoff=cutoff)
    #R = uc.niggli_reduction(epsilon)
    R = fast_minimum_reduction(uc)
    #minimum reduction gets r_inv as a tuple instead of scitbx.matrix
    rinverse = matrix.sqr( R.r_inv() )
    #NIG = R.as_unit_cell().parameters()
    #MIN = fast_minimum_reduction(uc).as_unit_cell().parameters()
    UC.change_basis(rinverse.transpose().inverse().elems)
    return UC

if __name__=='__main__':
  uc = unit_cell((10,20,30,90,90,90))
  print uc.parameters()
  print rwgk_niggli(uc).parameters()
