from __future__ import absolute_import, division, print_function
from scitbx import matrix
from libtbx.math_utils import iround
from cctbx import sgtbx
from cctbx.uctbx import unit_cell

'''Test the calculation of point-group symmetry operations by comparison
to published matrices.'''

cubic = unit_cell((10.,10.,10.,90.,90.,90.))
hexag = unit_cell((10.,10.,16.,90.,90.,120.))

ersatz_tests = """\n\n\n\n
2,1,0 0 1,-1 0 0 0 -1 0 0 0 1
2,1,0 1 -1,-1 0 0 0 0 -1 0 -1 0"""

IT96_Table_11_2 = """
Matrices for point-group symmetry operations and orientation of corresponding
symmetry elements, referred to a cubic, tetragonal, orthorhombic, monoclinic,
triclinic, or rhombohedral coordinate system, International Tables for
Crystallography, Volume A, 4th revised edition, 1996, p. 797.
1,1,1 1 1,1 0 0 0 1 0 0 0 1
3,1,1 1 1,0 0 1 1 0 0 0 1 0
3,-1,1 1 1,0 1 0 0 0 1 1 0 0
-1,1,0 0 0,-1 0 0 0 -1 0 0 0 -1
-3,1,1 1 1,0 0 -1 -1 0 0 0 -1 0
-3,-1,1 1 1,0 -1 0 0 0 -1 -1 0 0
2,1,0 0 1,-1 0 0 0 -1 0 0 0 1
3,1,1 -1 -1,0 0 -1 -1 0 0 0 1 0
3,-1,1 -1 -1,0 -1 0 0 0 1 -1 0 0
2,1,1 1 0,0 1 0 1 0 0 0 0 -1
2,1,1 -1 0,0 -1 0 -1 0 0 0 0 -1
4,1,0 0 1,0 -1 0 1 0 0 0 0 1
4,-1,0 0 1,0 1 0 -1 0 0 0 0 1
-2,1,0 0 1,1 0 0 0 1 0 0 0 -1
-3,1,1 -1 -1,0 0 1 1 0 0 0 -1 0
-3,-1,1 -1 -1,0 1 0 0 0 -1 1 0 0
-2,1,1 1 0,0 -1 0 -1 0 0 0 0 1
-2,1,1 -1 0,0 1 0 1 0 0 0 0 1
-4,1,0 0 1,0 1 0 -1 0 0 0 0 -1
-4,-1,0 0 1,0 -1 0 1 0 0 0 0 -1
2,1,0 1 0,-1 0 0 0 1 0 0 0 -1"""

IT96_Table_11_3 = """
Matrices for point-group symmetry operations and orientation of corresponding
symmetry elements, referred to a hexagonal coordinate system,
International Tables for Crystallography, Volume A, 4th revised edition,
1996, p. 798.
1,1,1 1 1,1 0 0 0 1 0 0 0 1
2,1,1 -1 0,0 -1 0 -1 0 0 0 0 -1
-2,1,1 1 0,0 -1 0 -1 0 0 0 0 1
6,1,0 0 1,1 -1 0 1 0 0 0 0 1
-3,1,0 0 1,0 1 0 -1 1 0 0 0 -1
-2,1,1 2 0,1 -1 0 0 -1 0 0 0 1
2,1,0 1 0,-1 0 0 -1 1 0 0 0 -1
-6,-1,0 0 1,0 -1 0 1 -1 0 0 0 -1"""

class CompareToInternationalTables:
  def __init__(self):
    self.coord_system_lookup = {
      'ersatz_tests':cubic,'IT96_Table_11_2':cubic,'IT96_Table_11_3':hexag}

  def test_list(self,key):
    niggli_cell=self.coord_system_lookup[key]
    frac = matrix.sqr(niggli_cell.fractionalization_matrix())
    orth = matrix.sqr(niggli_cell.orthogonalization_matrix())
    for item in globals()[key].split('\n')[5:]:
      type_plus_minus = int(item.split(',')[0])
      sense = int(item.split(',')[1])
      direct_space_axis = [int(x) for x in item.split(',')[2].split(' ')]
      reference_W = tuple([int(x) for x in item.split(',')[3].split(' ')])
      D = matrix.col(direct_space_axis)
      cartesian_axis = orth*D
      W_cart = matrix.sqr(
        sgtbx.n_fold_operator_from_axis_direction(
          ev_cart=cartesian_axis, n=abs(type_plus_minus), sense=sense))
      W_frac = frac*W_cart*orth
      W_as_int = matrix.sqr([iround(e) for e in W_frac.elems])
      if type_plus_minus<0: W_as_int = -1 * W_as_int
      assert W_as_int.elems == reference_W

if __name__=='__main__':
  C = CompareToInternationalTables()
  C.test_list('ersatz_tests')
  C.test_list('IT96_Table_11_2') #test non-hexagonal cells
  C.test_list('IT96_Table_11_3') #test hexagonal cells
  print("OK")
