from __future__ import absolute_import, division, print_function
from cctbx import miller
from cctbx import crystal
from cctbx import sgtbx
from cctbx.array_family import flex
from cmath import exp, pi
from scitbx.matrix import row, col
from libtbx.test_utils import approx_equal
from six.moves import zip

def exercise():
  ma = miller.array(
    miller.set(crystal.symmetry(unit_cell=(5,5,5, 90, 90, 90),
                                space_group=sgtbx.space_group('P 2x')),
               indices=flex.miller_index(
                 [(1,0,0), (0,1,0), (0,0,1),
                  (-1,0,0), (0,-1,0), (0,0,-1),
                  (1,1,0), (1,0,1), (0,1,1),
                  (-1,-1,0), (-1,0,-1), (0,-1,-1),
                  (1,-1,0), (1,0,-1), (0,1,-1),
                  (-1,1,0), (-1,0,1), (0,-1,1),
                  (1,1,1), (-1,1,1), (1,-1,1), (1,1,-1),
                  (-1,-1,-1), (1,-1,-1), (-1,1,-1), (-1,-1,1)])),
    data=flex.complex_double(flex.random_double(26), flex.random_double(26)))
  f_at_h = dict(zip(ma.indices(), ma.data()))
  for op in ("-x, y+1/2, -z", "x+1/2, -y, z-1/2"):
    op = sgtbx.rt_mx(op)
    original, transformed = ma.common_sets(
      ma.change_basis(sgtbx.change_of_basis_op(op.inverse())))
    for h, f in original:
      assert f == f_at_h[h]
    for h, op_f in transformed:
      assert approx_equal(
        op_f,
        f_at_h[h*op.r()]*exp(1j*2*pi*row(h).dot(col(op.t().as_double()))))

def run():
  exercise()
  print('OK')

if __name__ == '__main__':
  run()
