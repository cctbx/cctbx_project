from scitbx.array_family import flex
from scitbx import lstbx, sparse
from libtbx.test_utils import approx_equal

def exercise_basic_normal_equations():
  py_eqs = [ ( 1, (-1,  0,  0),  1),
             ( 2, ( 2, -1,  0),  3),
             (-1, ( 0,  2,  1),  2),
             (-2, ( 0,  1,  0), -2),
             ]

  eqs_0 = lstbx.normal_equations(3)
  for b, a, w in py_eqs:
    eqs_0.add_equation(right_hand_side=b,
                       design_matrix_row=flex.double(a),
                       weight=w)

  eqs_1 = lstbx.normal_equations(3)
  b = flex.double()
  w = flex.double()
  a = sparse.matrix(len(py_eqs), 3)
  for i, (b_, a_, w_) in enumerate(py_eqs):
    b.append(b_)
    w.append(w_)
    for j in xrange(3):
      if a_[j]: a[i, j] = a_[j]
  eqs_1.add_equations(right_hand_side=b, design_matrix=a, weights=w)

  assert approx_equal(
    eqs_0.normal_matrix_packed_u, eqs_1.normal_matrix_packed_u, eps=1e-15)
  assert approx_equal(
    eqs_0.right_hand_side, eqs_1.right_hand_side, eps=1e-15)
  assert approx_equal(
    list(eqs_0.normal_matrix_packed_u), [ 13, -6, 0, 9, 4, 2 ], eps=1e-15)
  assert approx_equal(
    list(eqs_0.right_hand_side), [ 11, -6, -2 ], eps=1e-15)

def exercise_normal_equations_separating_scale_factor():
  eqs = lstbx.normal_equations_separating_scale_factor(3)
  eqs.add_equation(y_calc=1,
                   grad_y_calc=flex.double((1, 2, 3)),
                   y_obs=1,
                   weight=1)
  eqs.finalise()
  a, b = eqs.reduced_equations()
  assert a.size() == 6
  assert list(b) == [0, 0, 0]

def run():
  exercise_basic_normal_equations()
  exercise_normal_equations_separating_scale_factor()
  print 'OK'

if __name__ == '__main__':
  run()
