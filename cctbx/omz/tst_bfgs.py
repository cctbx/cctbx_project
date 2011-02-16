from __future__ import division
from cctbx.omz import bfgs
import scitbx.lbfgs
from scitbx.array_family import flex
from scitbx.linalg import eigensystem
from scitbx import matrix
from libtbx.test_utils import approx_equal
from math import cos, sin

def exercise_two_loop_recursion(fgh):
  x0 = flex.double([3.0, -4.0])
  g0 = fgh.g(x0)
  h0 = flex.double([[1,0],[0,1]])
  memory = []
  hg0 = bfgs.hg_two_loop_recursion(memory, hk0=h0, gk=g0)
  assert approx_equal(hg0, h0.matrix_multiply(g0))
  #
  x1 = x0 - 1/3 * hg0
  g1 = fgh.g(x1)
  h1 = bfgs.h_update(hk=h0, sk=x1-x0, yk=g1-g0)
  memory.append(bfgs.memory_element(s=x1-x0, y=g1-g0))
  hg1 = bfgs.hg_two_loop_recursion(memory, hk0=h0, gk=g1)
  assert approx_equal(hg1, h1.matrix_multiply(g1))
  #
  x2 = x1 - 1/5 * hg1
  g2 = fgh.g(x2)
  h2 = bfgs.h_update(hk=h1, sk=x2-x1, yk=g2-g1)
  memory.append(bfgs.memory_element(s=x2-x1, y=g2-g1))
  hg2 = bfgs.hg_two_loop_recursion(memory, hk0=h0, gk=g2)
  assert approx_equal(hg2, h2.matrix_multiply(g2))
  #
  x3 = x2 - 3/8 * hg2
  g3 = fgh.g(x3)
  h3 = bfgs.h_update(hk=h2, sk=x3-x2, yk=g3-g2)
  memory.append(bfgs.memory_element(s=x3-x2, y=g3-g2))
  hg3 = bfgs.hg_two_loop_recursion(memory, hk0=h0, gk=g3)
  assert approx_equal(hg3, h3.matrix_multiply(g3))

class refinery(object):

  def __init__(O, fgh):
    O.fgh = fgh
    O.x = flex.double([3,-4])
    O.prev_x = O.x.deep_copy()
    O.memory = []
    scitbx.lbfgs.run(target_evaluator=O)

  def compute_functional_and_gradients(O):
    return O.fgh.f(O.x), O.fgh.g(O.x)

  def callback_after_step(O, minimizer):
    xk = O.prev_x
    xl = O.x
    gk = O.fgh.g(xk)
    gl = O.fgh.g(xl)
    def check(bk, hk):
      hl = bfgs.h_update(hk, xl-xk, gl-gk)
      bl = bfgs.b_update(bk, xl-xk, gl-gk)
      assert approx_equal(matrix.sqr(hl).inverse(), bl)
      es = eigensystem.real_symmetric(bl)
      assert es.values().all_gt(0)
      assert bfgs.h0_scaling(sk=xl-xk, yk=gl-gk) > 0
    #
    bk = flex.double([[1,0],[0,1]])
    hk = bk
    check(bk, hk)
    #
    bk = O.fgh.h(xk)
    es = eigensystem.real_symmetric(bk)
    if (es.values().all_gt(0)):
      hk = bk.deep_copy()
      hk.matrix_inversion_in_place()
      check(bk, hk)
    #
    h0 = flex.double([[0.9,0.1],[-0.2,0.8]])
    hg_tlr = bfgs.hg_two_loop_recursion(memory=O.memory, hk0=h0, gk=gl)
    h_upd = h0
    for m in O.memory:
      h_upd = bfgs.h_update(hk=h_upd, sk=m.s, yk=m.y)
    hg_upd = h_upd.matrix_multiply(gl)
    assert approx_equal(hg_tlr, hg_upd)
    #
    O.memory.append(bfgs.memory_element(s=xl-xk, y=gl-gk))
    O.prev_x = O.x.deep_copy()

class fgh1(object):

  def f(O, x):
    assert x.size() == 2
    a,b = x
    return (a-1)**2 + (b-2)**2

  def g(O, x):
    a,b = x
    return flex.double([2*(a-1), 2*(b-2)])

  def h(O, x):
    return flex.double([
      [2, 0],
      [0, 2]])

class fgh2(object):

  def f(O, x):
    assert x.size() == 2
    a,b = x
    return (a-1)**2 + sin(b-2)

  def g(O, x):
    a,b = x
    return flex.double([2*(a-1), cos(b-2)])

  def h(O, x):
    a,b = x
    return flex.double([
      [2, 0],
      [0, -sin(b-2)]])

def run(args):
  assert len(args) == 0
  for fgh in [fgh1, fgh2]:
    exercise_two_loop_recursion(fgh=fgh())
    refinery(fgh=fgh())
  print "OK"

if (__name__ == "__main__"):
  import sys
  run(args=sys.argv[1:])
