from __future__ import absolute_import, division, print_function
from scitbx.math import tetrahedron
from libtbx.test_utils import approx_equal
from itertools import permutations
from random import random
from six.moves import range

def exercise_volume():
  for (a,b,c) in ((1,2,3), (2,3,1), (3,1,2)):
    for vertices in permutations(((1,1,1), (1+a,1,1), (1,1+b,1), (1,1,1+c))):
      t = tetrahedron(vertices)
      assert t.vertices == vertices
      vol = a*b*c/6
      assert approx_equal(t.volume(), vol)
  # regular tetrahedron
  for vertices in permutations(((1,1,1), (1,-1,-1), (-1,1,-1), (-1,-1,1))):
    t = tetrahedron(vertices)
    assert approx_equal(t.volume(), 8/3)

def exercise_gradients(n_points):
  eta = 1e-6
  for i in range(n_points):
    v = tuple((random(), random(), random()) for j in range(4))
    t = tetrahedron(v)
    g = t.gradients()
    for j in range(4):
      for k in range(3):
        (tm, tp) = tuple(tetrahedron(
          v[:j] + (v[j][:k] + (v[j][k] + h,) + v[j][k+1:],) + v[j+1:])
          for h in (-eta, eta))
        dv = (tp.volume() - tm.volume())/(2*eta)
        assert approx_equal(dv, g[j][k]), (j,k)

def run():
  exercise_volume()
  exercise_gradients(n_points=100)
  print('OK')

if __name__ == '__main__':
  run()
