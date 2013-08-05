# -*- coding: utf-8; py-indent-offset: 2 -*-
"""
Methods for examining bond angles around an ion and categorizing
the coordination geometry.
"""
from __future__ import division
from collections import OrderedDict
from math import sqrt
from scitbx.matrix import col, rec

def _tetrahedron():
  return rec(
    [1, 1, 1,
     -1, -1, 1,
     -1, 1, -1,
     1, -1, -1],
     n = (4, 3)
    )

def _octahedron():
  return rec(
    [1, 0, 0,
     -1, 0, 0,
     0, 1, 0,
     0, -1, 0,
     0, 0, 1,
     0, 0, -1],
     n = (6, 3)
    )

def _trigonal_plane():
  return rec(
    [0, 1, 0,
     -sqrt(3) / 2, -1 / 2, 0,
     sqrt(3) / 2, -1 / 2, 0],
     n = (3, 3)
    )

def _square_plane():
  return rec(
    [0, 1, 0,
     0, -1, 0,
     1, 0, 0,
     -1, 0, 0],
     n = (4, 3)
    )

def _concatenate(*args):
  lst = [i for arg in args for i in arg]
  return rec(lst, n = (len(lst) // 3, 3))

def _pyramid(base):
  return _concatenate(base, col([0, 0, 1]))

def _bipyramid(base):
  return _concatenate(base, col([0, 0, 1]), col([0, 0, -1]))

def _square_pyramid():
  return _pyramid(_square_plane())

def _square_bipyramid():
  return _bipyramid(_square_plane())

def _trigonal_pyramid():
  return _pyramid(_trigonal_plane())

def _trigonal_bipyramid():
  return _bipyramid(_trigonal_plane())

def _trigonal_prism():
  return _concatenate(
    _trigonal_plane() + rec([0, 0, 1] * 3, n = (3, 3)),
    _trigonal_plane() + rec([0, 0, -1] * 3, n = (3, 3)),
    )

def _pentagon():
  # Taken from http://mathworld.wolfram.com/Pentagon.html
  c_1 = (sqrt(5) - 1) / 4
  c_2 = (sqrt(5) - 1) / 4
  s_1 = sqrt(10 + 2 * sqrt(5)) / 4
  s_2 = sqrt(10 - 2 * sqrt(5)) / 4

  return rec(
    [1, 0, 0,
     c_1, s_1, 0,
     -c_2, s_2, 0,
     -c_2, -s_2, 0,
     c_1, -s_1, 0],
     n = (5, 3)
    )

def _pentagonal_pyramid():
  return _pyramid(_pentagon())

def _pentagonal_bipyramid():
  return _bipyramid(_pentagon())

def _see_saw():
  """
  An octahedron missing two adjacent points.
  """
  return rec(
    [1, 0, 0,
     -1, 0, 0,
     0, 1, 0,
     0, 0, 1,],
     n = (4, 3)
    )

def _three_legs():
  """
  Better name? Imagine 3 orthogonal vectors pointed in the x, y, and z
  directions.
  """
  return rec(
    [1, 0, 0,
     0, 1, 0,
     0, 0, 1],
     n = (3, 3)
    )

SUPPORTED_GEOMETRIES = OrderedDict([
  (3, [
    ("three_legs", _three_legs, [1, 1, 1], 15),
    ("trigonal_plane", _trigonal_plane, [3, 1, 1], 15),
    ]),
  (4, [
    ("tetrahedron", _tetrahedron, [1, 1, 1], 15),
    ("square_plane", _square_plane, [4, 1, 1], 15),
    ("trigonal_pyramid", _trigonal_pyramid, [3, 1, 1], 15),
    ("see_saw", _see_saw, [1, 1, 1], 15),
    ]),
  (5, [
    ("square_pyramid", _square_pyramid, [4, 1, 1], 15),
    ("trigonal_bipyramid", _trigonal_bipyramid, [3, 2, 2], 15),
    ("pentagon", _pentagon, [5, 1, 1], 15)
    ]),
  (6, [
    ("octahedron", _octahedron, [4, 4, 4], 20),
    ("trigonal_prism", _trigonal_prism, [3, 1, 1], 15),
    ("pentagonal_pyramid", _pentagonal_pyramid, [5, 1, 1], 15),
    ]),
  (7, [
    ("pentagonal_bipyramid", _pentagonal_bipyramid, [5, 2, 2], 15),
    ]),
  ])

SUPPORTED_GEOMETRY_NAMES = \
  [lst[0] for vals in SUPPORTED_GEOMETRIES.values() for lst in vals]

def _norm(vec):
  return sqrt(sum([i ** 2 for i in vec]))

def _rmsd_shape(a, b, power = 2):
  assert a.n == b.n
  c = list(b)
  distances = []
  for i in a:
    c.sort(key = lambda x: _norm(i - x))
    distances.append(_norm(i - c.pop(0)))
  distances = [i ** power for i in distances]
  return sqrt(sum(distances) / a.n[0])

def _angle(a, b, degree = False):
  from math import acos, pi
  angle = acos(a.dot(b) / (_norm(a) * _norm(b)))
  if degree: angle *= 180 / pi
  return angle

def _rmsa_shape(a, b, power = 2, degree = False):
  assert a.n == b.n
  c = b.as_list_of_lists()
  angles = []
  for i in a.as_list_of_lists():
    i = col(i)
    c.sort(key = lambda x: _angle(i, col(x), degree = degree))
    angles.append(_angle(i, col(c.pop(0)), degree = degree))

  angles = [i ** power for i in angles]
  rmsa = sqrt(sum(angles) / a.n[0])

  return rmsa

def _transform_shape(shape, scale, rotation):
  from math import sin, cos
  rot_x, rot_y, rot_z = rotation
  mat_x = rec([1, 0, 0,
               0, cos(rot_x), -sin(rot_x),
               0, sin(rot_x), cos(rot_x)],
              n = (3, 3))
  mat_y = rec([cos(rot_y), 0, sin(rot_y),
               0, 1, 0,
               -sin(rot_y), 0, cos(rot_y)],
              n = (3, 3))
  mat_z = rec([cos(rot_z), -sin(rot_z), 0,
               sin(rot_z), cos(rot_z), 0,
               0, 0, 1],
              n = (3, 3))
  mat_r = mat_x * mat_y * mat_z

  return rec([elem * scale for column in shape.as_list_of_lists()
              for elem in mat_r * col(column)], n = shape.n)

def sphere_minimizer(f, symmetry, step = 45):
  def sphere_iterate():
    for x in xrange(0, 360 // symmetry[0], step):
      for y in xrange(0, 360 // symmetry[1], step):
        for z in xrange(0, 360 // symmetry[2], step):
          yield col([x, y, z])
  return min(((x, f(x)) for x in sphere_iterate()), key = lambda x: x[-1])

def _minimize_points(shape, target, symmetry = None):
  def minimize_me(x):
    return _rmsa_shape(shape, _transform_shape(target, 1, x), degree = True)
  class minimize_me_class:
    def __init__(self): pass
    def target(self, vector): return minimize_me(vector)

  if symmetry is None:
    symmetry = [1, 1, 1]
  m = sphere_minimizer(minimize_me, symmetry)
  x = m[0]
  from scitbx.direct_search_simulated_annealing import dssa
  from scitbx.array_family import flex
  x = flex.double([i for i in x])
  starting_simplex = [
    x + flex.double((0, 0, 0)),
    x + flex.double((90, 0, 0)),
    x + flex.double((0, 90, 0)),
    x + flex.double((0, 0, 90)),
  ]
  dsa = dssa(
    dimension = 3,
    matrix = starting_simplex,
    evaluator = minimize_me_class(),
    further_opt = True,
    )
  x = dsa.get_solution()
  f_x = minimize_me(x)
  return f_x

def find_coordination_geometry(nearby_atoms, cutoff = 2.9):
  """
  Returns a list of tuples of potential geometries for the vectors given by
  nearby_atoms, along with the average deviation from those geometries.
  """

  # Filter out overlapping atoms, we just want an idea of the coordinating
  # geometry, even if it is two different atoms are occupying the same spot.
  nearby_atoms = [atom for index, atom in enumerate(nearby_atoms)
                  if not [other for other in nearby_atoms[index + 1:]
                          if atom.distance_from(other) < 0.5]]

  filtered = [contact for contact in nearby_atoms
              if contact.distance() < cutoff]
  vectors = rec(
    [i for contact in filtered for i in contact.vector],
    n = (len(filtered), 3))

  n_vectors = vectors.n[0]
  if n_vectors not in SUPPORTED_GEOMETRIES:
    return []

  rmsas = []
  for name, func, symmetry, rmsa_cutoff in SUPPORTED_GEOMETRIES[n_vectors]:
    rmsa = _minimize_points(vectors, func(), symmetry = symmetry)
    if rmsa < rmsa_cutoff:
      rmsas.append((name, rmsa))

  rmsas.sort(key = lambda x: x[-1])
  rmsas = [i for i in rmsas if i[-1] < rmsas[0][-1] * 1.5]
  return rmsas
