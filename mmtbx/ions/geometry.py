# -*- coding: utf-8; py-indent-offset: 2 -*-
"""
Methods for examining bond angles around an ion and categorizing
the coordination geometry.
"""
from __future__ import division
from collections import OrderedDict
from math import sqrt
from scitbx.matrix import col, rec

def _bond_angles(vectors):
  """
  Returns a list of angles (In degrees) between all two-element combinations
  in vectors.
  """

  return [(v1, v2, v1.angle(v2, deg = True))
          for index, v1 in enumerate(vectors)
          for v2 in vectors[index + 1:]]

def _is_tetrahedron(vectors, dev_cutoff = 20):
  """
  Tetrahedrons have four vertices, with angles between all pairs of vertices
  uniformly about 104.5 degrees.
  """

  if len(vectors) > 4 or len(vectors) < 3:
    return

  angles = _bond_angles(vectors)
  deviation = sqrt(sum(abs(i[2] - 104.5) ** 2 for i in angles) / len(vectors))

  if deviation <= dev_cutoff:
    return deviation, 4 - len(vectors)

def _is_trigonal_plane(vectors, dev_cutoff = 20):
  """
  Triangular planar geometry has three vertices (By definition all on the same
  equatorial plane). The expected angles are 120 degrees between neighboring
  vertices.
  """
  if len(vectors) != 3:
    return

  angles = _bond_angles(vectors)
  a_120s = []

  for angle in angles:
    a_120s.append(angle[2] - 120)

  deviation = sqrt(sum(i ** 2 for i in a_120s) / len(angles))

  if deviation <= dev_cutoff:
    return deviation, 3 - len(vectors)

def _is_square_plane(vectors, dev_cutoff = 20):
  """
  Square planar geometry has four vertices, all on the same equatorial plane.
  The expected angles are 90 degrees between neighboring vertices and 180
  degrees between vertices across from one another.
  """

  if len(vectors) != 4:
    return

  angles = _bond_angles(vectors)

  # Expect 2x 180 degrees and 4x 90 degrees
  a_90s = []
  a_180s = []

  for angle in angles:
    if abs(angle[2] - 90) < abs(angle[2] - 180):
      a_90s.append(angle[2] - 90)
    else:
      a_180s.append(angle[2] - 180)

  # With up to one atom missing, we must have 2 to 4 90 degree angles and 1 to 2
  # 180 degree angles
  if len(a_90s) < 2 or len(a_90s) > 4 or len(a_180s) < 1 or len(a_180s) > 2:
    return

  deviation = sqrt(sum(i ** 2 for i in a_90s + a_180s) / len(angles))

  if deviation <= dev_cutoff:
    return deviation, 4 - len(vectors)

def _is_square_pyramid(vectors, dev_cutoff = 20):
  """
  Square bipyramids have five vertices, four on the same equatorial plane with
  one above. The expected angles are all either 90 degrees or 180 degrees.
  """
  if len(vectors) != 5:
    return

  angles = _bond_angles(vectors)
  a_90s, a_180s = [], []

  for angle in angles:
    if abs(angle[2] - 90) < abs(angle[2] - 180):
      a_90s.append(angle[2] - 90)
    else:
      a_180s.append(angle[2] - 180)

  if len(a_90s) != 8 or len(a_180s) != 2:
    return

  deviation = sqrt(sum(i ** 2 for i in a_90s + a_180s) / len(angles))

  if deviation <= dev_cutoff:
    return deviation, 5 - len(vectors)

def _is_octahedron(vectors, dev_cutoff = 20):
  """
  Octahedrons have six vertices (Their name comes from their eight faces).
  The expected angles are all either 90 degrees (Next to each other),
  or 180 degrees (Across from each other).

  Another name for this shape is square bipyramidal.
  """
  if len(vectors) != 6:
    return

  angles = _bond_angles(vectors)
  a_90s, a_180s = [], []

  for angle in angles:
    if abs(angle[-1] - 90) < abs(angle[-1] - 180):
      a_90s.append(angle[-1] - 90)
    else:
      a_180s.append(angle[-1] - 180)

  if len(a_180s) > 3 or len(a_180s) < 2 or len(a_90s) < 8 or len(a_90s) > 12:
    return

  deviation = sqrt(sum(i ** 2 for i in a_90s + a_180s) / len(angles))

  if deviation <= dev_cutoff:
    return deviation, 6 - len(vectors)

def _is_trigonal_pyramid(vectors, dev_cutoff = 15):
  """
  Trigional pyramids have four vertices. Three vertices form a plane with
  angles of 120 degrees between each pair. The last vertex resides axial
  to the plane, at 90 degrees from all of the equatorial vertices.
  """
  if len(vectors) != 4:
    return

  angles = _bond_angles(vectors)
  a_90s, a_120s = [], []

  for angle in angles:
    if abs(angle[2] - 90) < abs(angle[2] - 120):
      a_90s.append(angle[2] - 90)
    else:
      a_120s.append(angle[2] - 120)

  if len(a_90s) < 2 or len(a_90s) > 4 or len(a_120s) < 2 or len(a_120s) > 4:
    return

  deviation = sqrt(sum(i ** 2 for i in a_90s + a_120s) / len(angles))

  if deviation <= dev_cutoff:
    return deviation, 4 - len(vectors)

def _is_trigonal_bipyramid(vectors, dev_cutoff = 15):
  """
  Trigonal bipyramids have five vertices. Three vertices form a plane in the
  middle and the angles between all three are 120 degrees. The two other
  vertices reside axial to the plane, at 90 degrees from all the equatorial
  vertices.
  """
  if len(vectors) > 5 or len(vectors) < 4:
    return

  angles = _bond_angles(vectors)

  # Grab the two axial vectors
  ax1, ax2, axial_angle = max(angles, key = lambda x: abs(x[-1]))

  if axial_angle < 150:
    # Missing one of the two axial vectors, just quit
    return

  base_to_axials = []
  equatorial_angles = []

  for v1, v2, angle in angles:
    # Python has no boolean xor!
    # Grab the angles between the two endpoints of the bipyramid and the base
    if (v1 in [ax1, ax2]) != (v2 in [ax1, ax2]):
      base_to_axials += angle,
    elif (v1 not in [ax1, ax2]) and (v2 not in [ax1, ax2]):
      equatorial_angles += angle,

  deviants =  [axial_angle - 180]
  deviants += [i - 90 for i in base_to_axials]
  deviants += [i - 120 for i in equatorial_angles]
  deviation = sqrt(sum(i ** 2 for i in deviants) / len(deviants))

  if deviation <= dev_cutoff:
    return deviation, 5 - len(vectors)

def _is_pentagonal_bipyramid(vectors, dev_cutoff = 15):
  """
  Pentagonal bipyramids have seven vertices. Five vertices form a plane in the
  middle and the angles between all five are 72 degrees. The two other vertices
  reside axial to the plane, at 90 degrees from all the equatorial vertices.
  """
  if len(vectors) > 7 or len(vectors) < 6:
    return

  angles = _bond_angles(vectors)

  # Determine which two vectors define the axial angles
  axials = []
  for v1 in vectors:
    v_angles = [v1.angle(v2, deg = True) for v2 in vectors
                if v2 != v1]
    a_180s = len([i for i in v_angles if abs(i - 180) < 20])
    a_90s = len([i for i in v_angles if abs(i - 90) < 20])

    if a_180s > 0  and a_90s > 4:
      axials.append(v1)

  if len(axials) != 2:
    # Couldn't determine axial angles
    return

  ax1, ax2 = axials
  axial_angle = ax1.angle(ax2, deg = True)

  base_to_axials = []
  equatorial_angles = []

  for v1, v2, angle in angles:
    # Python has no boolean xor!
    # Grab the angles between the two endpoints of the bipyramid and the base
    if (v1 in [ax1, ax2]) != (v2 in [ax1, ax2]):
      base_to_axials += angle,
    elif (v1 not in [ax1, ax2]) and (v2 not in [ax1, ax2]):
      equatorial_angles += angle,

  deviants =  [axial_angle - 180]
  deviants += [i - 90 for i in base_to_axials]
  deviants += [min(abs(i - 72), abs(i - 144)) for i in equatorial_angles]
  deviation = sqrt(sum(i ** 2 for i in deviants) / len(deviants))

  if deviation <= dev_cutoff:
    return deviation, 7 - len(vectors)

def _is_trigonal_prism(vectors, dev_cutoff = 15):
  """
  Triangular prisms are defined by 3 vertices in a triangular pattern on two
  aligned planes. Unfortunately, the angles are dependent on the length and
  width of the prism. Need more examples to come up with a better way of
  detecting this shape.

  For now, this code is experimental.
  """
  if len(vectors) != 6:
    return

  angles = _bond_angles(vectors)
  a_85s, a_135s = [], []

  for angle in angles:
    if abs(angle[-1] - 85) < abs(angle[-1] - 135):
      a_85s.append(angle[-1] - 85)
    else:
      a_135s.append(angle[-1] - 135)

  if len(a_85s) != 9 and len(a_135s) != 6:
    return

  deviation = sqrt(sum(i ** 2 for i in a_85s + a_135s) / len(angles))

  if deviation < dev_cutoff:
    return deviation, 6 - len(vectors)

SUPPORTED_GEOMETRIES_OLD = OrderedDict([
  ("tetrahedron", _is_tetrahedron),
  ("trigonal_plane", _is_trigonal_plane),
  ("square_plane", _is_square_plane),
  ("square_pyramid", _is_square_pyramid),
  ("octahedron", _is_octahedron),
  ("trigonal_pyramid", _is_trigonal_pyramid),
  ("trigonal_bipyramid", _is_trigonal_bipyramid),
  ("pentagonal_bipyramid", _is_pentagonal_bipyramid),
  ("trigonal_prism", _is_trigonal_prism),
  ])

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

def find_coordination_geometry(nearby_atoms, minimizer_method = False,
                               cutoff = 2.9):
  """
  Returns a list of tuples of potential geometries for the vectors given by
  nearby_atoms, along with the average deviation from those geometries.
  """

  # Filter out overlapping atoms, we just want an idea of the coordinating
  # geometry, even if it is two different atoms are occupying the same spot.
  nearby_atoms = [atom for index, atom in enumerate(nearby_atoms)
                  if not [other for other in nearby_atoms[index + 1:]
                          if atom.distance_from(other) < 0.5]]

  filtered = [contact.vector for contact in nearby_atoms
              if contact.distance() < cutoff]
  geometries = []
  if minimizer_method:
    vectors = rec(
      [i for vector in filtered for i in vector],
      n = (len(filtered), 3))

    n_vectors = vectors.n[0]
    if n_vectors not in SUPPORTED_GEOMETRIES:
      return geometries

    for name, func, symmetry, rmsa_cutoff in SUPPORTED_GEOMETRIES[n_vectors]:
      rmsa = _minimize_points(vectors, func(), symmetry = symmetry)
      if rmsa < rmsa_cutoff:
        geometries.append((name, rmsa))

    geometries.sort(key = lambda x: x[-1])
    geometries = [i for i in geometries if i[-1] < geometries[0][-1] * 1.5]
  else:
    for name, func in SUPPORTED_GEOMETRIES_OLD.items():
      val = func(filtered)
      if val:
        deviation, missing = val
        geometries.append((name, deviation))

  return geometries
