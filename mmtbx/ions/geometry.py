# -*- coding: utf-8; py-indent-offset: 2 -*-
"""
This module provides tools for examining a set of vectors and find the geometry
that best fits from a set of built in shapes.
"""

from __future__ import absolute_import, division, print_function
from scitbx.matrix import col
from collections import OrderedDict
try:
  from collections.abc import Iterable
except ImportError:
  from collections import Iterable
from math import sqrt
from six.moves import zip


def _bond_angles(vectors):
  """
  Creates a list of angles (In degrees) between all two-element combinations
  in vectors.

  Parameters
  ----------
  vectors : scitbx.matrix.col

  Returns
  -------
  list of float
  """

  return [(v1, v2, v1.angle(v2, deg=True))
          for index, v1 in enumerate(vectors)
          for v2 in vectors[index + 1:]]

def _is_tetrahedron(vectors, dev_cutoff=20):
  """
  Tetrahedrons have four vertices, with angles between all pairs of vertices
  uniformly about 104.5 degrees.

  Parameters
  ----------
  vectors : list scitbx.matrix.col
  dev_cutoff : float, optional

  Returns
  -------
  bool
  """

  if len(vectors) > 4 or len(vectors) < 3:
    return

  angles = _bond_angles(vectors)
  deviation = sqrt(sum(abs(i[2] - 104.5) ** 2 for i in angles) / len(vectors))

  if deviation <= dev_cutoff:
    return deviation, 4 - len(vectors)

def _is_trigonal_plane(vectors, dev_cutoff=20):
  """
  Triangular planar geometry has three vertices (By definition all on the same
  equatorial plane). The expected angles are 120 degrees between neighboring
  vertices.

  Parameters
  ----------
  vectors : list scitbx.matrix.col
  dev_cutoff : float, optional

  Returns
  -------
  bool
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

def _is_square_plane(vectors, dev_cutoff=20):
  """
  Square planar geometry has four vertices, all on the same equatorial plane.
  The expected angles are 90 degrees between neighboring vertices and 180
  degrees between vertices across from one another.

  Parameters
  ----------
  vectors : list scitbx.matrix.col
  dev_cutoff : float, optional

  Returns
  -------
  bool
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

def _is_square_pyramid(vectors, dev_cutoff=20):
  """
  Square bipyramids have five vertices, four on the same equatorial plane with
  one above. The expected angles are all either 90 degrees or 180 degrees.

  Parameters
  ----------
  vectors : list scitbx.matrix.col
  dev_cutoff : float, optional

  Returns
  -------
  bool
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

def _is_octahedron(vectors, dev_cutoff=20):
  """
  Octahedrons have six vertices (Their name comes from their eight faces).
  The expected angles are all either 90 degrees (Next to each other),
  or 180 degrees (Across from each other).

  Another name for this shape is square bipyramidal.

  Parameters
  ----------
  vectors : list scitbx.matrix.col
  dev_cutoff : float, optional

  Returns
  -------
  bool
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

def _is_trigonal_pyramid(vectors, dev_cutoff=15):
  """
  Trigional pyramids have four vertices. Three vertices form a plane with
  angles of 120 degrees between each pair. The last vertex resides axial
  to the plane, at 90 degrees from all of the equatorial vertices.

  Parameters
  ----------
  vectors : list scitbx.matrix.col
  dev_cutoff : float, optional

  Returns
  -------
  bool
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

def _is_trigonal_bipyramid(vectors, dev_cutoff=15):
  """
  Trigonal bipyramids have five vertices. Three vertices form a plane in the
  middle and the angles between all three are 120 degrees. The two other
  vertices reside axial to the plane, at 90 degrees from all the equatorial
  vertices.

  Parameters
  ----------
  vectors : list scitbx.matrix.col
  dev_cutoff : float, optional

  Returns
  -------
  bool
  """
  if len(vectors) > 5 or len(vectors) < 4:
    return

  angles = _bond_angles(vectors)

  # Grab the two axial vectors
  ax1, ax2, axial_angle = max(angles, key=lambda x: abs(x[-1]))

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

  deviants = [axial_angle - 180]
  deviants += [i - 90 for i in base_to_axials]
  deviants += [i - 120 for i in equatorial_angles]
  deviation = sqrt(sum(i ** 2 for i in deviants) / len(deviants))

  if deviation <= dev_cutoff:
    return deviation, 5 - len(vectors)

def _is_pentagonal_bipyramid(vectors, dev_cutoff=15):
  """
  Pentagonal bipyramids have seven vertices. Five vertices form a plane in the
  middle and the angles between all five are 72 degrees. The two other vertices
  reside axial to the plane, at 90 degrees from all the equatorial vertices.

  Parameters
  ----------
  vectors : list scitbx.matrix.col
  dev_cutoff : float, optional

  Returns
  -------
  bool
  """
  if len(vectors) > 7 or len(vectors) < 6:
    return

  angles = _bond_angles(vectors)

  # Determine which two vectors define the axial angles
  axials = []
  for v1 in vectors:
    v_angles = []
    for v2 in vectors:
      if v2 != v1:
        v_angles.append(v1.angle(v2, deg=True))

    a_180s = len([i for i in v_angles if abs(i - 180) < 20])
    a_90s = len([i for i in v_angles if abs(i - 90) < 20])

    if a_180s > 0  and a_90s > 4:
      axials.append(v1)

  if len(axials) != 2:
    # Couldn't determine axial angles
    return

  ax1, ax2 = axials
  axial_angle = ax1.angle(ax2, deg=True)

  base_to_axials = []
  equatorial_angles = []

  for v1, v2, angle in angles:
    # Python has no boolean xor!
    # Grab the angles between the two endpoints of the bipyramid and the base
    if (v1 in [ax1, ax2]) != (v2 in [ax1, ax2]):
      base_to_axials += angle,
    elif (v1 not in [ax1, ax2]) and (v2 not in [ax1, ax2]):
      equatorial_angles += angle,

  deviants = [axial_angle - 180]
  deviants += [i - 90 for i in base_to_axials]
  deviants += [min(abs(i - 72), abs(i - 144)) for i in equatorial_angles]
  deviation = sqrt(sum(i ** 2 for i in deviants) / len(deviants))

  if deviation <= dev_cutoff:
    return deviation, 7 - len(vectors)

def _is_trigonal_prism(vectors, dev_cutoff=15):
  """
  Triangular prisms are defined by 3 vertices in a triangular pattern on two
  aligned planes. Unfortunately, the angles are dependent on the length and
  width of the prism. Need more examples to come up with a better way of
  detecting this shape.

  For now, this code is experimental.

  Parameters
  ----------
  vectors : list scitbx.matrix.col
  dev_cutoff : float, optional

  Returns
  -------
  bool
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
  ("tetrahedral", _is_tetrahedron),
  ("trigonal_planar", _is_trigonal_plane),
  ("square_planar", _is_square_plane),
  ("square_pyramidal", _is_square_pyramid),
  ("octahedral", _is_octahedron),
  ("trigonal_pyramidal", _is_trigonal_pyramid),
  ("trigonal_bipyramidal", _is_trigonal_bipyramid),
  ("pentagonal_bipyramidal", _is_pentagonal_bipyramid),
  ("trigonal_prism", _is_trigonal_prism),
  ])

def _concatenate(*args):
  """
  Reduces a list of a mixture of elements and lists down to a single list.

  Parameters
  ----------
  args : tuple of (object or list of object)

  Returns
  -------
  list
  """
  lst = []
  for arg in args:
    if isinstance(arg, list):
      for elem in arg:
        lst.append(elem)
    else:
      lst.append(arg)
  return lst

def _tetrahedron():
  """
  Returns
  -------
  list of scitbx.matrix.col
  """
  return [
    col([1, 1, 1]),
    col([-1, -1, 1]),
    col([-1, 1, -1]),
    col([1, -1, -1]),
    ]

def _octahedron():
  """
  Returns
  -------
  list of scitbx.matrix.col
  """
  return _bipyramid(_square_plane())

def _trigonal_plane():
  """
  Returns
  -------
  list of scitbx.matrix.col
  """
  return [
    col([0, 1, 0]),
    col([-sqrt(3) / 2, -1 / 2, 0]),
    col([sqrt(3) / 2, -1 / 2, 0]),
    ]

def _square_plane():
  """
  Returns
  -------
  list of scitbx.matrix.col
  """
  return [
    col([0, 1, 0]),
    col([0, -1, 0]),
    col([1, 0, 0]),
    col([-1, 0, 0]),
    ]

def _pyramid(base):
  """
  Returns
  -------
  list of scitbx.matrix.col
  """
  return _concatenate(base, col([0, 0, 1]))

def _bipyramid(base):
  """
  Returns
  -------
  list of scitbx.matrix.col
  """
  return _concatenate(base, col([0, 0, 1]), col([0, 0, -1]))

def _square_pyramid():
  """
  Returns
  -------
  list of scitbx.matrix.col
  """
  return _pyramid(_square_plane())

def _square_bipyramid():
  """
  Returns
  -------
  list of scitbx.matrix.col
  """
  return _bipyramid(_square_plane())

def _trigonal_pyramid():
  """
  Returns
  -------
  list of scitbx.matrix.col
  """
  return _pyramid(_trigonal_plane())

def _trigonal_bipyramid():
  """
  Returns
  -------
  list of scitbx.matrix.col
  """
  return _bipyramid(_trigonal_plane())

def _trigonal_prism():
  """
  Returns
  -------
  list of scitbx.matrix.col
  """
  return _concatenate(
    [i + col([0, 0, 1]) for i in _trigonal_plane()],
    [i + col([0, 0, -1]) for i in _trigonal_plane()],
    )

def _pentagon():
  """
  Create a list of vectors in the shape of a planar pentagon.

  Returns
  -------
  list of scitbx.matrix.col

  See Also
  --------
  http://mathworld.wolfram.com/Pentagon.html
  """
  c_1 = (sqrt(5) - 1) / 4
  c_2 = (sqrt(5) + 1) / 4
  s_1 = sqrt(10 + 2 * sqrt(5)) / 4
  s_2 = sqrt(10 - 2 * sqrt(5)) / 4

  return [
    col([1, 0, 0]),
    col([c_1, s_1, 0]),
    col([-c_2, s_2, 0]),
    col([-c_2, -s_2, 0]),
    col([c_1, -s_1, 0]),
    ]

def _pentagonal_pyramid():
  """
  Returns
  -------
  list of scitbx.matrix.col
  """
  return _pyramid(_pentagon())

def _pentagonal_bipyramid():
  """
  Creates a square bipyramid shape.

  Returns
  -------
  list of scitbx.matrix.col
  """
  return _bipyramid(_pentagon())

def _square_pyramid_bidentate_miss_1():
  """
  Creates a square pyramid shape with one vertex replaced with a bidentate
  coordination group. One vertex is missing in this shape.

  Returns
  -------
  list of scitbx.matrix.col
  """
  return _concatenate(
    _square_plane(),
    col([sqrt(2) / 2, sqrt(2) / 2, -1]),
    col([-sqrt(2) / 2, -sqrt(2) / 2, -1]),
    )

def _square_pyramid_bidentate_miss_2():
  """
  Creates a square pyramid shape with one vertex replaced with a bidentate
  coordination group. One vertex is missing in this shape.

  Returns
  -------
  list of scitbx.matrix.col
  """
  return [
    col([0, 1, 0]),
    col([0, -1, 0]),
    col([-1, 0, 0]),
    col([0, 0, 1]),
    col([sqrt(2) / 2, sqrt(2) / 2, -1]),
    col([-sqrt(2) / 2, -sqrt(2) / 2, -1]),
    ]

def _square_pyramid_bidentate_miss_3():
  """
  Creates a square pyramid shape with one vertex replaced with a bidentate
  coordination group. One vertex is missing in this shape.

  Returns
  -------
  list of scitbx.matrix.col
  """
  return [
    col([0, 1, 0]),
    col([0, -1, 0]),
    col([-1, 0, 0]),
    col([1, 0, 0]),
    col([sqrt(2) / 2, sqrt(2) / 2, -1]),
    col([-sqrt(2) / 2, -sqrt(2) / 2, -1]),
    ]

def _square_pyramid_bidentate():
  """
  Creates a square pyramid shape with one vertex replaced with a bidentate
  coordination group.

  Returns
  -------
  list of scitbx.matrix.col
  """
  return _concatenate(
    _square_pyramid(),
    col([sqrt(2) / 2, sqrt(2) / 2, -1]),
    col([-sqrt(2) / 2, -sqrt(2) / 2, -1])
    )

def _pentagonal_pyramid_bidentate():
  """
  Creates a pentagonal pyramid shape with one vertex replaced with a bidentate
  coordination group.

  Returns
  -------
  list of scitbx.matrix.col
  """
  return _concatenate(
    _pentagonal_pyramid(),
    col([sqrt(2) / 2, sqrt(2) / 2, -1]),
    col([-sqrt(2) / 2, -sqrt(2) / 2, -1]),
    )

def _pentagonal_bibidentate_miss_1():
  """
  A planar pentagon with bidentate atoms coordinating directly above and
  below. One vertex from the plane is missing in this shape.

  Returns
  -------
  list of scitbx.matrix.col
  """
  return _concatenate(
    col([sqrt(2) / 2, sqrt(2) / 2, 1]),
    col([-sqrt(2) / 2, -sqrt(2) / 2, 1]),
    col([sqrt(2) / 2, sqrt(2) / 2, -1]),
    col([-sqrt(2) / 2, -sqrt(2) / 2, -1]),
    *_pentagon()[:-1]
    )

def _pentagonal_bibidentate_miss_2():
  """
  A planar pentagon with bidentate atoms coordinating directly above and
  below. One vertex from the plane is missing in this shape.

  Returns
  -------
  list of scitbx.matrix.col
  """
  return _concatenate(
    col([sqrt(2) / 2, sqrt(2) / 2, 1]),
    col([-sqrt(2) / 2, -sqrt(2) / 2, 1]),
    col([sqrt(2) / 2, sqrt(2) / 2, -1]),
    col([-sqrt(2) / 2, -sqrt(2) / 2, -1]),
    *_pentagon()[1:]
    )

def _pentagonal_bibidentate_miss_3():
  """
  A planar pentagon with bidentate atoms coordinating directly above and
  below. One vertex from the plane is missing in this shape.

  Returns
  -------
  list of scitbx.matrix.col
  """
  pentagon = _pentagon()
  return _concatenate(
    col([sqrt(2) / 2, sqrt(2) / 2, 1]),
    col([-sqrt(2) / 2, -sqrt(2) / 2, 1]),
    col([sqrt(2) / 2, sqrt(2) / 2, -1]),
    col([-sqrt(2) / 2, -sqrt(2) / 2, -1]),
    pentagon[0],
    pentagon[1],
    pentagon[3],
    pentagon[4],
    )

def _pentagonal_bibidentate_miss_4():
  """
  A planar pentagon with bidentate atoms coordinating directly above and
  below. One vertex from a bidentate coordinator is missing in this shape.

  Returns
  -------
  list of scitbx.matrix.col
  """
  return _concatenate(
    _pentagon(),
    col([sqrt(2) / 2, sqrt(2) / 2, 1]),
    col([-sqrt(2) / 2, -sqrt(2) / 2, 1]),
    col([-sqrt(2) / 2, -sqrt(2) / 2, -1]),
    )

def _pentagonal_bibidentate_miss_5():
  """
  A planar pentagon with bidentate atoms coordinating directly above and
  below. One vertex from a bidentate coordinator is missing in this shape.

  Returns
  -------
  list of scitbx.matrix.col
  """
  return _concatenate(
    _pentagon(),
    col([sqrt(2) / 2, sqrt(2) / 2, 1]),
    col([-sqrt(2) / 2, -sqrt(2) / 2, 1]),
    col([sqrt(2) / 2, sqrt(2) / 2, -1]),
    )

def _pentagonal_bibidentate():
  """
  A planar pentagon with bidentate atoms coordinating directly above and below.

  Returns
  -------
  list of scitbx.matrix.col
  """
  return _concatenate(
    _pentagon(),
    col([sqrt(2) / 2, sqrt(2) / 2, 1]),
    col([-sqrt(2) / 2, -sqrt(2) / 2, 1]),
    col([sqrt(2) / 2, sqrt(2) / 2, -1]),
    col([-sqrt(2) / 2, -sqrt(2) / 2, -1]),
    )

def _see_saw():
  """
  An octahedron missing two adjacent points.

  Returns
  -------
  list of scitbx.matrix.col
  """
  return [
    col([1, 0, 0]),
    col([-1, 0, 0]),
    col([0, 1, 0]),
    col([0, 0, 1]),
    ]

def _three_legs():
  """
  Better name? Imagine 3 orthogonal vectors pointed in the x, y, and z
  directions.

  Returns
  -------
  list of scitbx.matrix.col
  """
  return [
    col([1, 0, 0]),
    col([0, 1, 0]),
    col([0, 0, 1]),
    ]

SUPPORTED_GEOMETRIES = OrderedDict([
  (3, [
    ("three_legs", _three_legs, 15),
    ("trigonal_plane", _trigonal_plane, 15),
    ]),
  (4, [
    ("tetrahedron", _tetrahedron, 15),
    ("square_plane", _square_plane, 20),
    ("trigonal_pyramid", _trigonal_pyramid, 15),
    ("see_saw", _see_saw, 15),
    ]),
  (5, [
    ("square_pyramid", _square_pyramid, 15),
    ("trigonal_bipyramid", _trigonal_bipyramid, 15),
    ("pentagon", _pentagon, 15)
    ]),
  (6, [
    ("octahedron", _octahedron, 15),
    ("trigonal_prism", _trigonal_prism, 15),
    ("pentagonal_pyramid", _pentagonal_pyramid, 15),
    ("square_pyramid_bidentate_miss", [_square_pyramid_bidentate_miss_1,
                                       _square_pyramid_bidentate_miss_2,
                                       _square_pyramid_bidentate_miss_3], 15),
    ]),
  (7, [
    ("pentagonal_bipyramid", _pentagonal_bipyramid, 15),
    ("square_pyramid_bidentate", _square_pyramid_bidentate, 15),
    ]),
  (8, [
    ("pentagonal_pyramid_bidentate", _pentagonal_pyramid_bidentate, 15),
    ("pentagonal_bibidentate_miss", [_pentagonal_bibidentate_miss_1,
                                     _pentagonal_bibidentate_miss_2,
                                     _pentagonal_bibidentate_miss_3], 15),
  ]),
  (9, [
    ("pentagonal_bibidentate", _pentagonal_bibidentate, 15),
  ]),
  ])

SUPPORTED_GEOMETRY_NAMES = \
  [lst_i[0] for vals in SUPPORTED_GEOMETRIES.values() for lst_i in vals]

def _angles_deviation(vectors_a, vectors_b):
  """
  Calculates the root mean square of the angle deviation (in degrees) between
  two lists of vectors.

  Parameters
  ----------
  vectors_a : list of scitbx.matrix.col
  vectors_b : list of scitbx.matrix.col

  Returns
  -------
  float
  """
  assert len(vectors_a) == len(vectors_b)

  angles_a = [vec.angle(vec_o, deg=True)
              for index, vec in enumerate(vectors_a)
              for vec_o in vectors_a[index + 1:]]

  angles_b = [vec.angle(vec_o, deg=True)
              for index, vec in enumerate(vectors_b)
              for vec_o in vectors_b[index + 1:]]

  angles_a.sort()
  angles_b.sort()

  angle_deviation = sqrt(sum((i - j) ** 2 for i, j in zip(angles_a, angles_b))
                   / len(angles_a))

  return angle_deviation

def find_coordination_geometry(nearby_atoms, minimizer_method=False,
                               cutoff=2.9):
  """
  Searches through a list of geometries to find those that fit nearby_atom.

  Geometries are recognized by generating a list of all combinations of angles
  between the vectors and comparing them against the angles among the vectors
  of the ideal geometry.

  Parameters
  ----------
  nearby_atoms: list of mmtbx.ions.environment.atom_contact
      A list of atom contacts, indicating the vertices of the shape to be
      recognized.
  minimizer_method: bool, optional
      Optional parameter to use the new, more efficient version of geometry
      recognition. The old method will be depreciated in later versions of
      cctbx.
  cutoff: float, optional
      A cutoff distance, past which vectors are not included in geometry
      calculations.

  Returns
  -------
  list of tuples of str, float
      A list of found geometries. Each tuple contains the name of the geometry
      in string form followed by the deviation from ideal angles.

  See Also
  --------
  mmtbx.ions.geometry.SUPPORTED_GEOMETRY_NAMES,
  mmtbx.ions.geometry.SUPPORTED_GEOMETRIES_OLD
  """

  # Filter out overlapping atoms, we just want an idea of the coordinating
  # geometry, even if it is two different atoms are occupying the same spot.
  non_overlapping = []
  for index, contact in enumerate(nearby_atoms):
    if all(contact.distance_from(other) > 0.5
           for other in nearby_atoms[index + 1:]):
      non_overlapping.append(contact)

  # Filter out contacts > cutoff away
  filtered = []
  for contact in non_overlapping:
    if contact.distance() < cutoff:
      filtered.append(contact.vector)

  geometries = []
  if minimizer_method:
    n_vectors = len(filtered)

    if n_vectors not in SUPPORTED_GEOMETRIES:
      return geometries

    for name, func, rmsa_cutoff in SUPPORTED_GEOMETRIES[n_vectors]:
      if isinstance(func, Iterable):
        rmsa = min(_angles_deviation(filtered, i()) for i in func)
      else:
        rmsa = _angles_deviation(filtered, func())

      if rmsa < rmsa_cutoff:
        geometries.append((name, rmsa))

    if geometries:
      geometries.sort(key=lambda x: x[-1])
      geometries = [geometries[0]]
  else:
    for name, func in SUPPORTED_GEOMETRIES_OLD.items():
      val = func(filtered)
      if val:
        deviation = val[0]
        geometries.append((name, deviation))

  return geometries
