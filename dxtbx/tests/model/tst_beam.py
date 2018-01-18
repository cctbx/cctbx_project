from __future__ import absolute_import, division
from dxtbx.model import Beam

def tst_set_direction_wavelength():
  """Test setting direction and wavelength"""
  from scitbx import matrix
  direction = matrix.col((0.013142, 0.002200, 1.450476))
  unit_direction = direction.normalize()
  wavelength = 0.689400

  # Create the beam
  b = Beam(direction, wavelength)

  eps = 1e-7

  # Check direction is a unit vector
  assert(abs(matrix.col(b.get_direction()).length() - 1) <= eps)
  assert(abs(matrix.col(b.get_direction()) - unit_direction) <= eps)

  # Check wavelength is correct
  assert(abs(b.get_wavelength() - wavelength) <= eps)

  # Check s0 is in direction and has length 1/wavelength
  assert(abs(matrix.col(b.get_s0()).length() - 1.0 / wavelength) <= eps)
  assert(abs(-matrix.col(b.get_s0()).normalize() - unit_direction) <= eps)

  # Test passed
  print "OK"

def tst_set_s0():
  """Test setting s0"""
  from scitbx import matrix
  direction = matrix.col((0.013142, 0.002200, 1.450476))
  unit_direction = direction.normalize()
  wavelength = 0.689400
  s0 = -unit_direction * 1.0 / wavelength

  # Create the beam
  b = Beam(s0)

  eps = 1e-7

  # Check direction is a unit vector
  assert(abs(matrix.col(b.get_direction()).length() - 1) <= eps)
  assert(abs(matrix.col(b.get_direction()) - unit_direction) <= eps)

  # Check wavelength is correct
  assert(abs(b.get_wavelength() - wavelength) <= eps)

  # Check s0 is in direction and has length 1/wavelength
  assert(abs(matrix.col(b.get_s0()).length() - 1.0 / wavelength) <= eps)
  assert(abs(-matrix.col(b.get_s0()).normalize() - unit_direction) <= eps)
  assert(abs(matrix.col(b.get_s0()) - s0) <= eps)

  # Test passed
  print "OK"

def tst_from_phil():

  from libtbx.phil import parse
  from dxtbx.model.beam import beam_phil_scope, BeamFactory
  from scitbx import matrix
  direction = matrix.col((0.013142, 0.002200, 1.450476))
  unit_direction = direction.normalize()
  wavelength = 0.689400

  reference = Beam(direction, wavelength)

  params1 = beam_phil_scope.fetch(parse("""
    beam {
      wavelength = 1.0
      direction = (0, 0, 1)
    }
  """)).extract()

  params2 = beam_phil_scope.fetch(parse("""
    beam {
      wavelength = 1.0
    }
  """)).extract()

  # Create the beam
  b1 = BeamFactory.from_phil(params1)
  b2 = BeamFactory.from_phil(params2, reference)
  try:
    b3 = BeamFactory.from_phil(params2)
    passed = True
  except Exception:
    passed = False
  assert passed == False

  print 'OK'

def tst_scan_varying():

  from scitbx import matrix
  direction = matrix.col((0.013142, 0.002200, 1.450476))
  unit_direction = direction.normalize()
  wavelength = 0.689400
  s0 = -unit_direction * 1.0 / wavelength

  # Create the beam
  b = Beam(s0)

  assert b.get_num_scan_points() == 0
  assert b.get_s0_at_scan_points().size() == 0
  try:
    b.get_s0_at_scan_point(0) # should raise RuntimeError
  except RuntimeError:
    pass

  # set varying beam
  num_scan_points = 11
  s0_static = matrix.col(b.get_s0())
  s0_as_scan_points = [s0_static]
  axis = matrix.col.random(3, -1., 1.).normalize()
  for i in range(num_scan_points-1):
    s0_as_scan_points.append(
      s0_as_scan_points[-1].rotate_around_origin(axis, angle=0.01, deg=True))
  b.set_s0_at_scan_points(s0_as_scan_points)
  assert b.get_num_scan_points() == 11
  assert b.get_s0_at_scan_points().size() == 11

  for t in range(num_scan_points):
    s0_t = matrix.col(b.get_s0_at_scan_point(t))
    assert s0_t == s0_as_scan_points[t]

  # also test setting as tuple
  b.set_s0_at_scan_points(tuple(s0_as_scan_points))
  assert b.get_num_scan_points() == 11
  assert b.get_s0_at_scan_points().size() == 11

  # test resetting
  b.reset_scan_points()
  assert b.get_num_scan_points() == 0
  assert b.get_s0_at_scan_points().size() == 0

  print 'OK'


def run():
  """Test the beam object"""
  tst_set_direction_wavelength()
  tst_set_s0()
  tst_from_phil()
  tst_scan_varying()

if __name__ == '__main__':
  run()
