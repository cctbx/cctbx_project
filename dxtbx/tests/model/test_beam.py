from __future__ import absolute_import, division, print_function

from dxtbx.model import Beam
import pytest

def test_setting_direction_and_wavelength():
  from scitbx import matrix
  direction = matrix.col((0.013142, 0.002200, 1.450476))
  unit_direction = direction.normalize()
  wavelength = 0.689400

  # Create the beam
  b = Beam(direction, wavelength)

  eps = 1e-7

  # Check direction is a unit vector
  assert matrix.col(b.get_direction()).length() == pytest.approx(1)
  assert abs(matrix.col(b.get_direction()) - unit_direction) <= eps

  # Check wavelength is correct
  assert b.get_wavelength() == pytest.approx(wavelength)

  # Check s0 is in direction and has length 1/wavelength
  assert matrix.col(b.get_s0()).length() == pytest.approx(1.0 / wavelength)
  assert abs(-matrix.col(b.get_s0()).normalize() - unit_direction) <= eps

def test_setting_s0():
  from scitbx import matrix
  direction = matrix.col((0.013142, 0.002200, 1.450476))
  unit_direction = direction.normalize()
  wavelength = 0.689400
  s0 = -unit_direction * 1.0 / wavelength

  # Create the beam
  b = Beam(s0)

  eps = 1e-7

  # Check direction is a unit vector
  assert matrix.col(b.get_direction()).length() == pytest.approx(1)
  assert abs(matrix.col(b.get_direction()) - unit_direction) <= eps

  # Check wavelength is correct
  assert b.get_wavelength() == pytest.approx(wavelength)

  # Check s0 is in direction and has length 1/wavelength
  assert matrix.col(b.get_s0()).length() == pytest.approx(1.0 / wavelength)
  assert abs(-matrix.col(b.get_s0()).normalize() - unit_direction) <= eps
  assert abs(matrix.col(b.get_s0()) - s0) <= eps

def test_from_phil():
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
  with pytest.raises(RuntimeError):
    b3 = BeamFactory.from_phil(params2)

  params3 = beam_phil_scope.fetch(parse("""
    beam {
      polarization_normal = 1,0,0
      polarization_fraction = 0.5
    }
  """)).extract()
  b3 = BeamFactory.from_phil(params3, reference)
  assert b3.get_polarization_fraction() == 0.5
  assert b3.get_polarization_normal() == (1.0, 0.0, 0.0)

def test_scan_varying():
  from scitbx import matrix
  direction = matrix.col((0.013142, 0.002200, 1.450476))
  unit_direction = direction.normalize()
  wavelength = 0.689400
  s0 = -unit_direction * 1.0 / wavelength

  # Create the beam
  b = Beam(s0)

  assert b.get_num_scan_points() == 0
  assert b.get_s0_at_scan_points().size() == 0
  with pytest.raises(RuntimeError):
    b.get_s0_at_scan_point(0) # should raise RuntimeError

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

def test_beam_object_comparison():
  from scitbx import matrix
  direction = matrix.col((0.013142, 0.002200, 1.450476))
  unit_direction = direction.normalize()
  wavelength = 0.689400
  s0 = -unit_direction * 1.0 / wavelength

  # Equal beams with scan-points set
  b1 = Beam(s0)
  b1.set_s0_at_scan_points([s0] * 5)
  b2 = Beam(s0)
  b2.set_s0_at_scan_points([s0] * 5)

  assert b1 == b2
  assert b1.is_similar_to(b2)

  # Different direction
  b3 = Beam(-s0)
  b3.set_s0_at_scan_points([-s0] * 5)
  assert b1 != b3
  assert not b1.is_similar_to(b3)

  # Different wavelength
  b4 = Beam(s0 * 1.5)
  b4.set_s0_at_scan_points([s0 * 1.5] * 5)
  assert b1 != b4
  assert not b1.is_similar_to(b4)
