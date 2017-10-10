from __future__ import absolute_import, division
from __future__ import print_function
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
  print("OK")

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
  print("OK")

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

  print('OK')

def run():
  """Test the beam object"""
  tst_set_direction_wavelength()
  tst_set_s0()
  tst_from_phil()

if __name__ == '__main__':
  run()
