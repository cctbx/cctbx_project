from __future__ import absolute_import, division, print_function

import six.moves.cPickle as pickle
from dxtbx.model import Detector, Panel
from libtbx.test_utils import approx_equal
from scitbx.array_family import flex

def tst_get_gain(detector):
  detector[0].set_gain(2.0)
  assert abs(detector[0].get_gain() - 2.0) < 1e-7

def tst_get_identifier(detector):
  detector[0].set_identifier("HELLO")
  assert detector[0].get_identifier() == "HELLO"
  detector2 = pickle.loads(pickle.dumps(detector))
  assert detector[0].get_identifier() == detector2[0].get_identifier()

def tst_get_pixel_lab_coord(detector):
  from scitbx import matrix
  eps = 1e-7

  # Check lab coordinates at the origin
  orig = detector[0].get_pixel_lab_coord((0, 0))
  dorig = abs(matrix.col(orig) - matrix.col(detector[0].get_origin()))
  assert(dorig < eps)

  # Check lab coordinate at the opposite corner
  corner = detector[0].get_pixel_lab_coord((512, 512))
  corner2 = (512 * 0.172, 512 * 0.172, 200)
  dcorner = abs(matrix.col(corner) - matrix.col(corner))
  assert(dcorner < eps)

def tst_get_image_size_mm(detector):
  from scitbx import matrix
  eps = 1e-7
  size = detector[0].get_image_size_mm()
  size2 = (512 * 0.172, 512 * 0.172)
  dsize = abs(matrix.col(size) - matrix.col(size2))
  assert(dsize < eps)

def tst_is_value_in_trusted_range(detector):
  """Check values are either inside or outside trusted range."""
  assert(detector[0].is_value_in_trusted_range(-1) == False)
  assert(detector[0].is_value_in_trusted_range(0) == True)
  assert(detector[0].is_value_in_trusted_range(999) == True)
  assert(detector[0].is_value_in_trusted_range(1000) == False)

def tst_is_coord_valid(detector):
  """Check points are either inside or outside detector range."""
  assert(detector[0].is_coord_valid((-1, 256)) == False)
  assert(detector[0].is_coord_valid((256, 256)) == True)
  assert(detector[0].is_coord_valid((512, 256)) == False)
  assert(detector[0].is_coord_valid((256, -1)) == False)
  assert(detector[0].is_coord_valid((256, 256)) == True)
  assert(detector[0].is_coord_valid((256, 513)) == False)

def tst_pixel_to_millimeter_to_pixel(detector):
  from scitbx import matrix
  from random import random
  eps = 1e-7

  # Pick some random pixels and check that px -> mm -> px give px == px
  w, h = detector[0].get_image_size()
  random_pixel = lambda: (random() * w, random() * h)
  pixels = flex.vec2_double(random_pixel() for i in range(100))
  xy_mm = detector[0].pixel_to_millimeter(pixels)
  xy_px = detector[0].millimeter_to_pixel(xy_mm)
  assert approx_equal(xy_px, pixels, eps=eps)
  for xy in pixels:
    xy_mm = detector[0].pixel_to_millimeter(xy)
    xy_px = detector[0].millimeter_to_pixel(xy_mm)
    assert(abs(matrix.col(xy_px) - matrix.col(xy)) < eps)

def tst_parallax_correction(detector):
  from random import uniform
  from scitbx import matrix
  random_coord = lambda: (
      uniform(-1000, 1000),
      uniform(-1000, 1000))
  for i in range(10000):
    mm = random_coord()
    px = detector[0].millimeter_to_pixel(mm)
    mm2 = detector[0].pixel_to_millimeter(px)
    assert(abs(matrix.col(mm) - matrix.col(mm2)) < 1e-3)

def tst_get_names(detector):
  names = detector.get_names()
  assert(len(names) == 1)
  assert(names[0] == 'Panel')

def tst_get_thickness(detector):
  for panel in detector:
    assert(panel.get_thickness() == 0.1)

def tst_get_material(detector):
  for panel in detector:
    assert(panel.get_material() == 'Si')

def tst_set_mosflm_beam_centre(detector):
  from scitbx import matrix
  from dxtbx.model import Beam
  wavelength = 1
  panel = detector[0]
  detector_normal = matrix.col(panel.get_normal())
  origin = matrix.col(panel.get_origin())
  fast_axis = matrix.col(panel.get_fast_axis())
  slow_axis = matrix.col(panel.get_slow_axis())
  image_size = panel.get_image_size_mm()

  s0 = (1.0/wavelength) * detector_normal
  beam = Beam(-s0.normalize(), wavelength)

  beam_centre = matrix.col(panel.get_beam_centre(beam.get_s0()))
  origin_shift = matrix.col((1, 0.5))
  new_beam_centre = beam_centre + origin_shift

  new_mosflm_beam_centre = tuple(reversed(new_beam_centre))

  from dxtbx.model.detector_helpers import set_mosflm_beam_centre
  set_mosflm_beam_centre(detector, beam, new_mosflm_beam_centre)

  assert (matrix.col(panel.get_beam_centre(beam.get_s0())) -
          matrix.col(tuple(reversed(new_mosflm_beam_centre)))).length() < 1e-6

def tst_detectors_are_same(detA, detB):
  '''Equality operator on detector objects must identify identical detectors'''
  assert(detA == detB)

def tst_detectors_are_different(detA, detB):
  '''Equality operator on detector objects must find differences in origin'''
  assert(detA != detB)

def tst_resolution(detector):
  from dxtbx.model import Beam
  beam = Beam(direction=(0,0,1), wavelength=1.0)
  d_min1 = detector.get_max_resolution(beam.get_s0())
  d_min2 = detector.get_max_inscribed_resolution(beam.get_s0())
  assert d_min1 < d_min2

def tst_panel_mask():
  from dxtbx.model import Panel

  panel = Panel()
  panel.set_image_size((100, 100))
  panel.add_mask(40,0,60,100)
  panel.add_mask(0,40,100,60)
  panel.set_trusted_range((-1, 10))

  data = flex.double(flex.grid(100,100))
  data[10,10] = -1
  data[20,20] = 10
  data[30,30] = 100
  data[40,40] = -10

  m1 = panel.get_untrusted_rectangle_mask()
  m2 = panel.get_trusted_range_mask(data)

  count = 0
  for j in range(100):
    for i in range(40,60):
      assert(m1[j,i] == False)
      count += 1
  for i in range(100):
    for j in range(40,60):
      if i >= 40 and i < 60:
        continue
      assert(m1[j,i] == False)
      count += 1
  assert m1.count(False) == count, "%d, %d" % (m1.count(False), count)

  assert m2.count(False) == 4
  assert m2[10,10] == False
  assert m2[20,20] == False
  assert m2[30,30] == False
  assert m2[40,40] == False

def test_detector():
  from dxtbx.model import ParallaxCorrectedPxMmStrategy

  def create_detector(offset = 0):
    # Create the detector
    detector = Detector(Panel(
      "",                 # Type
      "Panel",            # Name
      (10, 0, 0),         # Fast axis
      (0, 10, 0),         # Slow axis
      (0 + offset, 0 + offset, 200 - offset),
                          # Origin
      (0.172, 0.172),     # Pixel size
      (512, 512),         # Image size
      (0, 1000),          # Trusted range
      0.1,                # Thickness
      "Si",               # Material
      identifier="123"))  # Identifier
    return detector

  detector = create_detector()

  # Perform some tests
  tst_get_identifier(detector)
  tst_get_gain(detector)
  tst_set_mosflm_beam_centre(detector)
  tst_get_pixel_lab_coord(detector)
  tst_get_image_size_mm(detector)
  tst_is_value_in_trusted_range(detector)
  tst_is_coord_valid(detector)
  tst_pixel_to_millimeter_to_pixel(detector)
  tst_get_names(detector)
  tst_get_thickness(detector)
  tst_get_material(detector)
  tst_resolution(detector)
  tst_panel_mask()

  # Attenuation length
  from cctbx.eltbx import attenuation_coefficient
  table = attenuation_coefficient.get_table("Si")
  mu = table.mu_at_angstrom(1) / 10.0
  t0 = 0.320

  # Create another detector with different origin
  detector_moved = create_detector(offset=100)
  tst_detectors_are_different(detector, detector_moved)

  detector_moved_copy = create_detector(offset=100)
  tst_detectors_are_same(detector_moved, detector_moved_copy)

  # Create the detector
  detector = Detector(Panel(
      "",                 # Type
      "",                 # Name
      (10, 0, 0),         # Fast axis
      (0, 10, 0),         # Slow axis
      (0, 0, 200),        # Origin
      (0.172, 0.172),     # Pixel size
      (512, 512),         # Image size
      (0, 1000),          # Trusted range
      0.0,                # Thickness
      "",                 # Material
      ParallaxCorrectedPxMmStrategy(mu, t0)))

  tst_parallax_correction(detector)
