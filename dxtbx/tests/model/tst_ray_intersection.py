from __future__ import division

def tst_intersection_at_origin(intersection, wavelength, origin):

  eps = 1e-7

  # Set the beam vector to directed towards the detector origin
  s = (1.0 / wavelength) * origin

  # Get the x, y pixel coordinates
  x, y = intersection(s)

  # Check that they are equal to the detector origin
  assert(abs(x) < eps)
  assert(abs(y) < eps)

  # Test Passed
  print "OK"

def tst_intersection_at_corners(intersection, panel, wavelength):

  from scitbx import matrix
  eps = 1e-7

  # Get parameters from detector
  origin = matrix.col(panel.get_origin())
  fast_axis = matrix.col(panel.get_fast_axis())
  slow_axis = matrix.col(panel.get_slow_axis())
  image_size = panel.get_image_size_mm()

  # Set the beam vector to directed towards the detector origin
  s1 = (1.0 / wavelength) * origin
  s2 = (1.0 / wavelength) * (origin + fast_axis * image_size[0])
  s3 = (1.0 / wavelength) * (origin + slow_axis * image_size[1])
  s4 = (1.0 / wavelength) * (origin + fast_axis * image_size[0] +
                                      slow_axis * image_size[1])

  # Get the x, y pixel coordinates
  x1, y1 = intersection(s1)
  x2, y2 = intersection(s2)
  x3, y3 = intersection(s3)
  x4, y4 = intersection(s4)

  # Check that they are equal to the detector origin
  assert(abs(x1 - 0) < eps and abs(y1 - 0) < eps)
  assert(abs(x2 - image_size[0]) < eps and abs(y2 - 0) < eps)
  assert(abs(x3 - 0) < eps and abs(y3 - image_size[1]) < eps)
  assert(abs(x4 - image_size[0]) < eps and abs(y4 - image_size[1]) < eps)

  # Test Passed
  print "OK"

def tst_intersection_away_from_panel(intersection, panel, wavelength):

  from scitbx import matrix

  # Set the beam vector to directed away from the detector origin
  s = - (1.0 / wavelength) * matrix.col(panel.get_origin())

  # Check that they are equal to None
  try:
    intersection(s)
  except(RuntimeError):
    # Test Passed
    print "OK"
    return

  # Test Failed
  assert(False)

def tst_beam_plane_intersection():
  from dxtbx.model import Panel
  from scitbx import matrix

  # The input parameters (from a GXPARM.XDS file)
  fast_axis = (0.000000, -0.939693, -0.342020)
  slow_axis = (1.000000,  0.000000,  0.000000)
  normal = (0.000000, -0.342020,  0.939693)
  pixel_size = (0.172000, 0.172000)
  pixel_origin = (244.836136, 320.338531)
  image_size = (487, 619)
  distance = 122.124901
  wavelength = 0.689400

  # Calculate the vector to the detector (0, 0) pixel
  origin = ((matrix.col(fast_axis).normalize() * pixel_size[0]) *
              (0 - pixel_origin[0]) +
            (matrix.col(slow_axis).normalize() * pixel_size[1]) *
              (0 - pixel_origin[1]) +
            (distance * matrix.col(normal).normalize()))

  # Create a detector object
  panel = Panel("", "", fast_axis, slow_axis, origin, pixel_size,
                image_size, (0, 0), 0.0, "")

  # Create the intersection object
  intersection = lambda x: panel.get_ray_intersection(x)

  # Perform a load of tests
  tst_intersection_at_origin(intersection, wavelength, origin)
  tst_intersection_at_corners(intersection, panel, wavelength)
  tst_intersection_away_from_panel(intersection, panel, wavelength)

def tst_transform_at_origin(transform, panel):

  eps = 1e-7

  # Get lab coordinate at point
  x, y, z = transform((0, 0))

  # Check coordinate is at origin
  assert(abs(x - panel.get_origin()[0]) < eps)
  assert(abs(y - panel.get_origin()[1]) < eps)
  assert(abs(z - panel.get_origin()[2]) < eps)

  # Test Passed
  print "OK"

def tst_transform_at_corners(transform, panel):
  from scitbx import matrix
  eps = 1e-7

  # Get parameters from detector
  origin = matrix.col(panel.get_origin())
  fast_axis = matrix.col(panel.get_fast_axis())
  slow_axis = matrix.col(panel.get_slow_axis())
  image_size = panel.get_image_size_mm()

  # Set the beam vector to directed towards the detector origin
  xyz1 = origin
  xyz2 = origin + fast_axis * image_size[0]
  xyz3 = origin + slow_axis * image_size[1]
  xyz4 = origin + fast_axis * image_size[0] + slow_axis * image_size[1]

  # Get the x, y pixel coordinates
  xyz11 = matrix.col(transform((0, 0)))
  xyz22 = matrix.col(transform((image_size[0], 0)))
  xyz33 = matrix.col(transform((0, image_size[1])))
  xyz44 = matrix.col(transform((image_size[0], image_size[1])))

  # Check that they are equal to the detector origin
  assert(abs(xyz11 - xyz1) < eps)
  assert(abs(xyz22 - xyz2) < eps)
  assert(abs(xyz33 - xyz3) < eps)
  assert(abs(xyz44 - xyz4) < eps)

  # Test Passed
  print "OK"

def tst_plane_to_lab_transform():
  from dxtbx.model import Panel
  from scitbx import matrix

  # The input parameters (from a GXPARM.XDS file)
  fast_axis = (0.000000, -0.939693, -0.342020)
  slow_axis = (1.000000,  0.000000,  0.000000)
  normal = (0.000000, -0.342020,  0.939693)
  pixel_size = (0.172000, 0.172000)
  pixel_origin = (244.836136, 320.338531)
  image_size = (487, 619)
  distance = 122.124901
  wavelength = 0.689400

  # Calculate the vector to the detector (0, 0) pixel
  origin = ((matrix.col(fast_axis).normalize() * pixel_size[0]) *
              (0 - pixel_origin[0]) +
            (matrix.col(slow_axis).normalize() * pixel_size[1]) *
              (0 - pixel_origin[1]) +
            (distance * matrix.col(normal).normalize()))

  # Create a detector object
  panel = Panel("", "", fast_axis, slow_axis, origin, pixel_size,
               image_size, (0, 0), 0.0, "")

  # Create the intersection object
  transform = lambda x: panel.get_lab_coord(x)

  # Perform a load of tests
  tst_transform_at_origin(transform, panel)
  tst_transform_at_corners(transform, panel)

def tst_fwd_rev_random(intersection, transform, panel):
  from scitbx import matrix
  from random import random

  eps = 1e-7

  # Set the shift parameters
  image_size = panel.get_image_size_mm()
  random_xy = lambda: (random() * image_size[0], random() * image_size[1])

  # Loop a number of times
  num = 1000
  for i in range(num):

    # Create a detector coordinate
    xy = matrix.col(random_xy())

    # Calculate the lab coordinate of the vector
    s = transform(xy)

    # Calculate the detector coordinate
    xy_2 = intersection(s)

    # Check the vectors are almost equal
    assert(abs(xy - matrix.col(xy_2)) < eps)

  print "OK"

def tst_forward_and_reverse_transform():
  from dxtbx.model import Panel
  from scitbx import matrix

  # The input parameters (from a GXPARM.XDS file)
  fast_axis = (0.000000, -0.939693, -0.342020)
  slow_axis = (1.000000,  0.000000,  0.000000)
  normal = (0.000000, -0.342020,  0.939693)
  pixel_size = (0.172000, 0.172000)
  pixel_origin = (244.836136, 320.338531)
  image_size = (487, 619)
  distance = 122.124901
  wavelength = 0.689400

  # Calculate the vector to the detector (0, 0) pixel
  origin = ((matrix.col(fast_axis).normalize() * pixel_size[0]) *
              (0 - pixel_origin[0]) +
            (matrix.col(slow_axis).normalize() * pixel_size[1]) *
              (0 - pixel_origin[1]) +
            (distance * matrix.col(normal).normalize()))

  # Create a detector object
  panel = Panel("", "", fast_axis, slow_axis, origin, pixel_size,
               image_size, (0, 0), 0.0, "")

  # Create the intersection object and transform object
  intersection = lambda x: panel.get_ray_intersection(x)
  transform = lambda x: panel.get_lab_coord(x)

  # Do a test
  tst_fwd_rev_random(intersection, transform, panel)

def run():
  tst_beam_plane_intersection()
  tst_plane_to_lab_transform()
  tst_forward_and_reverse_transform()

if __name__ == '__main__':
  run()
