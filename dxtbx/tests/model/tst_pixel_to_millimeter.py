
from __future__ import division

class Test(object):

  def __init__(self):
    import os
    import libtbx.load_env
    try:
      dials_regression = libtbx.env.dist_path('dials_regression')
    except KeyError, e:
      print 'FAIL: dials_regression not configured'
      exit(0)

    filename = os.path.join(dials_regression, 'image_examples', 'XDS',
        'XPARM.XDS')

    import dxtbx
    models = dxtbx.load(filename)
    self.detector = models.get_detector()
    self.beam = models.get_beam()
    assert(len(self.detector) == 1)
    self.t0 = 0.320
    from cctbx.eltbx import attenuation_coefficient
    table = attenuation_coefficient.get_table("Si")
    self.mu = table.mu_at_angstrom(self.beam.get_wavelength()) / 10
    self.distance = self.detector[0].get_distance()
    self.origin = self.detector[0].get_ray_intersection(
        self.detector[0].get_normal())[1]
    self.pixel_size = self.detector[0].get_pixel_size()

  def run(self):
    from random import uniform

    # Generate some random coordinates and do the correction
    random_coord = lambda: (
        uniform(-1000, 1000),
        uniform(-1000, 1000))
    for i in range(10000):
      xy = random_coord()
      self.tst_single(xy)

    from scitbx.array_family import flex
    xy = flex.vec2_double([random_coord() for i in range(100)])
    self.tst_array(xy)

    self.tst_inverted_axis()

    # Test passed
    print 'OK'

  def tst_single(self, xy):

    from scitbx import matrix
    xy = matrix.col(xy)

    # Do the forward and reverse corrections
    xy_corr_gold = matrix.col(self.correct_gold(xy))
    xy_corr      = matrix.col(self.correct(xy))
    xy_corr_inv  = matrix.col(self.correct_inv(xy_corr))

    # Check the values
    assert(abs(xy_corr_gold - xy_corr) < 1e-7)
    assert(abs(xy_corr_inv - xy) < 1e-3)

  def tst_array(self, xy):
    from libtbx.test_utils import approx_equal
    xy_corr = self.detector[0].get_lab_coord(xy)
    xy_corr_panel = self.detector[0].get_lab_coord(xy)
    xy_corr_gold = [self.detector[0].get_lab_coord(xy_single)
                    for xy_single in xy]
    assert approx_equal(xy_corr, xy_corr_gold)
    assert approx_equal(xy_corr_panel, xy_corr_gold)

    print 'OK'

  def tst_inverted_axis(self):

    def get_values(invert_y):

      from dxtbx.model.beam import beam_factory
      beam = beam_factory.simple(wavelength = 1)

      if invert_y:
        y_direction = "-y"
      else:
        y_direction = "+y"

      from dxtbx.model.detector import detector_factory
      detector = detector_factory.simple(
        sensor = detector_factory.sensor("PAD"),
        distance = 100,
        beam_centre = [50,50],
        fast_direction = "+x",
        slow_direction = y_direction,
        pixel_size = [0.1,0.1],
        image_size = [1000,1000]
        )


      from dxtbx.model import ParallaxCorrectedPxMmStrategy
      from cctbx.eltbx import attenuation_coefficient
      wavelength = beam.get_wavelength()
      thickness = 0.5
      table = attenuation_coefficient.get_table("Si")
      mu = table.mu_at_angstrom(wavelength) / 10.0
      t0 = thickness

      for panel in detector:
        panel.set_px_mm_strategy(ParallaxCorrectedPxMmStrategy(mu, t0))
      v1 = detector[0].pixel_to_millimeter((0,0))
      v2 = detector[0].pixel_to_millimeter((1000,1000))

      return v1, v2

    v11, v12 = get_values(False)
    v21, v22 = get_values(False)

    assert abs(v11[0] - v21[0]) < 1e-7
    assert abs(v11[1] - v21[1]) < 1e-7
    assert abs(v12[0] - v22[0]) < 1e-7
    assert abs(v12[1] - v22[1]) < 1e-7

    print 'OK'

  def correct_gold(self, xy):
    from scitbx import matrix
    from math import exp
    mu = self.mu
    t0 = self.t0
    s1 = matrix.col(self.detector[0].get_lab_coord(xy)).normalize()
    d0 = matrix.col(self.detector[0].get_origin())
    d1 = matrix.col(self.detector[0].get_fast_axis())
    d2 = matrix.col(self.detector[0].get_slow_axis())
    dn = d1.cross(d2)
    cos_theta = s1.dot(dn)
    t = t0 / cos_theta
    o = (1.0 / mu) - (t + 1.0 / mu) * exp(-mu * t)
    cx = xy[0] + s1.dot(d1) * o
    cy = xy[1] + s1.dot(d2) * o
    return (cx / self.pixel_size[0], cy / self.pixel_size[1])

  def correct(self, xy):
    from dxtbx.model import ParallaxCorrectedPxMmStrategy
    convert = ParallaxCorrectedPxMmStrategy(self.mu, self.t0)
    return convert.to_pixel(self.detector[0], xy)

  def correct_inv(self, xy):
    from dxtbx.model import ParallaxCorrectedPxMmStrategy
    convert = ParallaxCorrectedPxMmStrategy(self.mu, self.t0)
    return convert.to_millimeter(self.detector[0], xy)


if __name__ == '__main__':

  test = Test()
  test.run()
