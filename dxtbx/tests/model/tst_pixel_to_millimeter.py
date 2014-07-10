
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
    self.t0 = 0.0
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
