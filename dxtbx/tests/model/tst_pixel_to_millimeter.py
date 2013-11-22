
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
    assert(len(self.detector) == 1)
    self.attlen = 0.252500934883
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
    s1 = matrix.col(self.detector[0].get_lab_coord(xy))
    lab = self.attlen * s1 / s1.length()
    d = self.detector[0].get_d_matrix()
    d0 = matrix.col((d[0], d[3], d[6]))
    d1 = matrix.col((d[1], d[4], d[7]))
    mm0 = d0.dot(lab) / d0.length()
    mm1 = d1.dot(lab) / d1.length()
    mmcal = matrix.col((xy[0] + mm0, xy[1] + mm1))
    return (mmcal[0] / self.pixel_size[0], mmcal[1] / self.pixel_size[1])

  def correct(self, xy):
    from dxtbx.model import ParallaxCorrectedPxMmStrategy
    convert = ParallaxCorrectedPxMmStrategy(self.attlen)
    return convert.to_pixel(self.detector[0], xy)

  def correct_inv(self, xy):
    from dxtbx.model import ParallaxCorrectedPxMmStrategy
    convert = ParallaxCorrectedPxMmStrategy(self.attlen)
    return convert.to_millimeter(self.detector[0], xy)


if __name__ == '__main__':

  test = Test()
  test.run()
