from __future__ import absolute_import, division

import libtbx.load_env

class Test(object):

  def __init__(self):
    import os
    dials_regression = libtbx.env.dist_path('dials_regression')
    filename = os.path.join(dials_regression, 'image_examples', 'XDS',
        'XPARM.XDS')

    import dxtbx
    models = dxtbx.load(filename)
    self.detector = models.get_detector()
    assert(len(self.detector) == 1)
    self.attlen = 0.252500934883
    self.distance = self.detector[0].get_distance()
    self.origin = self.detector[0].get_ray_intersection(
        self.detector[0].get_normal())

  def run(self):
    from random import uniform

    # Generate some random coordinates and do the correction
    random_coord = lambda: (
        uniform(-1000, 1000),
        uniform(-1000, 1000))
    for i in range(10000):
      xy = random_coord()
      self.tst_single(xy)

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
    return mmcal

  def correct(self, xy):
    from dxtbx.model import parallax_correction
    return parallax_correction(self.distance, self.attlen, self.origin, xy)

  def correct_inv(self, xy):
    from dxtbx.model import parallax_correction_inv
    return parallax_correction_inv(self.distance, self.attlen,
        self.origin, xy)


if __name__ == '__main__':
  if not libtbx.env.has_module("dials"):
    print "Skipping test: dials not present"
  elif not libtbx.env.has_module("dials_regression"):
    print "Skipping test: dials_regression not present"
  else:
    test = Test()
    test.run()
