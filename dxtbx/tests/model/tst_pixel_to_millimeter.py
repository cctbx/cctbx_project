
from __future__ import division

class Test(object):

    def __init__(self):
        import os
        import libtbx.load_env
        try:
            dials_regression = libtbx.env.dist_path('dials_regression')
        except KeyError, e:
            print 'FAIL: dials_regression not configured'
            return

        filename = os.path.join(dials_regression, 'image_examples', 'XDS',
            'XPARM.XDS')

        import dxtbx
        models = dxtbx.load(filename)
        self.detector = models.get_detector()
        self.attlen = 0.252500934883
        self.distance = self.detector.get_distance()
        self.origin = self.detector.get_ray_intersection(
            self.detector.get_normal())[1]
        self.pixel_size = self.detector.get_pixel_size()

    def run(self):
        from random import uniform

        # Generate some random coordinates and do the correction
        random_coord = lambda: (
            uniform(-1000, 1000),
            uniform(-1000, 1000))
        for i in range(10000):
            xy = random_coord()
            self.tst_single(xy)

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

    def correct_gold(self, xy):
        from scitbx import matrix
        s1 = matrix.col(self.detector.get_lab_coord(xy))
        lab = self.attlen * s1 / s1.length()
        d = self.detector.get_d_matrix()
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
