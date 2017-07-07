from __future__ import absolute_import, division
import libtbx.load_env
import os

class Test(object):

  def __init__(self):
    pass

  def run(self):
    self.tst_beam()
    self.tst_detector()
    self.tst_goniometer()
    self.tst_scan()
    # self.tst_sweep()
    # self.tst_cspad_hierarchy()
    # self.tst_no_image_data()

  def tst_beam(self):
    from dxtbx.model import Beam, BeamFactory
    b1 = Beam((1, 0, 0), 2, 0.1, 0.1)
    d = b1.to_dict()
    b2 = BeamFactory.from_dict(d)
    assert(d['direction'] == (1, 0, 0))
    assert(d['wavelength'] == 2)
    assert(abs(d['divergence'] - 0.1) <= 1e-7)
    assert(abs(d['sigma_divergence'] - 0.1) <= 1e-7)
    assert(b1 == b2)

    # Test with a template and partial dictionary
    d2 = {'direction' : (0, 1, 0), 'divergence' : 0.2 }
    b3 = BeamFactory.from_dict(d2, d)
    assert(b3.get_direction() == (0, 1, 0))
    assert(b3.get_wavelength() == 2)
    assert(abs(b3.get_divergence() - 0.2) <= 1e-7)
    assert(abs(b3.get_sigma_divergence() - 0.1) <= 1e-7)
    assert(b2 != b3)
    print 'OK'

  def tst_detector(self):
    print 'OK'

  def tst_goniometer(self):
    from dxtbx.model import Goniometer, GoniometerFactory
    g1 = Goniometer((1, 0, 0), (1, 0, 0, 0, 1, 0, 0, 0, 1))
    d = g1.to_dict()
    g2 = GoniometerFactory.from_dict(d)
    assert(d['rotation_axis'] == (1, 0, 0))
    assert(d['fixed_rotation'] == (1, 0, 0, 0, 1, 0, 0, 0, 1))
    assert(g1 == g2)

    # Test with a template and partial dictionary
    d2 = { 'rotation_axis' : (0, 1, 0) }
    g3 = GoniometerFactory.from_dict(d2, d)
    assert(g3.get_rotation_axis() == (0, 1, 0))
    assert(g3.get_fixed_rotation() == (1, 0, 0, 0, 1, 0, 0, 0, 1))
    assert(g2 != g3)

    print 'OK'

  def tst_scan(self):
    from dxtbx.model import Scan, ScanFactory
    from scitbx.array_family import flex
    s1 = Scan((1, 3), (1.0, 0.2), flex.double([0.1, 0.1, 0.1]), flex.double([0.1, 0.2, 0.3]))
    d = s1.to_dict()
    s2 = ScanFactory.from_dict(d)
    assert(d['image_range'] == (1, 3))
    assert(d['oscillation'] == (1.0, 0.2))
    assert(d['exposure_time'] == [0.1, 0.1, 0.1])
    assert(d['epochs'] == [0.1, 0.2, 0.3])
    assert(s1 == s2)

    # Test with a template and partial dictionary
    d2 = { 'exposure_time' : [0.2, 0.2, 0.2] }
    s3 = ScanFactory.from_dict(d2, d)
    assert(s3.get_image_range() == (1, 3))
    assert(s3.get_oscillation() == (1.0, 0.2))
    assert(list(s3.get_exposure_times()) == [0.2, 0.2, 0.2])
    assert(list(s3.get_epochs()) == [0.1, 0.2, 0.3])
    assert(s2 != s3)

    # Test with a partial epoch
    d3 = { 'image_range' : (1, 10), 'epochs' : [0.1, 0.2] }
    s4 = ScanFactory.from_dict(d3, d)
    assert(abs(s4.get_epochs()[2] - 0.3) < 1e-7)
    assert(abs(s4.get_epochs()[9] - 1.0) < 1e-7)

    print 'OK'

  # def tst_sweep(self):
  #   from dxtbx.serialize import load
  #   from dxtbx.serialize.helpers import tuple_almost_equal
  #   path = libtbx.env.dist_path('dials_regression')

  #   filename = os.path.join(path, "centroid_test_data", "test_sweep.json")
  #   sweep = load.imageset(filename)
  #   b = sweep.get_beam()
  #   d = sweep.get_detector()
  #   g = sweep.get_goniometer()
  #   s = sweep.get_scan()
  #   eps = 1e-7
  #   assert(len(d) == 1)
  #   assert(tuple_almost_equal(b.get_direction(), (0.0, 0.0, 1.0)))
  #   assert(abs(b.get_wavelength() - 0.9795) < eps)
  #   assert(abs(b.get_divergence() - 0.0) < eps)
  #   assert(abs(b.get_sigma_divergence() - 0.058) < eps)
  #   assert(d[0].get_type() == "SENSOR_PAD")
  #   assert(d[0].get_name() == "Panel")
  #   assert(tuple_almost_equal(d[0].get_fast_axis(), (1.0, 0.0, 0.0)))
  #   assert(tuple_almost_equal(d[0].get_slow_axis(), (0.0, -1.0, 0.0)))
  #   assert(tuple_almost_equal(d[0].get_origin(), (-211.5, 219.5, -192.7)))
  #   assert(tuple_almost_equal(d[0].get_pixel_size(), (0.172, 0.172)))
  #   assert(tuple_almost_equal(d[0].get_image_size(), (2463, 2527)))
  #   assert(tuple_almost_equal(d[0].get_trusted_range(), (-1.0, 495976.0)))
  #   assert(tuple_almost_equal(g.get_rotation_axis(), (1.0, 0.0, 0.0)))
  #   assert(tuple_almost_equal(g.get_fixed_rotation(),
  #       (1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0)))
  #   assert(tuple_almost_equal(s.get_image_range(), (1, 3)))
  #   assert(tuple_almost_equal(s.get_oscillation(), (0.0, 0.2)))
  #   assert(abs(s.get_exposure_times()[0] - 0.2) < eps)

  #   print 'OK'

  # def tst_cspad_hierarchy(self):
  #   from dxtbx.serialize import dump, load
  #   from dxtbx.imageset import ImageSetFactory
  #   import os
  #   from glob import glob
  #   path = libtbx.env.dist_path('dials_regression')

  #   # Get the imageset
  #   filename = os.path.join(path, "spotfinding_test_data", "idx*.cbf")
  #   imageset = ImageSetFactory.new(glob(filename))
  #   assert(len(imageset) == 1)
  #   imageset = imageset[0]

  #   # Dump and reload
  #   from uuid import uuid4
  #   filename = '%s.json' % uuid4().hex
  #   dump.imageset(imageset, filename)
  #   imageset2 = load.imageset(filename)

  #   # Check they're are the same
  #   assert(imageset2.get_beam() == imageset.get_beam())
  #   d1 = imageset.get_detector()
  #   d2 = imageset2.get_detector()
  #   assert(len(d1) == len(d2))
  #   for i, (p1, p2) in enumerate(zip(d1, d2)):
  #     assert(p1 == p2)
  #   assert(imageset2.get_detector() == imageset.get_detector())
  #   assert(imageset2 == imageset)

  #   # Test passed
  #   print 'OK'

  # def tst_no_image_data(self):
  #   from dxtbx.serialize import load
  #   from dxtbx.imageset import ImageSweep, NullReader
  #   path = libtbx.env.dist_path('dials_regression')
  #   filename = os.path.join(
  #     path,
  #     'centroid_test_data',
  #     'sweep_no_image_data.json')

  #   imageset = load.imageset(filename)
  #   assert(isinstance(imageset, ImageSweep))
  #   assert(isinstance(imageset.reader(), NullReader))
  #   try:
  #     imageset[0]
  #     assert(False)
  #   except Exception:
  #     pass

  #   print 'OK'

if __name__ == '__main__':
  if libtbx.env.has_module("dials") and \
     libtbx.env.has_module("dials_regression"):
    test = Test()
    test.run()
  else:
    print "Skipping test: dials or dials_regression not present"
