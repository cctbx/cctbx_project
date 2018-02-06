from __future__ import absolute_import, division, print_function

import libtbx.load_env

class Test(object):

  def __init__(self):
    pass

  def run(self):
    self.tst_beam()
    self.tst_beam_with_scan_points()
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
    assert 's0_at_scan_points' not in d

    # Test with a template and partial dictionary
    d2 = {'direction' : (0, 1, 0), 'divergence' : 0.2 }
    b3 = BeamFactory.from_dict(d2, d)
    assert(b3.get_direction() == (0, 1, 0))
    assert(b3.get_wavelength() == 2)
    assert(abs(b3.get_divergence() - 0.2) <= 1e-7)
    assert(abs(b3.get_sigma_divergence() - 0.1) <= 1e-7)
    assert(b2 != b3)

  def tst_beam_with_scan_points(self):
    from dxtbx.model import Beam, BeamFactory
    b1 = Beam((1, 0, 0), 2, 0.1, 0.1)
    from scitbx import matrix
    s0_static = matrix.col(b1.get_s0())
    s0_as_scan_points = [s0_static] * 5
    b1.set_s0_at_scan_points([s0_static] * 5)
    d = b1.to_dict()
    b2 = BeamFactory.from_dict(d)

    for s0comp in d['s0_at_scan_points']:
      assert matrix.col(s0comp) == s0_static

    for s0comp in b2.get_s0_at_scan_points():
      assert matrix.col(s0comp) == s0_static

    assert(b1 == b2)

  def tst_detector(self):
    pass

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

  def tst_scan(self):
    from dxtbx.model import Scan, ScanFactory
    from scitbx.array_family import flex
    s1 = Scan((1, 3), (1.0, 0.2), flex.double([0.1, 0.1, 0.1]), flex.double([0.1, 0.2, 0.3]), 0)
    d = s1.to_dict()
    s2 = ScanFactory.from_dict(d)
    assert(d['image_range'] == (1, 3))
    assert(d['oscillation'] == (1.0, 0.2))
    assert(d['exposure_time'] == [0.1, 0.1, 0.1])
    assert(d['epochs'] == [0.1, 0.2, 0.3])
    assert(d['batch_offset'] == 0)
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
    d3 = { 'image_range' : (1, 10), 'epochs' : [0.1, 0.2],}
    s4 = ScanFactory.from_dict(d3, d)
    assert(abs(s4.get_epochs()[2] - 0.3) < 1e-7)
    assert(abs(s4.get_epochs()[9] - 1.0) < 1e-7)

    d4 = {'batch_offset': 100}
    s5 = ScanFactory.from_dict(d4, d)
    assert(s5.get_batch_offset() == 100)
    assert(s5.get_batch_range() == (101, 103))

if __name__ == '__main__':
  if libtbx.env.has_module("dials") and \
     libtbx.env.has_module("dials_regression"):
    test = Test()
    test.run()
  else:
    print("Skipping test: dials or dials_regression not present")
