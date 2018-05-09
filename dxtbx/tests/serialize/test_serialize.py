from __future__ import absolute_import, division, print_function

import pytest

def test_beam():
  from dxtbx.model import Beam, BeamFactory
  b1 = Beam((1, 0, 0), 2, 0.1, 0.1)
  d = b1.to_dict()
  b2 = BeamFactory.from_dict(d)
  assert d['direction'] == (1, 0, 0)
  assert d['wavelength'] == 2
  assert d['divergence'] == pytest.approx(0.1)
  assert d['sigma_divergence'] == pytest.approx(0.1)
  assert b1 == b2
  assert 's0_at_scan_points' not in d

  # Test with a template and partial dictionary
  d2 = {'direction' : (0, 1, 0), 'divergence' : 0.2 }
  b3 = BeamFactory.from_dict(d2, d)
  assert b3.get_direction() == (0, 1, 0)
  assert b3.get_wavelength() == 2
  assert b3.get_divergence() == pytest.approx(0.2)
  assert b3.get_sigma_divergence() == pytest.approx(0.1)
  assert b2 != b3

def test_beam_with_scan_points():
  from dxtbx.model import Beam, BeamFactory
  b1 = Beam((1, 0, 0), 2, 0.1, 0.1)
  from scitbx import matrix
  s0_static = matrix.col(b1.get_s0())
  b1.set_s0_at_scan_points([s0_static] * 5)
  d = b1.to_dict()
  b2 = BeamFactory.from_dict(d)

  for s0comp in d['s0_at_scan_points']:
    assert matrix.col(s0comp) == s0_static

  for s0comp in b2.get_s0_at_scan_points():
    assert matrix.col(s0comp) == s0_static

  assert b1 == b2

def test_goniometer():
  from dxtbx.model import Goniometer, GoniometerFactory
  g1 = Goniometer((1, 0, 0), (1, 0, 0, 0, 1, 0, 0, 0, 1))
  d = g1.to_dict()
  g2 = GoniometerFactory.from_dict(d)
  assert d['rotation_axis'] == (1, 0, 0)
  assert d['fixed_rotation'] == (1, 0, 0, 0, 1, 0, 0, 0, 1)
  assert g1 == g2
  assert 'setting_rotation_at_scan_points' not in d

  # Test with a template and partial dictionary
  d2 = { 'rotation_axis' : (0, 1, 0) }
  g3 = GoniometerFactory.from_dict(d2, d)
  assert g3.get_rotation_axis() == (0, 1, 0)
  assert g3.get_fixed_rotation() == (1, 0, 0, 0, 1, 0, 0, 0, 1)
  assert g2 != g3

def test_multi_axis_goniometer():
  from dxtbx.model import GoniometerFactory
  from scitbx.array_family import flex
  g1 = GoniometerFactory.multi_axis(flex.vec3_double(((1,0,0),)),
      flex.double((0,)), flex.std_string(('PHI',)), 0)
  d = g1.to_dict()
  g2 = GoniometerFactory.from_dict(d)
  assert d['axes'] == [(1, 0, 0)]
  assert d['angles'] == [0.0]
  assert d['names'] == ['PHI']
  assert d['scan_axis'] == 0
  assert g1 == g2
  assert 'setting_rotation_at_scan_points' not in d

  # Test with a template and partial dictionary
  d2 = { 'axes' : [(0, 1, 0)] }
  g3 = GoniometerFactory.from_dict(d2, d)
  assert g3.get_rotation_axis() == (0, 1, 0)
  assert g3.get_fixed_rotation() == (1, 0, 0, 0, 1, 0, 0, 0, 1)
  assert g2 != g3

def test_goniometer_with_scan_points():
  from dxtbx.model import Goniometer, GoniometerFactory
  from scitbx.array_family import flex
  simple_g = Goniometer((1, 0, 0), (1, 0, 0, 0, 1, 0, 0, 0, 1))
  multi_ax_g = GoniometerFactory.multi_axis(flex.vec3_double(((1,0,0),)),
      flex.double((0,)), flex.std_string(('PHI',)), 0)

  from scitbx import matrix
  for g1 in [simple_g, multi_ax_g]:

    S_static = matrix.sqr(g1.get_setting_rotation())
    g1.set_setting_rotation_at_scan_points([S_static] * 5)
    d = g1.to_dict()
    g2 = GoniometerFactory.from_dict(d)

    for Scomp in d['setting_rotation_at_scan_points']:
      assert matrix.sqr(Scomp) == S_static

    for Scomp in g2.get_setting_rotation_at_scan_points():
      assert matrix.sqr(Scomp) == S_static

    assert g1 == g2

def test_scan():
  from dxtbx.model import Scan, ScanFactory
  from scitbx.array_family import flex
  s1 = Scan((1, 3), (1.0, 0.2), flex.double([0.1, 0.1, 0.1]), flex.double([0.1, 0.2, 0.3]), 0)
  d = s1.to_dict()
  s2 = ScanFactory.from_dict(d)
  assert d['image_range'] == (1, 3)
  assert d['oscillation'] == (1.0, 0.2)
  assert d['exposure_time'] == [0.1, 0.1, 0.1]
  assert d['epochs'] == [0.1, 0.2, 0.3]
  assert d['batch_offset'] == 0
  assert s1 == s2

  # Test with a template and partial dictionary
  d2 = { 'exposure_time' : [0.2, 0.2, 0.2] }
  s3 = ScanFactory.from_dict(d2, d)
  assert s3.get_image_range() == (1, 3)
  assert s3.get_oscillation() == (1.0, 0.2)
  assert list(s3.get_exposure_times()) == [0.2, 0.2, 0.2]
  assert list(s3.get_epochs()) == [0.1, 0.2, 0.3]
  assert s2 != s3

  # Test with a partial epoch
  d3 = { 'image_range' : (1, 10), 'epochs' : [0.1, 0.2],}
  s4 = ScanFactory.from_dict(d3, d)
  assert s4.get_epochs()[2] == pytest.approx(0.3)
  assert s4.get_epochs()[9] == pytest.approx(1.0)

  d4 = {'batch_offset': 100}
  s5 = ScanFactory.from_dict(d4, d)
  assert s5.get_batch_offset() == 100
  assert s5.get_batch_range() == (101, 103)
