from __future__ import absolute_import, division, print_function

import os

def test_multi_axis_goniometer(dials_regression):
  from libtbx.test_utils import approx_equal

  data_dir = os.path.join(dials_regression, 'image_examples', 'dials-190')

  from dxtbx.imageset import ImageSetFactory
  imgset = ImageSetFactory.new(os.path.join(data_dir, "whatev1_01_00001.cbf"))[0]
  gonio = imgset.get_goniometer(0)
  assert approx_equal(gonio.get_fixed_rotation(), (1,0,0,0,1,0,0,0,1))
  assert approx_equal(gonio.get_setting_rotation(), (1,0,0,0,1,0,0,0,1))

  from dxtbx.imageset import ImageSetFactory
  imgset = ImageSetFactory.new(os.path.join(data_dir, "whatev1_02_00001.cbf"))[0]
  gonio = imgset.get_goniometer(0)
  assert approx_equal(gonio.get_fixed_rotation(), (1,0,0,0,1,0,0,0,1))
  assert approx_equal(
    gonio.get_setting_rotation(),
    (1, 0, 0, 0, -0.5, 0.866, 0.0, -0.866, -0.5), eps=1e-4)

  from dxtbx.imageset import ImageSetFactory
  imgset = ImageSetFactory.new(os.path.join(data_dir, "whatev1_03_00001.cbf"))[0]
  gonio = imgset.get_goniometer(0)
  assert approx_equal(gonio.get_fixed_rotation(), (1,0,0,0,1,0,0,0,1))
  assert approx_equal(
    gonio.get_setting_rotation(),
    (1, 0, 0, 0, 0.5, 0.866, 0.0, -0.866, 0.5), eps=1e-4)

def test_still(dials_regression):
  from libtbx.test_utils import approx_equal

  data_dir = os.path.join(dials_regression, 'image_examples', 'DLS_I04')

  from dxtbx.imageset import ImageSetFactory
  imgset = ImageSetFactory.new(os.path.join(data_dir, "grid_full_cbf_0005.cbf"))[0]
  if imgset.get_scan(0):
    scan = imgset.get_scan(0)
    assert approx_equal(scan.get_image_oscillation(0), (30.0, 0.0), eps=1e-5)
  beam = imgset.get_beam(0)
  beam.get_s0()
  assert approx_equal(beam.get_s0(),
                      (-0.0, -0.0, -1.0209290454313424))
