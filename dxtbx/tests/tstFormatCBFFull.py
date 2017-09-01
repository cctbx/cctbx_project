from __future__ import absolute_import, division
def exercise_multi_axis_goniometer():
  import libtbx.load_env
  from libtbx.test_utils import approx_equal
  import os

  if not libtbx.env.has_module("dials"):
    print "Skipping tstFormatCBFFull.py: dials not present"
    return
  if not libtbx.env.has_module("dials_regression"):
    print "Skipping tstFormatCBFFull.py: dials_regression not present"
    return

  data_dir = libtbx.env.find_in_repositories(
    relative_path="dials_regression/image_examples/dials-190",
    test=os.path.isdir)

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

def exercise_still():
  import libtbx.load_env
  from libtbx.test_utils import approx_equal
  import os

  if not libtbx.env.has_module("dials"):
    print "Skipping tstFormatCBFFull.py: dials not present"
    return
  if not libtbx.env.has_module("dials_regression"):
    print "Skipping tstFormatCBFFull.py: dials_regression not present"
    return

  data_dir = libtbx.env.find_in_repositories(
    relative_path="dials_regression/image_examples/DLS_I04",
    test=os.path.isdir)

  from dxtbx.imageset import ImageSetFactory
  imgset = ImageSetFactory.new(os.path.join(data_dir, "grid_full_cbf_0005.cbf"))[0]
  if imgset.get_scan(0):
    scan = imgset.get_scan(0)
    assert approx_equal(scan.get_image_oscillation(0), (30.0, 0.0), eps=1e-5)
  beam = imgset.get_beam(0)
  beam.get_s0()
  assert approx_equal(beam.get_s0(),
                      (-0.0, -0.0, -1.0209290454313424))

def run():
  exercise_still()
  exercise_multi_axis_goniometer()
  print "OK"

if __name__ == '__main__':
  run()
