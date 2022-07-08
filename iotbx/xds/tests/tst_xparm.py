from __future__ import absolute_import, division, print_function

class Test(object):

  def __init__(self):
    pass

  def run(self):

    from iotbx.xds import xparm
    import os
    import libtbx.load_env
    from libtbx.test_utils import open_tmp_file

    iotbx_dir = libtbx.env.dist_path('iotbx')
    filename = os.path.join(iotbx_dir, 'xds', 'tests', 'XPARM.XDS')
    handle = xparm.reader()
    assert handle.find_version(filename) == 1
    handle.read_file(filename)

    print('OK')

    filename = os.path.join(iotbx_dir, 'xds', 'tests', 'NEW_XPARM.XDS')
    handle = xparm.reader()
    assert handle.find_version(filename) == 2
    handle.read_file(filename)
    print('OK')

    f = open_tmp_file(suffix='XPARM.XDS', mode='wb')
    f.close()
    xds_str = xparm.write(
      handle.starting_frame,
      handle.starting_angle,
      handle.oscillation_range,
      handle.rotation_axis,
      handle.wavelength,
      handle.beam_vector,
      handle.space_group,
      handle.unit_cell,
      handle.unit_cell_a_axis,
      handle.unit_cell_b_axis,
      handle.unit_cell_c_axis,
      handle.num_segments,
      handle.detector_size,
      handle.pixel_size,
      handle.detector_origin,
      handle.detector_distance,
      handle.detector_x_axis,
      handle.detector_y_axis,
      handle.detector_normal,
      handle.segments,
      handle.orientation)
    with open(f.name, 'w') as f:
      f.write(xds_str)
    handle_recycled = xparm.reader()
    # make sure we wrote out version 2
    assert handle_recycled.find_version(f.name) == 2
    handle_recycled.read_file(f.name)

    for handle in (handle, handle_recycled):

      # Scan and goniometer stuff
      assert handle.starting_frame == 1
      assert handle.starting_angle == 82.0
      assert handle.oscillation_range == 0.1500
      assert handle.rotation_axis == (0.999997, -0.001590, -0.001580)

      # Beam stuff
      assert handle.wavelength == 0.976250
      assert handle.beam_vector == (0.001608, 0.004392, 1.024317)

      # Detector stuff
      assert handle.detector_size == (2463, 2527)
      assert handle.pixel_size == (0.172, 0.172)

      assert handle.detector_distance == 264.928955
      assert handle.detector_origin == (1224.856812, 1187.870972)
      assert handle.detector_x_axis == (1.0, 0.0, 0.0)
      assert handle.detector_y_axis == (0.0, 1.0, 0.0)
      assert handle.detector_normal == (0.0, 0.0, 1.0)

      # Crystal stuff
      assert handle.space_group == 75
      assert handle.unit_cell == (57.7831, 57.7831, 150.0135, 90.000, 90.000, 90.000)
      assert handle.unit_cell_a_axis == (-14.918090, -22.358297, 51.151196)
      assert handle.unit_cell_b_axis == (-19.858326, 51.608330, 16.766487)
      assert handle.unit_cell_c_axis == (-135.447952, -34.400188, -54.539391)

    # segment stuff
    assert handle_recycled.num_segments == 1
    assert handle_recycled.segments == [(1, 1, 2463, 1, 2527)]
    assert handle_recycled.orientation == [
      (0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0)]

    print('OK')

if __name__ == '__main__':

  test = Test()
  test.run()
