from __future__ import absolute_import, division, print_function
from libtbx.test_utils import open_tmp_directory

class Test(object):

  def __init__(self):
    pass

  def run(self):

    from iotbx.xds import xds_inp
    import os
    import libtbx.load_env
    iotbx_dir = libtbx.env.dist_path('iotbx')

    filename = os.path.join(iotbx_dir, 'xds', 'tests', 'XDS.INP')
    handle = xds_inp.reader()
    handle.read_file(filename)

    assert handle.detector == 'PILATUS'
    assert handle.minimum_valid_pixel_value == 0
    assert handle.overload == 1048500
    assert handle.corrections == 'ALL'
    assert handle.direction_of_detector_x_axis == [1.0, 0.0, 0.0]
    assert handle.direction_of_detector_y_axis == [0.0, 1.0, 0.0]
    assert handle.trusted_region == [0.0, 1.41]
    assert handle.sensor_thickness == 0.32
    assert handle.untrusted_rectangle == [
      [487, 495, 0, 2528], [981, 989, 0, 2528], [1475, 1483, 0, 2528],
      [1969, 1977, 0, 2528], [0, 2464, 195, 213], [0, 2464, 407, 425],
      [0, 2464, 619, 637], [0, 2464, 831, 849], [0, 2464, 1043, 1061],
      [0, 2464, 1255, 1273], [0, 2464, 1467, 1485], [0, 2464, 1679, 1697],
      [0, 2464, 1891, 1909], [0, 2464, 2103, 2121], [0, 2464, 2315, 2333]]
    assert handle.maximum_number_of_processor == 16
    assert (handle.nx, handle.ny, handle.px, handle.py) == (2463, 2527, 0.172, 0.172)
    assert (handle.orgx, handle.orgy) == (1279.1, 1235.3)
    assert handle.rotation_axis == [1.0, 0.0, 0.0]
    assert handle.detector_distance == 190.18
    assert handle.xray_wavelength == 0.9795
    assert handle.incident_beam_direction == [0.0, 0.0, 1.0]
    assert handle.fraction_of_polarization == 0.99
    assert handle.polarization_plane_normal == [0.0, 1.0, 0.0]
    assert handle.friedels_law
    assert len(handle.name_template_of_data_frames) == 1
    assert handle.name_template_of_data_frames[0].endswith(
      'X4_lots_M1S4_1_????.cbf')
    assert handle.starting_angle == 0
    assert handle.starting_frame == 1
    assert handle.include_resolution_range == [30.0, 1.27]
    assert handle.unit_cell_constants == [42.45, 42.45, 39.8, 90.0, 90.0, 90.0]
    assert handle.space_group_number == 89
    assert handle.max_fac_rmeas == 3.0
    assert handle.data_range == [1,900]

    filename_1 = os.path.join(iotbx_dir, 'xds', 'tests', 'XDS_2.INP')
    tmp_dir = open_tmp_directory()
    filename_2 = os.path.join(tmp_dir, 'XDS.INP')
    with open(filename_2, 'w') as f2, open(filename_1, 'r') as f1:
      print(f1.read(), file=f2)
    handle = xds_inp.reader()
    handle.read_file(filename_2)
    assert handle.corrections is None
    assert handle.incident_beam_direction == [-0.003, 0.001, 1.032]
    assert len(handle.name_template_of_data_frames) == 2
    assert handle.name_template_of_data_frames[0] == '/blah/xtal1_1_????.cbf'
    assert handle.name_template_of_data_frames[1] == 'CBF'

    print('OK')

if __name__ == '__main__':

  test = Test()
  test.run()
