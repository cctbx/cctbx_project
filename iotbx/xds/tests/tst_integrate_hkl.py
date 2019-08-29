from __future__ import absolute_import, division, print_function

class Test(object):

  def __init__(self):
    pass

  def run(self):

    from iotbx.xds import integrate_hkl
    import os
    import libtbx.load_env
    iotbx_dir = libtbx.env.dist_path('iotbx')
    filename = os.path.join(iotbx_dir, 'xds', 'tests', 'INTEGRATE.HKL')
    handle = integrate_hkl.reader()
    handle.read_file(filename)

    handle.space_group
    handle.unit_cell
    handle.detector_size
    handle.pixel_size
    handle.starting_frame
    handle.starting_angle
    handle.oscillation_range
    handle.rotation_axis
    handle.wavelength
    handle.beam_vector
    handle.detector_x_axis
    handle.detector_y_axis
    handle.detector_origin
    handle.detector_distance
    handle.unit_cell_a_axis
    handle.unit_cell_b_axis
    handle.unit_cell_c_axis
    handle.sigma_divergence
    handle.sigma_mosaicity
    handle.template
    handle.detector_type
    handle.minpk
    handle.cut
    handle.variance_model

    from iotbx.reflection_file_reader import any_reflection_file
    file_in = any_reflection_file(filename)
    assert file_in.file_type() == 'xds_integrate_hkl'
    miller_arrays = file_in.as_miller_arrays()
    assert len(miller_arrays) == 10
    assert miller_arrays[0].space_group().type().number() == handle.space_group
    content = file_in.file_content()
    assert content.beam_vector == (-0.001316, 0.001644, 1.020927)
    assert content.wavelength == 0.9795

    print('OK')

if __name__ == '__main__':
  test = Test()
  test.run()
