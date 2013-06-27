from __future__ import division

class Test(object):

  def __init__(self):
    pass

  def run(self):

    from iotbx.xds import integrate_hkl
    import os
    import libtbx.load_env
    try:
        iotbx_dir = libtbx.env.dist_path('iotbx')
    except KeyError, e:
        print 'FAIL: dials_regression not configured'
        return

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

    print 'OK'

if __name__ == '__main__':

  test = Test()
  test.run()
