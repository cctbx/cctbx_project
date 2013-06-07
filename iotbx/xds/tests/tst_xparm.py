
class Test(object):

  def __init__(self):
    pass

  def run(self):

    from iotbx.xds import xparm
    import os
    import libtbx.load_env
    try:
        iotbx_dir = libtbx.env.dist_path('iotbx')
    except KeyError, e:
        print 'FAIL: dials_regression not configured'
        return

    filename = os.path.join(iotbx_dir, 'xds', 'tests', 'XPARM.XDS')
    handle = xparm.reader()
    handle.read_file(filename)

    # Scan and goniometer stuff
    handle.starting_frame
    handle.starting_angle
    handle.oscillation_range
    handle.rotation_axis

    # Beam stuff
    handle.wavelength
    handle.beam_vector

    # Detector stuff
    handle.detector_size
    handle.pixel_size
    handle.detector_distance
    handle.detector_origin
    handle.detector_x_axis
    handle.detector_y_axis
    handle.detector_normal

    # Crystal stuff
    handle.space_group
    handle.unit_cell
    handle.unit_cell_a_axis
    handle.unit_cell_b_axis
    handle.unit_cell_c_axis

    print 'OK'

if __name__ == '__main__':

  test = Test()
  test.run()
