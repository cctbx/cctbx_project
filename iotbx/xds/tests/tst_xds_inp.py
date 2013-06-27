from __future__ import division

class Test(object):

  def __init__(self):
    pass

  def run(self):

    from iotbx.xds import xds_inp
    import os
    import libtbx.load_env
    try:
        iotbx_dir = libtbx.env.dist_path('iotbx')
    except KeyError, e:
        print 'FAIL: dials_regression not configured'
        return

    filename = os.path.join(iotbx_dir, 'xds', 'tests', 'XDS.INP')
    handle = xds_inp.reader()
    handle.read_file(filename)

    handle.detector
    handle.minimum_valid_pixel_value
    handle.overload
    handle.corrections
    handle.direction_of_detector_x_axis
    handle.direction_of_detector_y_axis
    handle.trusted_region
    handle.sensor_thickness
    handle.untrusted_rectangle
    handle.maximum_number_of_processor
    handle.nx
    handle.ny
    handle.px
    handle.py
    handle.orgx
    handle.orgy
    handle.rotation_axis
    handle.detector_distance
    handle.xray_wavelength
    handle.incident_beam_direction
    handle.fraction_of_polarization
    handle.polarization_plane_normal
    handle.friedels_law
    handle.name_template_of_data_frames
    handle.starting_angle
    handle.starting_frane
    handle.include_resolution_range
    handle.unit_cell_constants
    handle.space_group_number
    handle.max_fac_rmeas
    handle.data_range

    print 'OK'

if __name__ == '__main__':

  test = Test()
  test.run()
