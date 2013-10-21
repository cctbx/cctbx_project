from __future__ import division

def test_dps_single_panel_labelit_input(process_dictionary,data,phil_set):
    from rstbx.indexing_api.lattice import DPS_primitive_lattice
    from scitbx.matrix import col

    sample_to_beamspot = col((0.,0.,float(process_dictionary['distance'])))
    detector_d1 = col((1., 0., 0.))
    detector_d2 = col((0., 1., 0.))
    detector_origin = sample_to_beamspot - \
      detector_d1 * float(process_dictionary['xbeam']) - \
      detector_d2 * float(process_dictionary['ybeam'])
    assert detector_d1.length() == 1.0
    assert detector_d2.length() == 1.0
    assert detector_d1.dot(detector_d2) == 0.0

    beam_vector = sample_to_beamspot.normalize() * (
                  1./float(process_dictionary['wavelength']))
    rot_axis = col(process_dictionary["endstation"].rot_axi) - col((0.0,0.0,0.0)) # coerce to float type

    DPS = DPS_primitive_lattice(max_cell = float(process_dictionary['ref_maxcel']),
          recommended_grid_sampling_rad = process_dictionary['recommended_grid_sampling'],
          horizon_phil = phil_set)
    DPS.set_beam_vector(beam = beam_vector)
    DPS.set_rotation_axis(axis = rot_axis)
    DPS.set_detector_position(origin = detector_origin, d1 = detector_d1, d2 = detector_d2)

    L = DPS.index(raw_spot_input = data)

    new_beam = DPS.new_beam

    DPS2= DPS_primitive_lattice(max_cell = float(process_dictionary['ref_maxcel']),
        recommended_grid_sampling_rad = process_dictionary['recommended_grid_sampling'],
        horizon_phil = phil_set)
    DPS2.set_beam_vector(beam = new_beam)
    DPS2.set_rotation_axis(axis = rot_axis)
    DPS2.set_detector_position(origin = detector_origin, d1 = detector_d1, d2 = detector_d2)
    L = DPS2.index(raw_spot_input = data)

