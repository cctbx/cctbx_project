from __future__ import absolute_import, division, print_function

def test_dps_single_panel_labelit_input_optimal_S0(process_dictionary,data,phil_set):
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
    DPS.set_beam_vector(beam = -beam_vector)
    DPS.set_rotation_axis(axis = rot_axis)
    DPS.set_detector_position(origin = detector_origin, d1 = detector_d1, d2 = detector_d2)

    DPS.index(raw_spot_input = data)
    L = DPS.get_basis_general() # can skip this first time around
    new_S0_vector = DPS.optimize_S0_local_scope()

    DPS2= DPS_primitive_lattice(max_cell = float(process_dictionary['ref_maxcel']),
        recommended_grid_sampling_rad = process_dictionary['recommended_grid_sampling'],
        horizon_phil = phil_set)
    DPS2.set_beam_vector(beam = -new_S0_vector)
    DPS2.set_rotation_axis(axis = rot_axis)
    DPS2.set_detector_position(origin = detector_origin, d1 = detector_d1, d2 = detector_d2)
    DPS2.index(raw_spot_input = data)
    L = DPS2.get_basis_general()
    #new_S0_vector = DPS2.optimize_S0_local_scope() #can skip this second time around

def test_dps_single_panel_labelit_input_optimal_origin(process_dictionary,data,phil_set):
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

    from dxtbx.model import DetectorFactory
    detector = DetectorFactory.make_detector(
      stype = "indexing",
      fast_axis = detector_d1,
      slow_axis = detector_d2,
      origin = detector_origin,
      pixel_size = (float(process_dictionary['pixel_size']),
                    float(process_dictionary['pixel_size'])),
      image_size = (int(process_dictionary['size1']),
                    int(process_dictionary['size2']))
      )

    beam_vector = sample_to_beamspot.normalize() * (
                  1./float(process_dictionary['wavelength']))
    rot_axis = col(process_dictionary["endstation"].rot_axi) - col((0.0,0.0,0.0)) # coerce to float type

    DPS = DPS_primitive_lattice(max_cell = float(process_dictionary['ref_maxcel']),
          recommended_grid_sampling_rad = process_dictionary['recommended_grid_sampling'],
          horizon_phil = phil_set)
    DPS.set_beam_vector(beam = -beam_vector)
    DPS.set_rotation_axis(axis = rot_axis)
    DPS.set_detector(detector)

    DPS.index(raw_spot_input = data)
    #L = DPS.get_basis_general() # can skip this first time around
    new_detector = DPS.optimize_origin_offset_local_scope()

    DPS2= DPS_primitive_lattice(max_cell = float(process_dictionary['ref_maxcel']),
        recommended_grid_sampling_rad = process_dictionary['recommended_grid_sampling'],
        horizon_phil = phil_set)
    DPS2.set_beam_vector(beam = -beam_vector)
    DPS2.set_rotation_axis(axis = rot_axis)
    DPS2.set_detector(new_detector)
    DPS2.index(raw_spot_input = data)
    L = DPS2.get_basis_general()

    #new_conforming_data = DPS2.get_outlier_rejected_subset()

def test_out(process_dictionary,data,phil_set):
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
    DPS.set_beam_vector(beam = -beam_vector)
    DPS.set_rotation_axis(axis = rot_axis)
    DPS.set_detector_position(origin = detector_origin, d1 = detector_d1, d2 = detector_d2)

    DPS.index(raw_spot_input = data)
    L = DPS.get_basis_general()


    new_detector = DPS.optimize_origin_offset_local_scope()

    DPS2= DPS_primitive_lattice(max_cell = float(process_dictionary['ref_maxcel']),
        recommended_grid_sampling_rad = process_dictionary['recommended_grid_sampling'],
        horizon_phil = phil_set)
    DPS2.set_beam_vector(beam = -beam_vector)
    DPS2.set_rotation_axis(axis = rot_axis)
    DPS2.set_detector(new_detector)
    DPS2.index(raw_spot_input = data)
    L = DPS2.get_basis_general()

    from rstbx.indexing_api.outlier_procedure import main_go
    main_go(index_engine=DPS2, phil_set=phil_set)

    print("Finishing")
    exit()

"""
Still to do:
1) Implement origin refinement instead of S0 (DONE)
2) Implement target cell (DONE)
3) Implement outlier rejection (in process)
4) Implement full-parameter refinement
5) Figure out how to evaluate the scoring function analytically & use 2nd-derivative LBFGS
6) Implement quick-refinement of direction vectors as in labelit

Nov. 4:
1) Refine the parameters--triclinic
2) Outlier rejection-essentially done; some refactoring needed.
2b) Rationalize iotbx converter.  Where is the .constrain() applied? Role for dials.crystal models vs. cctbx.crystal symmetry?
3) Refine again & return list of all settings with LABELIT-style output
4) Input can be either cctbx.spotfinder or dials.spotfinder
4b) Document the phil parameters that impinge on my API
5) json output and web-service
6) profile the code as to what is rate limiting
7) multipanel, refined as single block
8) support for Aaron's detector Format
9) multipanel with other parameterizations like quadrant & sensor
"""
