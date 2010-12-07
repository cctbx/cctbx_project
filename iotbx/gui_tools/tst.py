
from iotbx.gui_tools import reflections
import libtbx.load_env
from libtbx.test_utils import approx_equal, Exception_expected
from iotbx import file_reader
from libtbx.utils import Sorry
import os

def exercise_reflections () :
  phil_names = ["refinement.input.xray_data.file_name",
                "refinement.input.xray_data.r_free_flags.file_name",
                "refinement.input.neutron_data.file_name",
                "refinement.input.neutron_data.r_free_flags.file_name",
                "refinement.input.experimental_phases.file_name"]
  hkl_handler = reflections.reflections_handler(allowed_param_names=phil_names)
  file_name = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/reflection_files/xn_data_ricardo_leal.mtz",
    test=os.path.isfile)
  assert (hkl_handler.set_param_file(file_name=file_name,
          file_param_name="refinement.input.xray_data.file_name") == True)
  for i, phil_name in enumerate(phil_names[1:4]) :
    assert (hkl_handler.set_param_file(file_name=file_name,
            file_param_name=phil_name) == False)
  assert (hkl_handler.get_rfree_labels(
    file_param_name="refinement.input.xray_data.r_free_flags.file_name") ==
    hkl_handler.get_rfree_labels(file_name=file_name) ==
    ['R-free-flags', 'R-free-flags-neutron'])
  hkl_handler.check_symmetry(file_name=file_name)
  assert (hkl_handler.get_data_labels(
          file_param_name="refinement.input.xray_data.file_name") ==
          hkl_handler.get_data_labels(file_name=file_name) ==
          ["F-obs,SIGF-obs", "F-obs-neutron(+),SIGF-obs-neutron(+)," +
                             "F-obs-neutron(-),SIGF-obs-neutron(-)"])
  assert (hkl_handler.get_anomalous_data_labels(
          file_param_name="refinement.input.xray_data.file_name") ==
          ["F-obs-neutron(+),SIGF-obs-neutron(+)," +
           "F-obs-neutron(-),SIGF-obs-neutron(-)"])
  assert (hkl_handler.get_rfree_flag_value(array_name="F-obs,SIGF-obs",
    file_param_name="refinement.input.xray_data.file_name") is None)
  assert (hkl_handler.get_rfree_flag_value(array_name='R-free-flags',
    file_param_name="refinement.input.xray_data.r_free_flags.file_name") == 1)
  (d_max, d_min) = hkl_handler.d_max_min()
  assert approx_equal(d_max, 28.61, eps=0.01)
  assert approx_equal(d_min, 1.93, eps=0.01)
  assert hkl_handler.space_group_as_str() == "P 61 2 2"
  assert (hkl_handler.unit_cell_as_str() ==
          "33.034 33.034 78.440 90.000 90.000 120.000")
  assert (hkl_handler.unit_cell_as_str(separator=",") ==
          "33.034,33.034,78.440,90.000,90.000,120.000")
  cns_file = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/reflection_files/enk.hkl",
    test=os.path.isfile)
  hkl_handler.save_file(file_name=cns_file)
  try :
    hkl_handler.check_symmetry(file_name=cns_file)
  except Sorry :
    pass
  else :
    raise Exception_expected

  resolve_file = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/wizards/resolve_1_offset.mtz",
    test=os.path.isfile)
  hkl_handler = reflections.reflections_handler(
    allowed_param_names=["map_coeffs"])
  hkl_handler.set_param_file(
    file_name=resolve_file,
    file_param_name="map_coeffs")
  assert hkl_handler.has_data(file_name=resolve_file)
  l1 = hkl_handler.get_map_coeff_labels(file_name=resolve_file)
  assert (l1 == ['FP,PHIM,FOMM'])
  l2 = hkl_handler.get_phase_deg_labels(file_name=resolve_file)
  assert (l2 == ['HLAM,HLBM,HLCM,HLDM', 'PHIM', 'FOMM'])
  l3 = hkl_handler.get_experimental_phase_labels(file_name=resolve_file)
  #print l3
  l4 = hkl_handler.get_data_labels_for_wizard(file_name=resolve_file)
  assert (l4 == ['FP SIGFP'])
  l5 = hkl_handler.get_map_coeff_labels_for_build(file_name=resolve_file)
  assert (l5 == ['FP,PHIM,FOMM'])
  hkl_in = hkl_handler.get_file(file_name=resolve_file)
  assert (reflections.get_mtz_label_prefix(hkl_in) == "/NoName/NoName")
  map_coeffs = reflections.map_coeffs_from_mtz_file(resolve_file,
    f_label="FP,SIGFP")
  assert map_coeffs.is_complex_array()
  try :
    map_coeffs = reflections.map_coeffs_from_mtz_file(resolve_file)
  except Sorry :
    pass
  else :
    raise Exception_expected
  assert (hkl_handler.get_map_coeff_labels_for_fft(file_name=resolve_file) ==
    ['FP,SIGFP PHIM FOMM'])

  # miscellaneous utilities
  file_name = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/wizards/map_coeffs.mtz",
    test=os.path.isfile)
  hkl_in = file_reader.any_file(file_name)
  hkl_server = hkl_in.file_server
  assert approx_equal(reflections.get_high_resolution(hkl_server), 2.0,
    eps=0.0001)
  descriptions = []
  for miller_array in hkl_server.miller_arrays :
    (sg, uc) = reflections.get_miller_array_symmetry(miller_array)
    assert (uc == "60.832 38.293 42.211  90 90 90")
    assert str(sg) == "C 2 2 21"
    descriptions.append(reflections.get_array_description(miller_array))
  assert descriptions == [
    'Amplitude', 'Phases', 'Weights', 'HL coeffs', 'R-free flag']
  handler = reflections.reflections_handler()
  handler.save_file(input_file=hkl_in)
  assert (handler.get_resolution_range(file_name=file_name)=="(19.146 - 2.000)")
  assert (handler.get_resolution_limits(file_name=file_name) ==
          ('(19.146)', '(2.000)'))
  file_name = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/wizards/partial_refine_001_map_coeffs.mtz",
    test=os.path.isfile)
  hkl_server = file_reader.any_file(file_name).file_server
  map_labels = reflections.get_map_coeff_labels(hkl_server)
  assert (map_labels == ['2FOFCWT,PH2FOFCWT', 'FOFCWT,PHFOFCWT',
    '2FOFCWT_no_fill,PH2FOFCWT_no_fill', 'FOFCWT_no_fill,PHFOFCWT_no_fill'])
  map_labels = reflections.get_map_coeffs_for_build(hkl_server)
  assert map_labels == ['2FOFCWT,PH2FOFCWT','2FOFCWT_no_fill,PH2FOFCWT_no_fill']
  map_coeffs = reflections.extract_phenix_refine_map_coeffs(file_name)
  assert (len(map_coeffs) == 4)
  hkl_file = file_reader.any_file(file_name)
  assert reflections.get_mtz_label_prefix(hkl_file) == "/crystal/dataset"
  (fp, fpp) = reflections.get_fp_fpp_from_sasaki("Se", 0.979)
  assert fp is not None and fpp is not None

if (__name__ == "__main__") :
  hkl_dir = libtbx.env.find_in_repositories(
    relative_path="phenix_regression",
    test=os.path.isdir)
  if (hkl_dir is None) :
    print "phenix_regression/reflection_files not found, skipping tests."
  else :
    exercise_reflections()
    print "OK"
