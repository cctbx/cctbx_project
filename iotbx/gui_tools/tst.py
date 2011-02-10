
from iotbx import gui_tools
from iotbx.gui_tools import reflections, models
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
  sca_file = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/reflection_files/merge.sca",
    test=os.path.isfile)
  hkl_handler.save_file(file_name=sca_file)
  assert (hkl_handler.get_intensity_labels(file_name=sca_file) ==
          ['i_obs,sigma'])

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

def exercise_model () :
  model_handler = models.model_handler(
    allowed_param_names=["refinement.input.pdb.file_name"],
    allowed_multiple_params=["refinement.input.pdb.file_name"],
    cif_param_names=["refinement.input.monomers.file_name"],
    multiple_cif_params=["refinement.input.monomers.file_name"],
    tmp_dir=os.getcwd())
  pdb_file = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/pdb/ur0013.pdb",
    test=os.path.isfile)
  cif_file = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/cif_files/elbow.ur0013_ent.all.001.cif",
    test=os.path.isfile)
  model_handler.set_param_file(
    file_name=pdb_file,
    file_param_name="refinement.input.pdb.file_name")
  model_handler.set_param_cif_file(
    file_name=cif_file,
    file_param_name="refinement.input.monomers.file_name")
  assert (len(model_handler.get_cif_objects()) == 1)
  assert (model_handler.get_current_cif_file_names() == [cif_file])
  pdb_file2 = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/pdb/1ywf.pdb",
    test=os.path.isfile)
  model_handler.set_param_file(
    file_name=pdb_file2,
    file_param_name="refinement.input.pdb.file_name")
  assert (model_handler.get_current_file_names() == [pdb_file2, pdb_file])
  assert (model_handler.get_file_params(pdb_file) ==
          ["refinement.input.pdb.file_name"])
  assert (model_handler.get_param_files("refinement.input.pdb.file_name") ==
          [pdb_file, pdb_file2])
  pdb_hierarchy = model_handler.get_pdb_hierarchy(pdb_file2)
  xrs = model_handler.get_xray_structure(pdb_file2)
  assert (pdb_hierarchy.atoms().size() == xrs.scatterers().size() == 2127)
  model_handler.remove_file(pdb_file)
  assert (model_handler.get_param_files("refinement.input.pdb.file_name") ==
          [pdb_file2])
  atomic_bonds = model_handler.get_connectivity(pdb_file2)
  assert (atomic_bonds.size() == 2127)
  symm = model_handler.get_pdb_file_symmetry(pdb_file2)
  assert (str(symm.space_group_info()) == "I 41")
  assert (reflections.unit_cell_as_str(symm.unit_cell()) ==
          "113.068 113.068 53.292 90.000 90.000 90.000")
  f = model_handler.create_copy_with_fake_symmetry(pdb_file2,
    tmp_dir=os.getcwd())
  pdb_in = file_reader.any_file(f, force_type="pdb").file_object
  symm = pdb_in.crystal_symmetry()
  assert (str(symm.space_group_info()) == "P 1")
  assert (reflections.unit_cell_as_str(symm.unit_cell()) ==
          "59.227 55.922 60.264 90.000 90.000 90.000")

def exercise_symmetry () :
  from cctbx import sgtbx, uctbx
  m = gui_tools.symmetry_manager(prefer_pdb_space_group=True)
  (uc_mismatch, sg_mismatch) = m.add_reflections_file(
    file_name="data.mtz",
    space_group=sgtbx.space_group_info("P222"),
    unit_cell=uctbx.unit_cell("50 60 70 90 90 90"))
  assert (m.get_current_as_strings() == ('P 2 2 2', '50 60 70 90 90 90'))
  (uc_mismatch, sg_mismatch) = m.add_pdb_file(
    file_name="model.pdb",
    space_group=sgtbx.space_group_info("P212121"),
    unit_cell=uctbx.unit_cell("50 60 70 90 90 90"))
  assert (not (uc_mismatch or sg_mismatch))
  (uc_mismatch, sg_mismatch) = m.add_pdb_file(
    file_name="reference_model.pdb",
    space_group=sgtbx.space_group_info("P63"),
    unit_cell=uctbx.unit_cell("40 40 75 90 90 120"))
  assert ((uc_mismatch, sg_mismatch) == (True, True))
  assert (m.get_current_as_strings() == ('P 21 21 21', '50 60 70 90 90 90'))
  (uc_mismatch, sg_mismatch) = m.add_reflections_file(
    file_name="data_neutron.mtz",
    space_group=sgtbx.space_group_info("P222"),
    unit_cell=uctbx.unit_cell("50.1 60 70.1 90 90 90"))
  assert (not (uc_mismatch or sg_mismatch))
  (uc_mismatch, sg_mismatch) = m.add_reflections_file(
    file_name="data_rfree.hkl",
    space_group=None,
    unit_cell=None)
  assert (not (uc_mismatch or sg_mismatch))
  assert (m.get_current_as_strings() == ('P 21 21 21', '50 60 70 90 90 90'))
  assert (m.check_cell_compatibility("phenix.refine"))
  symm_choices = m.get_symmetry_choices()
  assert (symm_choices.space_group_files == [('model.pdb', 'P 21 21 21'),
    ('reference_model.pdb', 'P 63'), ('data.mtz', 'P 2 2 2'),
    ('data_neutron.mtz', 'P 2 2 2')])
  assert (symm_choices.unit_cell_files == [
    ('model.pdb', '(50, 60, 70, 90, 90, 90)'),
    ('reference_model.pdb', '(40, 40, 75, 90, 90, 120)'),
    ('data.mtz', '(50, 60, 70, 90, 90, 90)'),
    ('data_neutron.mtz', '(50.1, 60, 70.1, 90, 90, 90)')])
  m.set_current_as_strings("P63", "50 60 70 90 90 90")
  try :
    m.check_cell_compatibility(
      program_name="phenix.refine",
      raise_error_if_incomplete=True)
  except Sorry :
    pass
  else :
    raise Exception_expected

if (__name__ == "__main__") :
  exercise_symmetry() # this doesn't need files
  hkl_dir = libtbx.env.find_in_repositories(
    relative_path="phenix_regression",
    test=os.path.isdir)
  if (hkl_dir is None) :
    print "phenix_regression/reflection_files not found, skipping tests."
  else :
    exercise_model()
    exercise_reflections()
    print "OK"
