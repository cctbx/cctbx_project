
from iotbx.gui_tools import reflections, models, maps
from iotbx import file_reader
import libtbx.load_env
from libtbx.test_utils import approx_equal, Exception_expected
from libtbx.utils import Sorry
from libtbx import easy_run
import os

def find_file (file_name) :
  full_path = libtbx.env.find_in_repositories(
    relative_path=file_name,
    test=os.path.isfile)
  return full_path

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
          "33.0343 33.0343 78.4404 90 90 120")
  assert (hkl_handler.unit_cell_as_str(separator=",") ==
          "33.0343,33.0343,78.4404,90,90,120")
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
          ['I(+),SIGI(+),I(-),SIGI(-)'])

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
    assert (uc == "60.832 38.293 42.211 90 90 90")
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
    allowed_param_names=["refinement.input.pdb.file_name",
      "refinement.reference_model.file"],
    allowed_multiple_params=["refinement.input.pdb.file_name"],
    cif_param_names=["refinement.input.monomers.file_name"],
    multiple_cif_params=["refinement.input.monomers.file_name"],
    tmp_dir=os.getcwd())
  model_handler.set_viewable_params(["refinement.input.pdb.file_name"])
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
  pdb_file3 = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/pdb/1ywf_h.pdb",
    test=os.path.isfile)
  model_handler.set_param_file(
    file_name=pdb_file3,
    file_param_name="refinement.reference_model.file")
  assert (model_handler.get_files_for_viewing() == [pdb_file2])
  atomic_bonds = model_handler.get_connectivity(pdb_file2)
  assert (atomic_bonds.size() == 2127)
  symm = model_handler.get_pdb_file_symmetry(pdb_file2)
  assert (str(symm.space_group_info()) == "I 41")
  assert (reflections.unit_cell_as_str(symm.unit_cell()) ==
          "113.068 113.068 53.292 90 90 90")
  f = model_handler.create_copy_with_fake_symmetry(pdb_file2,
    tmp_dir=os.getcwd())
  pdb_in = file_reader.any_file(f, force_type="pdb").file_object
  symm = pdb_in.crystal_symmetry()
  assert (str(symm.space_group_info()) == "P 1")
  assert (reflections.unit_cell_as_str(symm.unit_cell()) ==
          "59.227 55.922 60.264 90 90 90")

def exercise_maps () :
  pdb_file = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/wizards/partial_refine_001.pdb",
    test=os.path.isfile)
  refine_file = pdb_file[:-4] + "_map_coeffs.mtz"
  easy_run.call("cp %s ." % refine_file)
  easy_run.call("rm -f *.ccp4")
  server = maps.server(file_name="partial_refine_001_map_coeffs.mtz")
  files = server.convert_phenix_maps(file_base="refine")
  assert (files == ['refine_mFo-DFc.ccp4', 'refine_mFo-DFc_2.ccp4',
                    'refine_2mFo-DFc.ccp4', 'refine_2mFo-DFc_no_fill.ccp4'])
  for fn in files :
    assert os.path.isfile(fn)
  resolve_file = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/wizards/resolve_1_offset.mtz",
    test=os.path.isfile)
  easy_run.call("cp %s ." % resolve_file)
  server = maps.server(file_name="resolve_1_offset.mtz")
  map_file = server.convert_resolve_map(pdb_file=None,
    force=False)
  assert (map_file == 'resolve_1_offset.ccp4')
  assert os.path.exists(map_file)

if (__name__ == "__main__") :
  hkl_dir = libtbx.env.find_in_repositories(
    relative_path="phenix_regression",
    test=os.path.isdir)
  if (hkl_dir is None) :
    print "phenix_regression/reflection_files not found, skipping tests."
  else :
    exercise_maps()
    exercise_model()
    exercise_reflections()
    print "OK"
