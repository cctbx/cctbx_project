from __future__ import absolute_import, division, print_function

from iotbx.gui_tools import reflections, models
from iotbx import file_reader
import iotbx.pdb
import libtbx.load_env
from libtbx.test_utils import approx_equal, Exception_expected
from libtbx.utils import Sorry
import os

def find_file(file_name):
  full_path = libtbx.env.find_in_repositories(
    relative_path=file_name,
    test=os.path.isfile)
  return full_path

regression_dir = libtbx.env.find_in_repositories(
  relative_path="phenix_regression",
  test=os.path.isdir)

phil_names = ["refinement.input.xray_data.file_name",
              "refinement.input.xray_data.r_free_flags.file_name",
              "refinement.input.neutron_data.file_name",
              "refinement.input.neutron_data.r_free_flags.file_name",
              "refinement.input.experimental_phases.file_name"]

def exercise_reflections():
  hkl_handler = reflections.reflections_handler(allowed_param_names=phil_names)
  from cctbx import miller
  from cctbx import crystal
  from cctbx.array_family import flex
  symm = crystal.symmetry(
    unit_cell=(30,30,40,90,90,120),
    space_group_symbol="P 61 2 2")
  miller_set = miller.build_set(
    crystal_symmetry=symm,
    anomalous_flag=True,
    d_min=1.5)
  n_refl = miller_set.indices().size()
  data = flex.random_double(n_refl)
  sigmas = flex.random_double(n_refl)
  f_obs = miller_set.array(data=data, sigmas=sigmas)
  f_obs_merged = f_obs.average_bijvoet_mates()
  flags = f_obs_merged.generate_r_free_flags()
  # single dataset
  mtz_dataset = f_obs_merged.as_mtz_dataset(
    column_root_label="F-obs",
    wavelength=1.54)
  mtz_dataset.add_miller_array(flags,
    column_root_label="R-free-flags")
  file_name = "tst_iotbs_gui_tools.mtz"
  mtz_dataset.mtz_object().write(file_name)
  assert (hkl_handler.set_param_file(file_name=file_name,
          file_param_name="refinement.input.xray_data.file_name") == True)
  assert (hkl_handler.set_param_file(file_name=file_name,
          file_param_name="refinement.input.xray_data.r_free_flags.file_name")
          == False)
  assert (hkl_handler.get_rfree_labels(
    file_param_name="refinement.input.xray_data.r_free_flags.file_name") ==
    hkl_handler.get_rfree_labels(file_name=file_name) == ['R-free-flags'])
  assert (hkl_handler.get_rfree_labels(file_name=file_name, neutron=False) ==
          hkl_handler.get_rfree_labels(file_name=file_name, neutron=True) ==
          ['R-free-flags'])
  assert approx_equal(1.54, hkl_handler.get_wavelength(file_name=file_name,
                              labels="F-obs,SIGF-obs"))
  # join X/N datasets
  hkl_handler = reflections.reflections_handler(allowed_param_names=phil_names)
  data_neutron = flex.random_double(n_refl)
  sigmas_neutron = flex.random_double(n_refl)
  f_obs_neutron = miller_set.array(data=data_neutron, sigmas=sigmas_neutron)
  mtz_dataset = f_obs_merged.as_mtz_dataset(
    column_root_label="F-obs-xray")
  mtz_dataset.add_miller_array(f_obs_neutron,
    column_root_label="F-obs-neutron")
  mtz_dataset.add_miller_array(flags,
    column_root_label="R-free-flags-xray")
  mtz_dataset.add_miller_array(flags.deep_copy(),
    column_root_label="R-free-flags-neutron")
  file_name = "tst_iotbs_gui_tools.mtz"
  mtz_dataset.mtz_object().write(file_name)
  assert (hkl_handler.set_param_file(file_name=file_name,
          file_param_name="refinement.input.xray_data.file_name") == True)
  for i, phil_name in enumerate(phil_names[1:4]):
    assert (hkl_handler.set_param_file(file_name=file_name,
            file_param_name=phil_name) == False)
  assert (hkl_handler.get_rfree_labels(
    file_param_name="refinement.input.xray_data.r_free_flags.file_name") ==
    hkl_handler.get_rfree_labels(file_name=file_name) ==
    ['R-free-flags-xray', 'R-free-flags-neutron'])
  assert (hkl_handler.get_rfree_labels(file_name=file_name, neutron=False) ==
          ['R-free-flags-xray'])
  assert (hkl_handler.get_rfree_labels(file_name=file_name, neutron=True) ==
          ['R-free-flags-neutron'])
  hkl_handler.check_symmetry(file_name=file_name)
  assert (hkl_handler.get_data_labels(
          file_param_name="refinement.input.xray_data.file_name") ==
          hkl_handler.get_data_labels(file_name=file_name) ==
          hkl_handler.get_amplitude_labels(file_name=file_name) ==
          ["F-obs-xray,SIGF-obs-xray", "F-obs-neutron(+),SIGF-obs-neutron(+),"+
                             "F-obs-neutron(-),SIGF-obs-neutron(-)"])
  assert (hkl_handler.get_anomalous_data_labels(
          file_param_name="refinement.input.xray_data.file_name") ==
          ["F-obs-neutron(+),SIGF-obs-neutron(+)," +
           "F-obs-neutron(-),SIGF-obs-neutron(-)"])
  assert (hkl_handler.has_anomalous_data(
          file_param_name="refinement.input.xray_data.file_name"))
  assert (hkl_handler.get_rfree_flag_value(
    array_name="F-obs-xray,SIGF-obs-xray",
    file_param_name="refinement.input.xray_data.file_name") is None)
  assert (hkl_handler.get_rfree_flag_value(array_name='R-free-flags-xray',
    file_param_name="refinement.input.xray_data.r_free_flags.file_name") == 1)
  (d_max, d_min) = hkl_handler.d_max_min()
  assert approx_equal(d_max, 25.98, eps=0.01)
  assert approx_equal(d_min, 1.5, eps=0.01)
  assert hkl_handler.space_group_as_str() == "P 61 2 2"
  assert (hkl_handler.unit_cell_as_str() ==
          "30 30 40 90 90 120")
  assert (hkl_handler.unit_cell_as_str(separator=",") ==
          "30,30,40,90,90,120")

  n_refl_merged = len(f_obs_merged.indices())
  phi_array = f_obs_merged.random_phases_compatible_with_phase_restrictions(
    deg=True)
  fom_array = phi_array.array(data=flex.random_double(n_refl_merged))
  hl_data = flex.hendrickson_lattman(n_refl_merged, (0,0,0,0))
  hl_coeffs = phi_array.array(data=hl_data)
  assert (hl_coeffs.is_hendrickson_lattman_array())
  from iotbx.mtz import label_decorator
  import iotbx.mtz
  class resolve_label_decorator(label_decorator):
    def phases(self, *args, **kwds):
      return label_decorator.phases(self, *args, **kwds) + "M"
    def hendrickson_lattman(self, *args, **kwds):
      return label_decorator.hendrickson_lattman(self, *args, **kwds) + "M"
  mtz_dataset = f_obs_merged.as_mtz_dataset(
    column_root_label="FP")
  mtz_dataset.add_miller_array(phi_array,
    column_root_label="PHIM",
    label_decorator=resolve_label_decorator(),
    column_types="P")
  mtz_dataset.add_miller_array(fom_array,
    column_root_label="FOMM",
    column_types="W")
  mtz_dataset.add_miller_array(hl_coeffs,
    column_root_label="HL",
    label_decorator=resolve_label_decorator())
  fwt_map = f_obs_merged.customized_copy(sigmas=None).phase_transfer(phi_array)
  mtz_dataset.add_miller_array(fwt_map,
    column_root_label="FWT",
    label_decorator=iotbx.mtz.ccp4_label_decorator())
  mtz_dataset.add_miller_array(flags,
    column_root_label="FreeR_flag")
  resolve_file = "tst_iotbx_gui_tools_resolve.mtz"
  mtz_dataset.mtz_object().write(resolve_file)
  hkl_handler = reflections.reflections_handler(
    allowed_param_names=["map_coeffs"])
  hkl_handler.set_param_file(
    file_name=resolve_file,
    file_param_name="map_coeffs")
  assert hkl_handler.has_data(file_name=resolve_file)
  l1 = hkl_handler.get_map_coeff_labels(file_name=resolve_file)
  assert (l1 == ['FWT,PHWT', 'FP,PHIM,FOMM'])
  l2 = hkl_handler.get_phase_deg_labels(file_name=resolve_file)
  assert (l2 == ['HLAM,HLBM,HLCM,HLDM', 'FWT,PHWT', 'PHIM',]), l2
  l3 = hkl_handler.get_experimental_phase_labels(file_name=resolve_file)
  #print l3
  l4 = hkl_handler.get_data_labels_for_wizard(file_name=resolve_file)
  assert (l4 == ['FP SIGFP'])
  l5 = hkl_handler.get_map_coeff_labels_for_build(file_name=resolve_file)
  assert (l5 == ['FWT,PHWT', 'FP,PHIM,FOMM']), l5
  hkl_in = hkl_handler.get_file(file_name=resolve_file)
  assert (reflections.get_mtz_label_prefix(hkl_in) == "/crystal/dataset")
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
    ['FWT,PHWT', 'FP,SIGFP PHIM FOMM'])

  # miscellaneous utilities
  file_name = resolve_file
  hkl_in = file_reader.any_file(file_name)
  hkl_server = hkl_in.file_server
  assert approx_equal(reflections.get_high_resolution(hkl_server), 1.5,
    eps=0.001)
  descriptions = []
  for miller_array in hkl_server.miller_arrays :
    (sg, uc) = reflections.get_miller_array_symmetry(miller_array)
    assert (uc == "30 30 40 90 90 120")
    assert (str(sg) == "P 61 2 2")
    descriptions.append(reflections.get_array_description(miller_array))
  assert (descriptions == [
    'Amplitude', 'Phases', 'Weights', 'HL coeffs', 'Map coeffs',
    'R-free flag']), descriptions
  handler = reflections.reflections_handler()
  handler.save_file(input_file=hkl_in)
  assert (not handler.has_anomalous_data())
  assert (handler.get_resolution_range(file_name=file_name)=="(25.981 - 1.500)"
          or handler.get_resolution_range(file_name=file_name)=="(25.981 - 1.501)")
  assert (handler.get_resolution_limits(file_name=file_name) == ('(25.981)', '(1.500)')
          or handler.get_resolution_limits(file_name=file_name) == ('(25.981)', '(1.501)'))
  fmodel = phi_array.array(data=flex.complex_double(n_refl_merged,
    complex(0.5,0.8)))
  m1 = phi_array.array(data=flex.complex_double(n_refl_merged, complex(1,0)))
  m2 = phi_array.array(data=flex.complex_double(n_refl_merged, complex(0.5,0)))
  m3 = phi_array.array(flex.complex_double(n_refl_merged, complex(1,1)))
  dec = label_decorator(phases_prefix="PH")
  mtz_dataset = fmodel.as_mtz_dataset(
    column_root_label="F-model",
    label_decorator=dec)
  mtz_dataset.add_miller_array(m1,
    column_root_label="2FOFCWT",
    label_decorator=dec)
  mtz_dataset.add_miller_array(m2,
    column_root_label="FOFCWT",
    label_decorator=dec)
  mtz_dataset.add_miller_array(m3,
    column_root_label="2FOFCWT_no_fill",
    label_decorator=dec)
  file_name = "tst_iotbx_gui_tools_map_coeffs.mtz"
  mtz_dataset.mtz_object().write(file_name)
  hkl_handler = reflections.reflections_handler(
    allowed_param_names=["fmodel", "map_coeffs"])
  hkl_handler.set_param_file(
    file_name=file_name,
    file_param_name="fmodel")
  assert (hkl_handler.get_fmodel_labels(file_name=file_name) ==
    ['F-model,PHF-model'])
  assert (hkl_handler.get_amplitude_labels(file_name=file_name) == [])
  phi_labels = hkl_handler.get_phase_deg_labels(file_name=file_name)
  assert (len(phi_labels)  == 4)
  phi_cols = hkl_handler.get_phase_column_labels(file_name=file_name)
  assert (phi_cols == ['PHF-model','PH2FOFCWT','PHFOFCWT','PH2FOFCWT_no_fill'])
  assert (len(hkl_handler.get_amplitude_column_labels(file_name=file_name,
              allow_conversion=True)) == 0)
  fc_cols = hkl_handler.get_fmodel_labels(file_name=file_name,
    first_column_only=True)
  assert (fc_cols == ['F-model'])
  hkl_server = file_reader.any_file(file_name).file_server
  map_labels = reflections.get_map_coeff_labels(hkl_server)
  assert (map_labels == ['F-model,PHF-model', '2FOFCWT,PH2FOFCWT', 'FOFCWT,PHFOFCWT',
    '2FOFCWT_no_fill,PH2FOFCWT_no_fill',]), map_labels
  map_labels = reflections.get_map_coeffs_for_build(hkl_server)
  assert map_labels == ['2FOFCWT,PH2FOFCWT', 'F-model,PHF-model', '2FOFCWT_no_fill,PH2FOFCWT_no_fill'], map_labels
  map_coeffs = reflections.extract_phenix_refine_map_coeffs(file_name)
  assert (len(map_coeffs) == 3)
  hkl_file = file_reader.any_file(file_name)
  assert reflections.get_mtz_label_prefix(hkl_file) == "/crystal/dataset"
  # other stuff
  (fp, fpp) = reflections.get_fp_fpp_from_sasaki("Se", 0.979)
  assert fp is not None and fpp is not None

def exercise_other_reflection_formats():
  hkl_handler = reflections.reflections_handler(allowed_param_names=phil_names)
  # test other file formats (requires phenix_regression)
  if (regression_dir is None):
    print("phenix_regression not found, skipping exercise_other_reflection_formats()")
    return
  cns_file = os.path.join(regression_dir, "reflection_files", "enk.hkl")
  hkl_handler.save_file(file_name=cns_file)
  try :
    hkl_handler.check_symmetry(file_name=cns_file)
  except Sorry :
    pass
  else :
    raise Exception_expected
  sca_file = os.path.join(regression_dir, "reflection_files", "merge.sca")
  hkl_handler.save_file(file_name=sca_file)
  assert (hkl_handler.get_intensity_labels(file_name=sca_file) ==
          ['I(+),SIGI(+),I(-),SIGI(-)'])
  assert (hkl_handler.get_amplitude_labels(file_name=sca_file) == [])
  # test handling of reconstructed (anomalous) amplitudes
  dano_file = os.path.join(regression_dir, "reflection_files", "dano.mtz")
  hkl_handler = reflections.reflections_handler(
    allowed_param_names=["labin"])
  hkl_handler.set_param_file(file_name=dano_file,
    file_param_name="labin")
  labels = hkl_handler.get_anomalous_data_labels(file_param_name="labin")
  assert (len(labels) == 3)
  labels = hkl_handler.get_anomalous_data_labels(file_param_name="labin",
    allow_reconstructed_amplitudes=False)
  assert (len(labels) == 0)

def exercise_model():
  # FIXME should be possible to run this independently of phenix_regression
  if (regression_dir is None):
    print("phenix_regression not found, skipping exercise_model()")
    return
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
  try :
    pdb_all = model_handler.get_combined_pdb_input(
      file_param_name="refinement.input.pdb.file_name")
  except Sorry :
    pass
  else :
    raise Exception_expected
  model_handler.set_param_file(
    file_name=pdb_file,
    file_param_name="refinement.input.pdb.file_name")
  assert (model_handler.get_file_type_label(file_name=pdb_file) == "PDB")
  model_handler.set_param_cif_file(
    file_name=cif_file,
    file_param_name="refinement.input.monomers.file_name")
  assert (len(model_handler.get_cif_objects()) == 1)
  assert (model_handler.get_current_cif_file_names() == [cif_file])
  hierarchy, xrs = model_handler.get_combined_pdb_input(
      file_param_name="refinement.input.pdb.file_name")
  assert (not None in [hierarchy, xrs])
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
  assert (pdb_hierarchy.atoms_size() == xrs.scatterers().size() == 2127)
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
  pdb_in = iotbx.pdb.input(f)
  symm = pdb_in.crystal_symmetry()
  assert (str(symm.space_group_info()) == "P 1")
  assert (reflections.unit_cell_as_str(symm.unit_cell()) ==
          "59.227 55.922 60.264 90 90 90")

if (__name__ == "__main__"):
  exercise_reflections()
  exercise_other_reflection_formats()
  exercise_model()
  print("OK")
