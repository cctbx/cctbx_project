
from __future__ import absolute_import, division, print_function
from libtbx.test_utils import Exception_expected
from libtbx.utils import null_out, Sorry
import iotbx.pdb
from six.moves import cStringIO as StringIO
import os.path

def exercise():
  from mmtbx.regression.make_fake_anomalous_data import generate_cd_cl_inputs
  from mmtbx.command_line import mtz2map
  import mmtbx.utils
  from iotbx import file_reader
  from scitbx.array_family import flex
  mtz_file, pdb_file = generate_cd_cl_inputs(file_base = "tst_mmtbx_mtz2map")
  pdb_in = iotbx.pdb.input(pdb_file)
  xrs = pdb_in.xray_structure_simple()
  mtz_in = file_reader.any_file(mtz_file)
  f_obs = mtz_in.file_server.miller_arrays[0]
  f_obs_mean = f_obs.average_bijvoet_mates()
  flags = mtz_in.file_server.miller_arrays[1]
  flags = flags.customized_copy(data=flags.data()==1)
  fmodel = mmtbx.utils.fmodel_simple(
    f_obs=f_obs,
    r_free_flags=flags,
    xray_structures=[xrs],
    scattering_table="n_gaussian",
    skip_twin_detection=True)
  assert f_obs.anomalous_flag()
  mtz_data = f_obs.as_mtz_dataset(column_root_label="F")
  #mtz_data.add_miller_array(
  #  miller_array=f_obs.average_bijvoet_mates(),
  #  column_root_label="F_mean")
  mtz_data.add_miller_array(
    miller_array=fmodel.f_model(),
    column_root_label="FMODEL")
  mtz_data.add_miller_array(
    miller_array=fmodel.f_model().average_bijvoet_mates().phases(deg=True),
    column_root_label="PHI",
    column_types="P")
  mtz_data.add_miller_array(
    miller_array=f_obs_mean.customized_copy(
      data=flex.double(f_obs_mean.data().size(), 0.95),
      sigmas=None).set_observation_type(None),
    column_root_label="FOM",
    column_types="W")
  two_fofc_map = fmodel.map_coefficients(map_type="2mFo-DFc")
  fofc_map = fmodel.map_coefficients(map_type="mFo-Dfc")
  anom_map = fmodel.map_coefficients(map_type="anom")
  mtz_data.add_miller_array(
    miller_array=two_fofc_map.average_bijvoet_mates(),
    column_root_label="2FOFCWT")
  mtz_data.add_miller_array(
    miller_array=fmodel.map_coefficients(map_type="mFo-DFc"),
    column_root_label="FOFCWT")
  mtz_data.add_miller_array(
    miller_array=fmodel.map_coefficients(map_type="anom"),
    column_root_label="ANOM")
  mtz_data.add_miller_array(flags,
    column_root_label="FreeR_flag")
  map_file = "tst_mmtbx_mtz2map_map_coeffs.mtz"
  mtz_data.mtz_object().write(map_file)
  # exercise defaults with PDB file
  file_info = mtz2map.run([pdb_file, map_file], log=null_out())
  file_info = [ (os.path.basename(fn), desc) for fn, desc in file_info ]
  assert (file_info == [
    ('tst_mmtbx_mtz2map_map_coeffs_2mFo-DFc.ccp4', 'CCP4 map'),
    ('tst_mmtbx_mtz2map_map_coeffs_mFo-DFc.ccp4', 'CCP4 map'),
    ('tst_mmtbx_mtz2map_map_coeffs_anom.ccp4', 'CCP4 map'),
    ('tst_mmtbx_mtz2map_map_coeffs_4.ccp4', 'CCP4 map')
  ])
  # without PDB file
  file_info_2 = mtz2map.run([map_file], log=null_out())
  file_info_2 = [ (os.path.basename(fn), desc) for fn, desc in file_info_2 ]
  assert file_info_2 == file_info, file_info_2
  # with FMODEL
  file_info_3 = mtz2map.run([pdb_file, map_file, "include_fmodel=True"],
    log=null_out())
  file_info_3 = [ (os.path.basename(fn), desc) for fn, desc in file_info_3 ]
  assert (file_info_3 == [
    ('tst_mmtbx_mtz2map_map_coeffs_fmodel.ccp4', 'CCP4 map'),
    ('tst_mmtbx_mtz2map_map_coeffs_2mFo-DFc.ccp4', 'CCP4 map'),
    ('tst_mmtbx_mtz2map_map_coeffs_mFo-DFc.ccp4', 'CCP4 map'),
    ('tst_mmtbx_mtz2map_map_coeffs_anom.ccp4', 'CCP4 map'),
    ('tst_mmtbx_mtz2map_map_coeffs_5.ccp4', 'CCP4 map')
  ])
  # exercise bad parameter
  try :
    file_info = mtz2map.run([pdb_file, "1yjp_mtz2map_map_coeffs.mtz",
      "output.directory=1yjp_mtz2map_map_coeffs.mtz"],
      log=null_out())
  except Sorry :
    pass
  else :
    raise Exception_expected
  # bad atom selection
  try :
    file_info = mtz2map.run([pdb_file, map_file, "selection=\"resname ZN\""],
      log=null_out())
  except Sorry as s  :
    assert (str(s) == "No atoms found matching the specified selection.")
  else :
    raise Exception_expected
  # remove R-free flags
  out = StringIO()
  mtz2map.run([pdb_file, "tst_mmtbx_mtz2map_map_coeffs.mtz",
    "r_free_flags.remove=True"], log=out)
  assert (out.getvalue().count("R-free flagged reflections")==4)
  print("OK")

if (__name__ == "__main__"):
  exercise()
