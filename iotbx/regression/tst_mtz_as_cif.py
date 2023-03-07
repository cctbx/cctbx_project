from __future__ import absolute_import, division, print_function
import os
import libtbx.load_env
if (libtbx.env.has_module("ccp4io")):
  from iotbx import mtz
else:
  mtz = None
from libtbx.test_utils import approx_equal
import iotbx.cif

def get_array_by_label(miller_arrays, label):
  for ma in miller_arrays:
    if label in list(ma.info().labels):
      return ma

def exercise():
  if mtz is None:
    print("Skipping iotbx/regression/tst_mtz_as_cif.py: ccp4io not available")
    return
  from iotbx.command_line import mtz_as_cif

  file_name = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/reflection_files/1akg.mtz",
    test=os.path.isfile)
  if file_name is None:
    print("Skipping 1akg.mtz test: input file not available")
  else:
    mtz_as_cif.run(args=[file_name])
    assert os.path.exists("1akg.reflections.cif")
    miller_arrays = iotbx.cif.reader(
      file_path="1akg.reflections.cif").as_miller_arrays()
    mtz_arrays = mtz.object(file_name=file_name).as_miller_arrays()
    assert '_refln.F_meas_au' in miller_arrays[0].info().labels
    assert '_refln.F_meas_sigma_au' in miller_arrays[0].info().labels
    assert approx_equal(
      get_array_by_label(miller_arrays, '_refln.F_meas_au').data(),
      get_array_by_label(mtz_arrays, 'FOBS').data())
    assert approx_equal(
      get_array_by_label(miller_arrays, '_refln.F_meas_au').sigmas(),
      get_array_by_label(mtz_arrays, 'FOBS').sigmas())
    flags1 = get_array_by_label(miller_arrays,'_refln.pdbx_r_free_flag').data()
    flags2 = get_array_by_label(mtz_arrays, 'R-free-flags').data()
    assert flags1.all_eq(flags2)

  file_name = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/reflection_files/1zff.mtz",
    test=os.path.isfile)
  if file_name is None:
    print("Skipping 1zff.mtz test: input file not available")
  else:
    mtz_as_cif.run(args=[file_name])
    assert os.path.exists("1zff.reflections.cif")
    miller_arrays = iotbx.cif.reader(
      file_path="1zff.reflections.cif").as_miller_arrays()
    mtz_arrays = mtz.object(file_name=file_name).as_miller_arrays()
    #assert approx_equal(
      #get_array_by_label(miller_arrays, '_refln.F_calc_au').data(),
      #get_array_by_label(mtz_arrays, 'FC').data())
    assert approx_equal(
      get_array_by_label(miller_arrays, '_refln.F_meas_au').data(),
      get_array_by_label(mtz_arrays, 'F').data(), eps=1e-3)
    assert approx_equal(
      get_array_by_label(miller_arrays, '_refln.F_meas_au').sigmas(),
      get_array_by_label(mtz_arrays, 'F').sigmas(), eps=1e-3)
    assert approx_equal(
      get_array_by_label(miller_arrays, '_refln.pdbx_FWT').data(),
      get_array_by_label(mtz_arrays, '2FOFCWT').data(), eps=1e-2)
    assert approx_equal(
      get_array_by_label(miller_arrays, '_refln.pdbx_DELFWT').data(),
      get_array_by_label(mtz_arrays, 'FOFCWT').data(), eps=1e-3)
    assert approx_equal(
      get_array_by_label(miller_arrays, '_refln.fom').data(),
      get_array_by_label(mtz_arrays, 'FOM').data(), eps=1e-3)

  file_name = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/reflection_files/l.mtz",
    test=os.path.isfile)
  if file_name is None:
    print("Skipping l.mtz test: input file not available")
  else:
    mtz_as_cif.run(args=[file_name])
    assert os.path.exists("l.reflections.cif")
    miller_arrays = iotbx.cif.reader(
      file_path="l.reflections.cif").as_miller_arrays()
    mtz_arrays = mtz.object(file_name=file_name).as_miller_arrays()
    f_anom_cif = get_array_by_label(miller_arrays, '_refln.pdbx_F_plus')
    f_anom_mtz = get_array_by_label(mtz_arrays, 'F-obs(+)')
    f_anom_cif, f_anom_mtz = f_anom_cif.common_sets(
      f_anom_mtz, assert_no_singles=True)
    rfree_cif = get_array_by_label(miller_arrays, '_refln.pdbx_r_free_flag')
    rfree_mtz = get_array_by_label(mtz_arrays, 'R-free-flags(+)')
    assert f_anom_cif.anomalous_flag()
    assert f_anom_mtz.anomalous_flag()
    assert not rfree_cif.anomalous_flag()
    assert rfree_mtz.anomalous_flag()
    assert approx_equal(f_anom_cif.data(), f_anom_mtz.data(), 1e-2)
    assert approx_equal(f_anom_cif.sigmas(), f_anom_mtz.sigmas(), 1e-2)
    assert rfree_mtz.data().all_eq(
      rfree_cif.generate_bijvoet_mates().common_set(rfree_mtz).data())

  file_name = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/reflection_files/outside_4_exptl_fobs_phases_freeR_flags.mtz",
    test=os.path.isfile)
  if file_name is None:
    print("Skipping outside_4_exptl_fobs_phases_freeR_flags.mtz test: input file not available")
  else:
    mtz_as_cif.run(
      args=[file_name,
            'mtz_labels=HLAM HLBM HLCM HLDM',
            'cif_labels=_refln.pdbx_HL_A_iso _refln.pdbx_HL_B_iso _refln.pdbx_HL_C_iso _refln.pdbx_HL_D_iso',
            'output_file=custom_cif_labels.cif'])
    assert os.path.exists("custom_cif_labels.cif")
    miller_arrays = iotbx.cif.reader(
      file_path="custom_cif_labels.cif").as_miller_arrays()
    mtz_arrays = mtz.object(file_name=file_name).as_miller_arrays()
    HL_coeffs_cif = get_array_by_label(miller_arrays, '_refln.pdbx_HL_A_iso')
    HL_coeffs_mtz = get_array_by_label(mtz_arrays, 'HLAM')
    assert HL_coeffs_cif.is_hendrickson_lattman_array()
    assert HL_coeffs_mtz.is_hendrickson_lattman_array()
    assert approx_equal(HL_coeffs_cif.data(), HL_coeffs_mtz.data())

  file_name = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/reflection_files/r1wqzsf.mtz",
    test=os.path.isfile)
  if file_name is None:
    print("Skipping r1wqzsf.mtz test: input file not available")
  else:
    mtz_as_cif.run(
      args=[file_name])
    assert os.path.exists("r1wqzsf.reflections.cif")
    miller_arrays = iotbx.cif.reader(
      file_path="r1wqzsf.reflections.cif").as_miller_arrays()
    mtz_arrays = mtz.object(file_name=file_name).as_miller_arrays()
    fobs_cif = get_array_by_label(miller_arrays, '_refln.F_meas_au')
    fobs_mtz = get_array_by_label(mtz_arrays, 'FOBS_N')
    assert approx_equal(fobs_cif.data(), fobs_mtz.data())
    assert approx_equal(fobs_cif.sigmas(), fobs_mtz.sigmas())
    rfree_cif = get_array_by_label(miller_arrays, '_refln.pdbx_r_free_flag')
    rfree_mtz = get_array_by_label(mtz_arrays, 'R-free-flags')
    assert rfree_cif.data().all_eq(rfree_mtz.data())

  file_name = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/reflection_files/ur0013.sf.mtz",
    test=os.path.isfile)
  if file_name is None:
    print("Skipping ur0013.sf.mtz test: input file not available")
  else:
    mtz_as_cif.run(
      args=[file_name])
    assert os.path.exists("ur0013.sf.reflections.cif")
    cif_reader = iotbx.cif.reader(
      file_path="ur0013.sf.reflections.cif")
    miller_arrays = cif_reader.as_miller_arrays()
    cif_object = cif_reader.model()
    assert list(cif_object.keys()) == ['ur0013.sf_neutron']
    mtz_arrays = mtz.object(file_name=file_name).as_miller_arrays()
    iobs_cif = get_array_by_label(miller_arrays, '_refln.pdbx_I_plus')
    iobs_mtz = get_array_by_label(mtz_arrays, 'IOBS_N(+)')
    assert iobs_cif.anomalous_flag()
    assert iobs_mtz.anomalous_flag()
    iobs_cif, iobs_mtz = iobs_cif.common_sets(iobs_mtz, assert_no_singles=True)
    assert approx_equal(iobs_cif.data(), iobs_mtz.data())
    assert approx_equal(iobs_cif.sigmas(), iobs_mtz.sigmas())

  file_name = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/reflection_files/1aba_test.mtz",
    test=os.path.isfile)
  if file_name is None:
    print("Skipping 1aba_test.mtz test: input file not available")
  else:
    mtz_as_cif.run(
      args=[file_name])

def run():
  exercise()
  print("OK")

if __name__ == '__main__':
  run()
