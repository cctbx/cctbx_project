from __future__ import division
from __future__ import print_function
import iotbx.cif
from iotbx.pdb import mmcif
from libtbx.test_utils import approx_equal
import libtbx.load_env
import os

def exercise_extract_f_model_core_constants():
  cif_input = """\
data_test
_refine.solvent_model_param_bsol               22.0810
_refine.solvent_model_param_ksol               0.3340
_refine.pdbx_solvent_vdw_probe_radii           1.2000
_refine.pdbx_solvent_ion_probe_radii           ?
_refine.pdbx_solvent_shrinkage_radii           0.9500
_refine.aniso_B[1][1]                          2.7593
_refine.aniso_B[2][2]                          2.7593
_refine.aniso_B[3][3]                          -5.5186
_refine.aniso_B[1][2]                          0.0000
_refine.aniso_B[1][3]                          -0.0000
_refine.aniso_B[2][3]                          0.0000
_refine.ls_R_factor_R_work                     0.1690
_refine.ls_R_factor_R_free                     0.2238
"""
  cif_block = iotbx.cif.reader(input_string=cif_input).model()["test"]
  constants = mmcif.extract_f_model_core_constants(cif_block)
  assert approx_equal(constants.k_sol, 0.3340)
  assert approx_equal(constants.b_sol, 22.0810)
  assert approx_equal(constants.r_solv, 1.2)
  assert approx_equal(constants.r_shrink, 0.95)
  assert approx_equal(constants.b_cart, (2.7593, 2.7593, -5.5186, 0, 0, 0,))
  assert approx_equal(constants.r_work, 0.1690)
  assert approx_equal(constants.r_free, 0.2238)

def exercise_extract_header_misc () :
  cif_file = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/pdb/3orl.cif",
    test=os.path.isfile)
  if (cif_file is None) :
    return
  cif_in = mmcif.cif_input(file_name=cif_file)
  assert (cif_in.file_type() == "mmcif")
  wavelength = cif_in.extract_wavelength()
  assert (approx_equal(wavelength, 1.8927))
  exptl_method = cif_in.get_experiment_type()
  assert (exptl_method == 'X-RAY DIFFRACTION')
  r_rfree_sigma = cif_in.get_r_rfree_sigma(cif_file)
  r_rfree_sigma.show()


def run():
  exercise_extract_f_model_core_constants()
  exercise_extract_header_misc()

if __name__ == '__main__':
  run()
  print("OK")
