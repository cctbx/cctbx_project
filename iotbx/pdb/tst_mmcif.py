"""Test working with mmCIF, calculating fmodel, reading headers"""
from __future__ import absolute_import, division, print_function
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

def exercise_extract_header_misc():
  cif_file = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/pdb/3orl.cif",
    test=os.path.isfile)
  if (cif_file is None):
    return
  cif_in = mmcif.cif_input(file_name=cif_file)
  assert (cif_in.file_type() == "mmcif")
  wavelength = cif_in.extract_wavelength()
  assert (approx_equal(wavelength, 1.8927))
  exptl_method = cif_in.get_experiment_type()
  assert exptl_method.is_xray()
  assert ("%s" % exptl_method == 'X-RAY DIFFRACTION')
  r_rfree_sigma = cif_in.get_r_rfree_sigma(cif_file)
  r_rfree_sigma.show()

def exercise_rows_splitted_by_newline():
  cif_txt = """
data_myblock
#
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_entity_id
_atom_site.label_seq_id
_atom_site.pdbx_PDB_ins_code
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
_atom_site.pdbx_formal_charge
_atom_site.auth_seq_id
_atom_site.auth_comp_id
_atom_site.auth_asym_id
_atom_site.auth_atom_id
_atom_site.pdbx_PDB_model_num
ATOM   1    N  N   ? MET A 1 1   ? 34.204 -29.806 -27.207 1.00 48.27 ? 0   MET
A N   1
ATOM   2    C  CA  ? MET A 1 1   ? 34.593 -28.347 -27.188 1.00 47.97 ? 0   MET
A CA  1
ATOM   3    C  C   ? MET A 1 1   ? 34.888 -27.889 -25.763 1.00 46.15 ? 0   MET
A C   1
ATOM   4    O  O   ? MET A 1 1   ? 35.022 -28.720 -24.852 1.00 46.47 ? 0   MET
A O   1
ATOM   5    C  CB  ? MET A 1 1   ? 33.519 -27.457 -27.852 1.00 48.41 ? 0   MET
A CB  1
ATOM   6    C  CG  ? MET A 1 1   ? 32.163 -27.412 -27.138 1.00 49.22 ? 0   MET
A CG  1
ATOM   7    S  SD  ? MET A 1 1   ? 31.161 -25.921 -27.465 1.00 50.44 ? 0   MET
A SD  1
ATOM   8    C  CE  ? MET A 1 1   ? 32.261 -24.545 -27.108 1.00 50.23 ? 0   MET
A CE  1
#
"""
  cif_block = iotbx.cif.reader(input_string=cif_txt).model()["myblock"]
  assert list(cif_block.get("_atom_site.group_PDB")) == ["ATOM"]*8, list(cif_block.get("_atom_site.group_PDB"))
  assert list(cif_block.get("_atom_site.id")) == ['1','2','3','4','5','6','7','8'], list(cif_block.get("_atom_site.id"))
  assert list(cif_block.get("_atom_site.type_symbol")) == ["N", "C", "C", "O", "C", "C", "S", "C"]
  assert list(cif_block.get("_atom_site.label_atom_id")) == ["N", "CA", "C", "O", "CB", "CG", "SD", "CE"]
  assert list(cif_block.get("_atom_site.label_alt_id")) == ["?"]*8
  assert list(cif_block.get("_atom_site.label_comp_id")) == ["MET"]*8
  assert list(cif_block.get("_atom_site.label_asym_id")) == ["A"]*8
  assert list(cif_block.get("_atom_site.label_entity_id")) == ["1"]*8
  assert list(cif_block.get("_atom_site.label_seq_id")) == ["1"]*8
  assert list(cif_block.get("_atom_site.pdbx_PDB_ins_code")) == ["?"]*8
  # assert list(cif_block.get("_atom_site.Cartn_x")) == []
  # assert list(cif_block.get("_atom_site.Cartn_y")) == []
  # assert list(cif_block.get("_atom_site.Cartn_z")) == []
  assert list(cif_block.get("_atom_site.occupancy")) == ['1.00']*8, list(cif_block.get("_atom_site.occupancy"))
  # assert list(cif_block.get("_atom_site.B_iso_or_equiv")) == []
  assert list(cif_block.get("_atom_site.pdbx_formal_charge")) == ["?"]*8
  assert list(cif_block.get("_atom_site.auth_seq_id")) == ["0"]*8
  assert list(cif_block.get("_atom_site.auth_comp_id")) == ["MET"]*8
  assert list(cif_block.get("_atom_site.auth_asym_id")) == ["A"]*8
  assert list(cif_block.get("_atom_site.auth_atom_id")) == ["N", "CA", "C", "O", "CB", "CG", "SD", "CE"]
  assert list(cif_block.get("_atom_site.pdbx_PDB_model_num")) == ["1"]*8

def run():
  exercise_extract_f_model_core_constants()
  exercise_extract_header_misc()
  exercise_rows_splitted_by_newline()

if __name__ == '__main__':
  run()
  print("OK")
