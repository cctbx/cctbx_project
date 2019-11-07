from __future__ import absolute_import, division, print_function
import iotbx.pdb
from cctbx import crystal
from libtbx.test_utils import approx_equal
import mmtbx.model
from six.moves import range

cif_str="""
data_5JUP
#
_atom_sites.entry_id                    5JUP
_atom_sites.fract_transf_matrix[1][1]   1.000000
_atom_sites.fract_transf_matrix[1][2]   0.000000
_atom_sites.fract_transf_matrix[1][3]   0.000000
_atom_sites.fract_transf_matrix[2][1]   0.000000
_atom_sites.fract_transf_matrix[2][2]   1.000000
_atom_sites.fract_transf_matrix[2][3]   0.000000
_atom_sites.fract_transf_matrix[3][1]   0.000000
_atom_sites.fract_transf_matrix[3][2]   0.000000
_atom_sites.fract_transf_matrix[3][3]   1.000000
_atom_sites.fract_transf_vector[1]      0.00000
_atom_sites.fract_transf_vector[2]      0.00000
_atom_sites.fract_transf_vector[3]      0.00000
#
loop_
_atom_type.symbol
C
MG
N
O
P
S
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
_atom_site.auth_seq_id
_atom_site.auth_comp_id
_atom_site.auth_asym_id
_atom_site.auth_atom_id
_atom_site.pdbx_PDB_model_num
ATOM   9      C  "C3'" . U   A  1  1    ? 175.229 258.091 229.127 1.00 113.67 1    U   A  "C3'" 1
ATOM   7      C  "C4'" . U   A  1  1    ? 175.388 257.913 227.616 1.00 113.21 1    U   A  "C4'" 1
ATOM   6      C  "C5'" . U   A  1  1    ? 174.124 257.916 226.769 1.00 114.11 1    U   A  "C5'" 1
ATOM   10     O  "O3'" . U   A  1  1    ? 176.253 257.344 229.761 1.00 112.89 1    U   A  "O3'" 1
ATOM   8      O  "O4'" . U   A  1  1    ? 176.291 258.974 227.185 1.00 112.70 1    U   A  "O4'" 1
ATOM   5      O  "O5'" . U   A  1  1    ? 173.735 259.270 226.430 1.00 115.09 1    U   A  "O5'" 1
ATOM   3      O  OP1   . U   A  1  1    ? 174.056 259.302 223.918 1.00 117.16 1    U   A  OP1   1
ATOM   4      O  OP2   . U   A  1  1    ? 172.661 261.090 225.058 1.00 104.87 1    U   A  OP2   1
ATOM   1      O  OP3   . U   A  1  1    ? 171.796 258.709 224.899 1.00 116.23 1    U   A  OP3   1
ATOM   2      P  P     . U   A  1  1    ? 173.022 259.610 224.994 1.00 116.34 1    U   A  P     1
#
"""

def run(prefix="tst_iotbx_tst_xray_scale_1"):
  """
  Exercise obtaining exactly the same xray structure from pdb and cif.
  Caveate is SCALE records (fract_transf_matrix). In pdb they explicitly
  suppressed when crystal_symmetry is supplied, see
  iotbx/pdb/__init__.py:939 (def xray_structure_simple):
  if(crystal_symmetry is not None): self._scale_matrix = None
  """
  fo = open("%s.cif" % prefix,"w")
  print(cif_str, file=fo)
  fo.close()
  # crystal symmetry from map
  cs = crystal.symmetry((419., 419., 419., 90.0, 90.0, 90.0), 1)
  # xrs from mmCIF
  cif_inp = iotbx.pdb.input(file_name="%s.cif" % prefix)
  # print "======== doing xrs from cif =============="
  xrs_cif = cif_inp.xray_structure_simple(crystal_symmetry = cs)
  # xrs from PDB
  cif_inp.construct_hierarchy().write_pdb_file(prefix+".pdb")
  pdb_inp = iotbx.pdb.input(file_name="%s.pdb" % prefix)
  # print "======== doing xrs from pdb =============="
  xrs_pdb = pdb_inp.xray_structure_simple(crystal_symmetry = cs)
  #
  for i in range(len(xrs_cif.sites_cart())):
    assert approx_equal(xrs_cif.sites_cart()[i], xrs_pdb.sites_cart()[i])
  print("OK")

def run2(prefix="iotbx_tst_xray_scale_2"):
  """ same as run(), but using mmtbx.model.
  Not clear why it is failing... Need to ask Pavel to look together."""
  cs = crystal.symmetry((419., 419., 419., 90.0, 90.0, 90.0), 1)
  model_cif = mmtbx.model.manager(
      model_input = iotbx.pdb.input(lines=cif_str.split("\n"), source_info=None),
      crystal_symmetry=cs)
  xrs_cif = model_cif.get_xray_structure()
  txt_pdb = model_cif.model_as_pdb()
  model_pdb = mmtbx.model.manager(
      model_input = iotbx.pdb.input(lines=txt_pdb.split("\n"), source_info=None),
      crystal_symmetry=cs)
  xrs_pdb = model_pdb.get_xray_structure()
  for i in range(len(xrs_cif.sites_cart())):
    assert approx_equal(xrs_cif.sites_cart()[i], xrs_pdb.sites_cart()[i])
  print("OK")

if (__name__ == "__main__"):
  run()
  run2()
