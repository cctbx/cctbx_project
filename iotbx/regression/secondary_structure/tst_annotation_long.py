from __future__ import absolute_import, division, print_function
import time
from iotbx.pdb.secondary_structure import annotation
import iotbx
import iotbx.cif
from libtbx.test_utils import show_diff

cif_one_helix = """\
#
_struct_conf.conf_type_id            HELX_P
_struct_conf.id                      HELX_P1
_struct_conf.pdbx_PDB_helix_id       AA1
_struct_conf.beg_label_comp_id       longGLN
_struct_conf.beg_label_asym_id       longA
_struct_conf.beg_label_seq_id        149
_struct_conf.pdbx_beg_PDB_ins_code   ?
_struct_conf.end_label_comp_id       longGLY
_struct_conf.end_label_asym_id       longA
_struct_conf.end_label_seq_id        160
_struct_conf.pdbx_end_PDB_ins_code   ?
_struct_conf.beg_auth_comp_id        alongGLN
_struct_conf.beg_auth_asym_id        alongA
_struct_conf.beg_auth_seq_id         212
_struct_conf.end_auth_comp_id        alongGLY
_struct_conf.end_auth_asym_id        alongA
_struct_conf.end_auth_seq_id         223
_struct_conf.pdbx_PDB_helix_class    1
_struct_conf.details                 ?
_struct_conf.pdbx_PDB_helix_length   12
#
_struct_conf_type.id          HELX_P
_struct_conf_type.criteria    ?
_struct_conf_type.reference   ?
"""

def tst_helix_interface():
  cif_str = "data_1UCS" + cif_one_helix
  cif_model = iotbx.cif.reader(input_string=cif_str).model()
  cif_block = list(cif_model.values())[0]
  ann = annotation.from_cif_block(cif_block)
  h = ann.helices[0]
  #
  # !!! Note, here we got info from *auth* fields.
  #
  assert h.serial == 1
  assert h.helix_id == 'AA1', h.helix_id
  assert h.start_resname == 'alongGLN', h.start_resname
  assert h.start_chain_id == 'alongA'
  assert h.start_resseq == ' 212', h.start_resseq
  assert h.start_icode == ' '
  assert h.end_resname == 'alongGLY'
  assert h.end_chain_id == 'alongA'
  assert h.end_resseq == ' 223'
  assert h.end_icode == ' '
  assert h.helix_class == 'alpha'
  assert h.comment == '', "'%s'" % h.comment
  assert h.length == 12, h.length
  assert h.hbond_list == []
  assert h.helix_selection == None
  assert h.enabled == True
  assert h.sigma ==0.05
  assert h.slack ==0
  assert h.top_out == False

  assert h.get_start_resseq_as_int() == 212
  assert h.get_end_resseq_as_int() == 223
  assert not h.fits_in_pdb_format()
  assert not ann.fits_in_pdb_format()
  # print (h.as_pdb_or_mmcif_str())
  assert not show_diff(h.as_mmcif_str(),
"""data_phenix
loop_
  _struct_conf.conf_type_id
  _struct_conf.id
  _struct_conf.pdbx_PDB_helix_id
  _struct_conf.beg_label_comp_id
  _struct_conf.beg_label_asym_id
  _struct_conf.beg_label_seq_id
  _struct_conf.pdbx_beg_PDB_ins_code
  _struct_conf.end_label_comp_id
  _struct_conf.end_label_asym_id
  _struct_conf.end_label_seq_id
  _struct_conf.pdbx_end_PDB_ins_code
  _struct_conf.pdbx_PDB_helix_class
  _struct_conf.details
  _struct_conf.pdbx_PDB_helix_length
  HELX_P  1  AA1  alongGLN  alongA  212  ?  alongGLY  alongA  223  ?  1  ?  12

loop_
  _struct_conf_type.id
  _struct_conf_type.criteria
  _struct_conf_type.reference
  HELX_P  ?  ?

""")

if (__name__ == "__main__"):
  t0 = time.time()
  tst_helix_interface()
  print("OK time =%8.3f"%(time.time() - t0))
