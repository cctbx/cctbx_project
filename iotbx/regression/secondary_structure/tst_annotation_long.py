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

cif_one_sheet = """\
#
_struct_sheet.id               A
_struct_sheet.type             ?
_struct_sheet.number_strands   2
_struct_sheet.details          ?
#
_struct_sheet_order.sheet_id     A
_struct_sheet_order.range_id_1   1
_struct_sheet_order.range_id_2   2
_struct_sheet_order.offset       ?
_struct_sheet_order.sense        anti-parallel
#
loop_
_struct_sheet_range.sheet_id
_struct_sheet_range.id
_struct_sheet_range.beg_label_comp_id
_struct_sheet_range.beg_label_asym_id
_struct_sheet_range.beg_label_seq_id
_struct_sheet_range.pdbx_beg_PDB_ins_code
_struct_sheet_range.end_label_comp_id
_struct_sheet_range.end_label_asym_id
_struct_sheet_range.end_label_seq_id
_struct_sheet_range.pdbx_end_PDB_ins_code
_struct_sheet_range.symmetry
_struct_sheet_range.beg_auth_comp_id
_struct_sheet_range.beg_auth_asym_id
_struct_sheet_range.beg_auth_seq_id
_struct_sheet_range.end_auth_comp_id
_struct_sheet_range.end_auth_asym_id
_struct_sheet_range.end_auth_seq_id
A 1 longSER longA 4  ? longALA longA 7  ? ? alongSER alongA 4  alongALA alongA 7
A 2 longMET longA 22 ? longGLU longA 25 ? ? alongMET alongA 22 alongGLU alongA 25
#
_pdbx_struct_sheet_hbond.sheet_id                A
_pdbx_struct_sheet_hbond.range_id_1              1
_pdbx_struct_sheet_hbond.range_id_2              2
_pdbx_struct_sheet_hbond.range_1_label_atom_id   N
_pdbx_struct_sheet_hbond.range_1_label_comp_id   longSER
_pdbx_struct_sheet_hbond.range_1_label_asym_id   longA
_pdbx_struct_sheet_hbond.range_1_label_seq_id    4
_pdbx_struct_sheet_hbond.range_1_PDB_ins_code    ?
_pdbx_struct_sheet_hbond.range_1_auth_atom_id    N
_pdbx_struct_sheet_hbond.range_1_auth_comp_id    alongSER
_pdbx_struct_sheet_hbond.range_1_auth_asym_id    alongA
_pdbx_struct_sheet_hbond.range_1_auth_seq_id     4
_pdbx_struct_sheet_hbond.range_2_label_atom_id   O
_pdbx_struct_sheet_hbond.range_2_label_comp_id   longGLU
_pdbx_struct_sheet_hbond.range_2_label_asym_id   longA
_pdbx_struct_sheet_hbond.range_2_label_seq_id    25
_pdbx_struct_sheet_hbond.range_2_PDB_ins_code    ?
_pdbx_struct_sheet_hbond.range_2_auth_atom_id    O
_pdbx_struct_sheet_hbond.range_2_auth_comp_id    alongGLU
_pdbx_struct_sheet_hbond.range_2_auth_asym_id    alongA
_pdbx_struct_sheet_hbond.range_2_auth_seq_id     25
#
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

def tst_sheet_interface():
  cif_str = "data_1UCS" + cif_one_sheet
  cif_model = iotbx.cif.reader(input_string=cif_str).model()
  cif_block = list(cif_model.values())[0]
  ann = annotation.from_cif_block(cif_block)

  sh = ann.sheets[0]
  st = sh.strands[1]
  reg = sh.registrations[1]

  #
  # !!! Note, here we got info from *auth* fields.
  #

  assert reg.cur_atom == ' O  '
  assert reg.cur_resname == 'alongGLU', reg.cur_resname
  assert reg.cur_chain_id == 'alongA'
  assert reg.cur_resseq == '  25'
  assert reg.cur_icode == ' '
  assert reg.prev_atom == ' N  '
  assert reg.prev_resname == 'alongSER'
  assert reg.prev_chain_id == 'alongA'
  assert reg.prev_resseq == '   4'
  assert reg.prev_icode  == ' '

  assert st.sheet_id == 'A', st.sheet_id
  assert st.strand_id == 2, st.strand_id
  assert st.start_resname == 'alongMET',st.start_resname
  assert st.start_chain_id == 'alongA'
  assert st.start_resseq == '  22', st.start_resseq
  assert st.start_icode == ' '
  assert st.end_resname == 'alongGLU'
  assert st.end_chain_id == 'alongA'
  assert st.end_resseq == '  25'
  assert st.end_icode == ' '
  assert st.sense == -1

  assert not sh.fits_in_pdb_format()
  assert not ann.fits_in_pdb_format()
  # print (sh.as_pdb_or_mmcif_str())
  assert not show_diff(sh.as_mmcif_str(),
"""data_phenix
loop_
  _struct_sheet.id
  _struct_sheet.type
  _struct_sheet.number_strands
  _struct_sheet.details
  A  ?  2  ?

loop_
  _struct_sheet_order.sheet_id
  _struct_sheet_order.range_id_1
  _struct_sheet_order.range_id_2
  _struct_sheet_order.offset
  _struct_sheet_order.sense
  A  1  2  ?  anti-parallel

loop_
  _struct_sheet_range.sheet_id
  _struct_sheet_range.id
  _struct_sheet_range.beg_label_comp_id
  _struct_sheet_range.beg_label_asym_id
  _struct_sheet_range.beg_label_seq_id
  _struct_sheet_range.pdbx_beg_PDB_ins_code
  _struct_sheet_range.end_label_comp_id
  _struct_sheet_range.end_label_asym_id
  _struct_sheet_range.end_label_seq_id
  _struct_sheet_range.pdbx_end_PDB_ins_code
  A  1  alongSER  alongA   4  ?  alongALA  alongA   7  ?
  A  2  alongMET  alongA  22  ?  alongGLU  alongA  25  ?

loop_
  _pdbx_struct_sheet_hbond.sheet_id
  _pdbx_struct_sheet_hbond.range_id_1
  _pdbx_struct_sheet_hbond.range_id_2
  _pdbx_struct_sheet_hbond.range_1_label_atom_id
  _pdbx_struct_sheet_hbond.range_1_label_comp_id
  _pdbx_struct_sheet_hbond.range_1_label_asym_id
  _pdbx_struct_sheet_hbond.range_1_label_seq_id
  _pdbx_struct_sheet_hbond.range_1_PDB_ins_code
  _pdbx_struct_sheet_hbond.range_2_label_atom_id
  _pdbx_struct_sheet_hbond.range_2_label_comp_id
  _pdbx_struct_sheet_hbond.range_2_label_asym_id
  _pdbx_struct_sheet_hbond.range_2_label_seq_id
  _pdbx_struct_sheet_hbond.range_2_PDB_ins_code
  A  1  2  N  alongSER  alongA  4  ?  O  alongGLU  alongA  25  ?

""")

if (__name__ == "__main__"):
  t0 = time.time()
  tst_helix_interface()
  tst_sheet_interface()
  print("OK time =%8.3f"%(time.time() - t0))
