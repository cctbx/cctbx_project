from __future__ import absolute_import, division, print_function
import time
from iotbx.pdb.secondary_structure import annotation, pdb_helix, pdb_strand
import iotbx
import iotbx.cif
from libtbx.test_utils import show_diff
from six.moves import cStringIO as StringIO
import libtbx.load_env
import os
from six.moves import range

def test_helix_interface():
  # helix_class_to_int.
  h_class_array = ['unknown','alpha', 'unknown', 'pi', 'unknown',
        '3_10', 'unknown', 'unknown', 'unknown', 'unknown', 'unknown']
  for n, hc in enumerate(h_class_array):
    # print pdb_helix.helix_class_to_int(hc), n
    if hc == 'unknown':
      assert pdb_helix.helix_class_to_int(hc) == 0
    else:
      assert pdb_helix.helix_class_to_int(hc) == n
    assert pdb_helix.helix_class_to_str(n) == hc

  h_string = \
"HELIX    1   1 ALA A   16  THR A   18  1                                   3"
  h_string_trunc = "HELIX    1   1 ALA A   16  THR A   18  1"
  got_exception = False
  try:
    h = pdb_helix.from_pdb_record(h_string_trunc)
  except ValueError:
    got_exception = True
  assert got_exception
  h = pdb_helix.from_pdb_record(h_string)
  assert h.serial == '  1'
  assert h.helix_id == '1'
  assert h.start_resname == 'ALA'
  assert h.start_chain_id == 'A'
  assert h.start_resseq == '  16'
  assert h.start_icode == ' '
  assert h.end_resname == 'THR'
  assert h.end_chain_id == 'A'
  assert h.end_resseq == '  18'
  assert h.end_icode == ' '
  assert h.helix_class == 'alpha'
  assert h.comment == '                              '
  assert h.length == 3
  assert h.hbond_list == []
  assert h.helix_selection == None
  assert h.enabled == True
  assert h.sigma ==0.05
  assert h.slack ==0
  assert h.top_out == False

  assert h.get_start_resseq_as_int() == 16
  assert h.get_end_resseq_as_int() == 18

  for resseq in ['15', ' 15', '  15', 15]:
    h.set_start_resseq(resseq)
    assert h.start_resseq == '  15'
    assert h.get_start_resseq_as_int() == 15
  for resseq in ['15', ' 15', '  15', 15]:
    h.set_end_resseq(resseq)
    assert h.end_resseq == '  15'
    assert h.get_end_resseq_as_int() == 15
  for resseq in ['AA15', 23001]:
    h.set_start_resseq(resseq)
    assert h.start_resseq == 'AA15'
    assert h.get_start_resseq_as_int() == 23001
  for resseq in ['AA19', 23005]:
    h.set_end_resseq(resseq)
    assert h.end_resseq == 'AA19'
    assert h.get_end_resseq_as_int() == 23005


def test_sheet_interface():
  sheet_str = """\
SHEET    1   A 2 ARG A  13  ASP A  14  0
SHEET    2   A 2 LEU A  27  SER A  30 -1  O  ARG A  29   N  ARG A  13
  """
  ann = annotation.from_records(sheet_str.split("\n"))
  sh = ann.sheets[0]
  st = sh.strands[1]
  reg = sh.registrations[1]
  assert reg.cur_atom == ' O  '
  assert reg.cur_resname == 'ARG'
  assert reg.cur_chain_id == 'A'
  assert reg.cur_resseq == '  29'
  assert reg.cur_icode == ' '
  assert reg.prev_atom == ' N  '
  assert reg.prev_resname == 'ARG'
  assert reg.prev_chain_id == 'A'
  assert reg.prev_resseq == '  13'
  assert reg.prev_icode  == ' '

  assert st.sheet_id == '  A', st.sheet_id
  assert st.strand_id == '  2', st.strand_id
  assert st.start_resname == 'LEU'
  assert st.start_chain_id == 'A'
  assert st.start_resseq == '  27', st.start_resseq
  assert st.start_icode == ' '
  assert st.end_resname == 'SER'
  assert st.end_chain_id == 'A'
  assert st.end_resseq == '  30'
  assert st.end_icode == ' '
  assert st.sense == -1

  for resseq in ['15', ' 15', '  15', 15]:
    reg.set_cur_resseq(resseq)
    assert reg.cur_resseq == '  15'
    assert reg.get_cur_resseq_as_int() == 15
  for resseq in ['AA15', 23001]:
    reg.set_cur_resseq(resseq)
    assert reg.cur_resseq == 'AA15'
    assert reg.get_cur_resseq_as_int() == 23001

  for resseq in ['15', ' 15', '  15', 15]:
    reg.set_prev_resseq(resseq)
    assert reg.prev_resseq == '  15'
    assert reg.get_prev_resseq_as_int() == 15
  for resseq in ['AA15', 23001]:
    reg.set_prev_resseq(resseq)
    assert reg.prev_resseq == 'AA15'
    assert reg.get_prev_resseq_as_int() == 23001

  for resseq in ['15', ' 15', '  15', 15]:
    st.set_start_resseq(resseq)
    assert st.start_resseq == '  15'
    assert st.get_start_resseq_as_int() == 15
  for resseq in ['AA15', 23001]:
    st.set_start_resseq(resseq)
    assert st.start_resseq == 'AA15'
    assert st.get_start_resseq_as_int() == 23001

  for resseq in ['15', ' 15', '  15', 15]:
    st.set_end_resseq(resseq)
    assert st.end_resseq == '  15'
    assert st.get_end_resseq_as_int() == 15
  for resseq in ['AA15', 23001]:
    st.set_end_resseq(resseq)
    assert st.end_resseq == 'AA15'
    assert st.get_end_resseq_as_int() == 23001

cif_one_helix = """\
#
_struct_conf.conf_type_id            HELX_P
_struct_conf.id                      HELX_P1
_struct_conf.pdbx_PDB_helix_id       AA1
_struct_conf.beg_label_comp_id       GLN
_struct_conf.beg_label_asym_id       A
_struct_conf.beg_label_seq_id        149
_struct_conf.pdbx_beg_PDB_ins_code   ?
_struct_conf.end_label_comp_id       GLY
_struct_conf.end_label_asym_id       A
_struct_conf.end_label_seq_id        160
_struct_conf.pdbx_end_PDB_ins_code   ?
_struct_conf.beg_auth_comp_id        GLN
_struct_conf.beg_auth_asym_id        A
_struct_conf.beg_auth_seq_id         212
_struct_conf.end_auth_comp_id        GLY
_struct_conf.end_auth_asym_id        A
_struct_conf.end_auth_seq_id         223
_struct_conf.pdbx_PDB_helix_class    1
_struct_conf.details                 ?
_struct_conf.pdbx_PDB_helix_length   12
#
_struct_conf_type.id          HELX_P
_struct_conf_type.criteria    ?
_struct_conf_type.reference   ?
"""

cif_4_helices = """\
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
_struct_conf.beg_auth_comp_id
_struct_conf.beg_auth_asym_id
_struct_conf.beg_auth_seq_id
_struct_conf.end_auth_comp_id
_struct_conf.end_auth_asym_id
_struct_conf.end_auth_seq_id
_struct_conf.pdbx_PDB_helix_class
_struct_conf.details
_struct_conf.pdbx_PDB_helix_length
HELX_P HELX_P1 AA1 SER A 24  ? LYS A 28  ? SER A 26  LYS A 30  5 ? 5
HELX_P HELX_P2 AA2 LEU A 194 ? GLU A 197 ? LEU A 196 GLU A 199 5 ? 4
HELX_P HELX_P3 AA3 SER B 24  ? LYS B 28  ? SER B 26  LYS B 30  5 ? 5
HELX_P HELX_P4 AA4 LEU B 194 ? GLU B 197 ? LEU B 196 GLU B 199 5 ? 4
#
"""

cif_many_sheets = """\
#
loop_
_struct_sheet.id
_struct_sheet.type
_struct_sheet.number_strands
_struct_sheet.details
AA1 ? 4 ?
AA2 ? 3 ?
AA3 ? 2 ?
AA4 ? 4 ?
AA5 ? 3 ?
AA6 ? 4 ?
AA7 ? 3 ?
AA8 ? 2 ?
AA9 ? 4 ?
AB1 ? 3 ?
#
loop_
_struct_sheet_order.sheet_id
_struct_sheet_order.range_id_1
_struct_sheet_order.range_id_2
_struct_sheet_order.offset
_struct_sheet_order.sense
AA1 1 2 ? parallel
AA1 2 3 ? anti-parallel
AA1 3 4 ? anti-parallel
AA2 1 2 ? anti-parallel
AA2 2 3 ? anti-parallel
AA3 1 2 ? anti-parallel
AA4 1 2 ? parallel
AA4 2 3 ? anti-parallel
AA4 3 4 ? anti-parallel
AA5 1 2 ? anti-parallel
AA5 2 3 ? anti-parallel
AA6 1 2 ? parallel
AA6 2 3 ? anti-parallel
AA6 3 4 ? anti-parallel
AA7 1 2 ? anti-parallel
AA7 2 3 ? anti-parallel
AA8 1 2 ? anti-parallel
AA9 1 2 ? parallel
AA9 2 3 ? anti-parallel
AA9 3 4 ? anti-parallel
AB1 1 2 ? anti-parallel
AB1 2 3 ? anti-parallel
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
_struct_sheet_range.beg_auth_comp_id
_struct_sheet_range.beg_auth_asym_id
_struct_sheet_range.beg_auth_seq_id
_struct_sheet_range.end_auth_comp_id
_struct_sheet_range.end_auth_asym_id
_struct_sheet_range.end_auth_seq_id
AA1 1 ILE A 5   ? PRO A 8   ? ILE A 7   PRO A 10
AA1 2 MET A 90  ? THR A 97  ? MET A 92  THR A 99
AA1 3 THR A 71  ? SER A 80  ? THR A 73  SER A 82
AA1 4 VAL A 32  ? THR A 37  ? VAL A 34  THR A 39
AA2 1 LYS A 17  ? GLN A 21  ? LYS A 19  GLN A 23
AA2 2 TRP A 57  ? VAL A 60  ? TRP A 59  VAL A 62
AA2 3 PHE A 49  ? ILE A 51  ? PHE A 51  ILE A 53
AA3 1 GLU A 105 ? PHE A 106 ? GLU A 107 PHE A 108
AA3 2 ALA A 130 ? THR A 131 ? ALA A 132 THR A 133
AA4 1 VAL A 110 ? MET A 116 ? VAL A 112 MET A 118
AA4 2 SER A 200 ? THR A 210 ? SER A 202 THR A 212
AA4 3 THR A 184 ? ALA A 192 ? THR A 186 ALA A 194
AA4 4 ALA A 145 ? ASP A 152 ? ALA A 147 ASP A 154
AA5 1 THR A 123 ? GLU A 127 ? THR A 125 GLU A 129
AA5 2 VAL A 169 ? VAL A 172 ? VAL A 171 VAL A 174
AA5 3 PHE A 161 ? ILE A 163 ? PHE A 163 ILE A 165
AA6 1 ILE B 5   ? PRO B 8   ? ILE B 7   PRO B 10
AA6 2 MET B 90  ? THR B 97  ? MET B 92  THR B 99
AA6 3 THR B 71  ? SER B 80  ? THR B 73  SER B 82
AA6 4 VAL B 32  ? GLY B 38  ? VAL B 34  GLY B 40
AA7 1 LYS B 17  ? GLN B 21  ? LYS B 19  GLN B 23
AA7 2 TRP B 57  ? VAL B 60  ? TRP B 59  VAL B 62
AA7 3 PHE B 49  ? ILE B 51  ? PHE B 51  ILE B 53
AA8 1 GLU B 105 ? PHE B 106 ? GLU B 107 PHE B 108
AA8 2 ALA B 130 ? THR B 131 ? ALA B 132 THR B 133
AA9 1 VAL B 110 ? MET B 116 ? VAL B 112 MET B 118
AA9 2 SER B 200 ? THR B 210 ? SER B 202 THR B 212
AA9 3 THR B 184 ? ALA B 192 ? THR B 186 ALA B 194
AA9 4 ALA B 145 ? ASP B 152 ? ALA B 147 ASP B 154
AB1 1 SER B 124 ? GLU B 127 ? SER B 126 GLU B 129
AB1 2 VAL B 169 ? VAL B 172 ? VAL B 171 VAL B 174
AB1 3 PHE B 161 ? ILE B 163 ? PHE B 163 ILE B 165
#
loop_
_pdbx_struct_sheet_hbond.sheet_id
_pdbx_struct_sheet_hbond.range_id_1
_pdbx_struct_sheet_hbond.range_id_2
_pdbx_struct_sheet_hbond.range_1_label_atom_id
_pdbx_struct_sheet_hbond.range_1_label_comp_id
_pdbx_struct_sheet_hbond.range_1_label_asym_id
_pdbx_struct_sheet_hbond.range_1_label_seq_id
_pdbx_struct_sheet_hbond.range_1_PDB_ins_code
_pdbx_struct_sheet_hbond.range_1_auth_atom_id
_pdbx_struct_sheet_hbond.range_1_auth_comp_id
_pdbx_struct_sheet_hbond.range_1_auth_asym_id
_pdbx_struct_sheet_hbond.range_1_auth_seq_id
_pdbx_struct_sheet_hbond.range_2_label_atom_id
_pdbx_struct_sheet_hbond.range_2_label_comp_id
_pdbx_struct_sheet_hbond.range_2_label_asym_id
_pdbx_struct_sheet_hbond.range_2_label_seq_id
_pdbx_struct_sheet_hbond.range_2_PDB_ins_code
_pdbx_struct_sheet_hbond.range_2_auth_atom_id
_pdbx_struct_sheet_hbond.range_2_auth_comp_id
_pdbx_struct_sheet_hbond.range_2_auth_asym_id
_pdbx_struct_sheet_hbond.range_2_auth_seq_id
AA1 1 2 N ILE A 5   ? N ILE A 7   O LEU A 93  ? O LEU A 95
AA1 2 3 O ILE A 92  ? O ILE A 94  N LEU A 74  ? N LEU A 76
AA1 3 4 O VAL A 79  ? O VAL A 81  N PHE A 33  ? N PHE A 35
AA2 1 2 N VAL A 20  ? N VAL A 22  O LEU A 58  ? O LEU A 60
AA2 2 3 O LYS A 59  ? O LYS A 61  N ILE A 50  ? N ILE A 52
AA3 1 2 N GLU A 105 ? N GLU A 107 O THR A 131 ? O THR A 133
AA4 1 2 N VAL A 115 ? N VAL A 117 O THR A 210 ? O THR A 212
AA4 2 3 O ALA A 205 ? O ALA A 207 N LEU A 187 ? N LEU A 189
AA4 3 4 O GLN A 190 ? O GLN A 192 N THR A 147 ? N THR A 149
AA5 1 2 N MET A 126 ? N MET A 128 O ILE A 170 ? O ILE A 172
AA5 2 3 O SER A 171 ? O SER A 173 N THR A 162 ? N THR A 164
AA6 1 2 N ILE B 5   ? N ILE B 7   O THR B 95  ? O THR B 97
AA6 2 3 O ILE B 92  ? O ILE B 94  N LEU B 74  ? N LEU B 76
AA6 3 4 O HIS B 77  ? O HIS B 79  N SER B 35  ? N SER B 37
AA7 1 2 N VAL B 20  ? N VAL B 22  O LEU B 58  ? O LEU B 60
AA7 2 3 O LYS B 59  ? O LYS B 61  N ILE B 50  ? N ILE B 52
AA8 1 2 N GLU B 105 ? N GLU B 107 O THR B 131 ? O THR B 133
AA9 1 2 N PHE B 111 ? N PHE B 113 O THR B 204 ? O THR B 206
AA9 2 3 O ALA B 205 ? O ALA B 207 N LEU B 187 ? N LEU B 189
AA9 3 4 O THR B 186 ? O THR B 188 N ASP B 152 ? N ASP B 154
AB1 1 2 N MET B 126 ? N MET B 128 O ILE B 170 ? O ILE B 172
AB1 2 3 O SER B 171 ? O SER B 173 N THR B 162 ? N THR B 164
# """

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
A 1 SER A 4  ? ALA A 7  ? ? SER A 4  ALA A 7
A 2 MET A 22 ? GLU A 25 ? ? MET A 22 GLU A 25
#
_pdbx_struct_sheet_hbond.sheet_id                A
_pdbx_struct_sheet_hbond.range_id_1              1
_pdbx_struct_sheet_hbond.range_id_2              2
_pdbx_struct_sheet_hbond.range_1_label_atom_id   N
_pdbx_struct_sheet_hbond.range_1_label_comp_id   SER
_pdbx_struct_sheet_hbond.range_1_label_asym_id   A
_pdbx_struct_sheet_hbond.range_1_label_seq_id    4
_pdbx_struct_sheet_hbond.range_1_PDB_ins_code    ?
_pdbx_struct_sheet_hbond.range_1_auth_atom_id    N
_pdbx_struct_sheet_hbond.range_1_auth_comp_id    SER
_pdbx_struct_sheet_hbond.range_1_auth_asym_id    A
_pdbx_struct_sheet_hbond.range_1_auth_seq_id     4
_pdbx_struct_sheet_hbond.range_2_label_atom_id   O
_pdbx_struct_sheet_hbond.range_2_label_comp_id   GLU
_pdbx_struct_sheet_hbond.range_2_label_asym_id   A
_pdbx_struct_sheet_hbond.range_2_label_seq_id    25
_pdbx_struct_sheet_hbond.range_2_PDB_ins_code    ?
_pdbx_struct_sheet_hbond.range_2_auth_atom_id    O
_pdbx_struct_sheet_hbond.range_2_auth_comp_id    GLU
_pdbx_struct_sheet_hbond.range_2_auth_asym_id    A
_pdbx_struct_sheet_hbond.range_2_auth_seq_id     25
#
"""

def tst_from_cif_block():
  test_str_1 = "data_4ZTE\n" + cif_4_helices + cif_many_sheets
  pdb_str = """\
HELIX    1 AA1 SER A   26  LYS A   30  5                                   5
HELIX    2 AA2 LEU A  196  GLU A  199  5                                   4
HELIX    3 AA3 SER B   26  LYS B   30  5                                   5
HELIX    4 AA4 LEU B  196  GLU B  199  5                                   4
SHEET    1 AA1 4 ILE A   7  PRO A  10  0
SHEET    2 AA1 4 MET A  92  THR A  99  1  O  LEU A  95   N  ILE A   7
SHEET    3 AA1 4 THR A  73  SER A  82 -1  N  LEU A  76   O  ILE A  94
SHEET    4 AA1 4 VAL A  34  THR A  39 -1  N  PHE A  35   O  VAL A  81
SHEET    1 AA2 3 LYS A  19  GLN A  23  0
SHEET    2 AA2 3 TRP A  59  VAL A  62 -1  O  LEU A  60   N  VAL A  22
SHEET    3 AA2 3 PHE A  51  ILE A  53 -1  N  ILE A  52   O  LYS A  61
SHEET    1 AA3 2 GLU A 107  PHE A 108  0
SHEET    2 AA3 2 ALA A 132  THR A 133 -1  O  THR A 133   N  GLU A 107
SHEET    1 AA4 4 VAL A 112  MET A 118  0
SHEET    2 AA4 4 SER A 202  THR A 212  1  O  THR A 212   N  VAL A 117
SHEET    3 AA4 4 THR A 186  ALA A 194 -1  N  LEU A 189   O  ALA A 207
SHEET    4 AA4 4 ALA A 147  ASP A 154 -1  N  THR A 149   O  GLN A 192
SHEET    1 AA5 3 THR A 125  GLU A 129  0
SHEET    2 AA5 3 VAL A 171  VAL A 174 -1  O  ILE A 172   N  MET A 128
SHEET    3 AA5 3 PHE A 163  ILE A 165 -1  N  THR A 164   O  SER A 173
SHEET    1 AA6 4 ILE B   7  PRO B  10  0
SHEET    2 AA6 4 MET B  92  THR B  99  1  O  THR B  97   N  ILE B   7
SHEET    3 AA6 4 THR B  73  SER B  82 -1  N  LEU B  76   O  ILE B  94
SHEET    4 AA6 4 VAL B  34  GLY B  40 -1  N  SER B  37   O  HIS B  79
SHEET    1 AA7 3 LYS B  19  GLN B  23  0
SHEET    2 AA7 3 TRP B  59  VAL B  62 -1  O  LEU B  60   N  VAL B  22
SHEET    3 AA7 3 PHE B  51  ILE B  53 -1  N  ILE B  52   O  LYS B  61
SHEET    1 AA8 2 GLU B 107  PHE B 108  0
SHEET    2 AA8 2 ALA B 132  THR B 133 -1  O  THR B 133   N  GLU B 107
SHEET    1 AA9 4 VAL B 112  MET B 118  0
SHEET    2 AA9 4 SER B 202  THR B 212  1  O  THR B 206   N  PHE B 113
SHEET    3 AA9 4 THR B 186  ALA B 194 -1  N  LEU B 189   O  ALA B 207
SHEET    4 AA9 4 ALA B 147  ASP B 154 -1  N  ASP B 154   O  THR B 188
SHEET    1 AB1 3 SER B 126  GLU B 129  0
SHEET    2 AB1 3 VAL B 171  VAL B 174 -1  O  ILE B 172   N  MET B 128
SHEET    3 AB1 3 PHE B 163  ILE B 165 -1  N  THR B 164   O  SER B 173"""

  cif_model = iotbx.cif.reader(input_string=test_str_1).model()
  cif_block = list(cif_model.values())[0]
  ann = annotation.from_cif_block(cif_block)
  assert ann.get_n_helices() == 4
  resnames = [x.start_resname for x in ann.helices]
  assert resnames == ["SER","LEU","SER","LEU"]
  resnames = [x.end_resname for x in ann.helices]
  assert resnames == ["LYS","GLU","LYS","GLU"]
  assert ann.get_n_sheets() == 10, ann.get_n_sheets()
  assert [len(x.strands) for x in ann.sheets] == [4, 3, 2, 4, 3, 4, 3, 2, 4, 3]
  # print ann.as_pdb_str()
  assert not show_diff(pdb_str, ann.as_pdb_str())

def tst_from_cif_block_2():
  cif_str = "data_1UCS\n" + cif_4_helices + cif_one_sheet
  cif_model = iotbx.cif.reader(input_string=cif_str).model()
  cif_block = list(cif_model.values())[0]
  ann = annotation.from_cif_block(cif_block)
  assert ann.get_n_helices() == 4
  resnames = [x.start_resname for x in ann.helices]
  assert resnames == ["SER","LEU","SER","LEU"]
  resnames = [x.end_resname for x in ann.helices]
  assert resnames == ["LYS","GLU","LYS","GLU"]
  assert ann.get_n_sheets() == 1
  assert [len(x.strands) for x in ann.sheets] == [2]
  # print ann.as_pdb_str()

def tst_from_cif_block_3():
  "Only one helix"
  cif_str = "data_1UCS" + cif_one_helix

  cif_model = iotbx.cif.reader(input_string=cif_str).model()
  cif_block = list(cif_model.values())[0]
  ann = annotation.from_cif_block(cif_block)
  assert ann.get_n_helices() == 1
  resnames = [x.start_resname for x in ann.helices]
  assert resnames == ["GLN"], resnames
  resnames = [x.end_resname for x in ann.helices]
  assert resnames == ["GLY"]
  assert ann.get_n_sheets() == 0
  # print ann.as_pdb_str()

def tst_from_cif_block_4():
  "Only one sheet"
  cif_str = "data_1UCS" + cif_one_sheet
  cif_model = iotbx.cif.reader(input_string=cif_str).model()
  cif_block = list(cif_model.values())[0]
  ann = annotation.from_cif_block(cif_block)
  assert ann.get_n_helices() == 0, ann.get_n_helices()
  assert ann.get_n_sheets() == 1
  assert [len(x.strands) for x in ann.sheets] == [2]
  # print ann.as_pdb_str()

def tst_from_minimal_cif_helix():
  """reading annotation that contains only required mmCIF records -
  http://mmcif.wwpdb.org/dictionaries/mmcif_pdbx_v40.dic/Categories/struct_conf.html
  """
  cif_answer_minimal_helix = """\
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
  HELX_P  1  1  ARG  A   87  ?  GLN  A   92  ?  1  ?  5
  HELX_P  2  2  ARG  B  287  ?  GLN  B  292  ?  1  ?  5
"""
  cif_answer_conf_type_loop = """\
loop_
  _struct_conf_type.id
  _struct_conf_type.criteria
  _struct_conf_type.reference
  HELX_P  ?  ?
"""
  pdb_answer_minimal_helix = """\
HELIX    1   1 ARG A   87  GLN A   92  1                                   5
HELIX    2   2 ARG B  287  GLN B  292  1                                   5"""
  minimal_helix = """\
    data_4ZTE
    loop_
    _struct_conf.id
    _struct_conf.conf_type_id
    _struct_conf.beg_label_comp_id
    _struct_conf.beg_label_asym_id
    _struct_conf.beg_label_seq_id
    _struct_conf.end_label_comp_id
    _struct_conf.end_label_asym_id
    _struct_conf.end_label_seq_id
      HELX1  HELX_RH_AL_P  ARG  A   87  GLN  A   92
      HELX2  HELX_RH_AL_P  ARG  B  287  GLN  B  292
      STRN1  STRN_P        PRO  A    1  LEU  A    5
  """
  cif_model = iotbx.cif.reader(input_string=minimal_helix).model()
  cif_block = list(cif_model.values())[0]
  ann = annotation.from_cif_block(cif_block)
  assert not show_diff(ann.as_pdb_str(), pdb_answer_minimal_helix)

  cif_loops = ann.as_cif_loops()
  assert len(cif_loops) == 2
  for i, answer in enumerate([cif_answer_minimal_helix, cif_answer_conf_type_loop]):
    out = StringIO()
    cif_loops[i].show(out)
    v = out.getvalue()
    # print ("\"%s\"" % v)
    assert not show_diff(out.getvalue(), answer)

def tst_from_minimal_cif_sheet():
  cif_minimal_sheet = """\
      data_4ZTE
      #
      _struct_sheet.id               A
      #
      _struct_sheet_order.sheet_id     A
      _struct_sheet_order.range_id_1   1
      _struct_sheet_order.range_id_2   2
      #
      loop_
      _struct_sheet_range.sheet_id
      _struct_sheet_range.id
      _struct_sheet_range.beg_label_comp_id
      _struct_sheet_range.beg_label_asym_id
      _struct_sheet_range.beg_label_seq_id
      _struct_sheet_range.end_label_comp_id
      _struct_sheet_range.end_label_asym_id
      _struct_sheet_range.end_label_seq_id
      A 1 SER A 4  ALA A 7
      A 2 MET A 22 GLU A 25
      #
      _pdbx_struct_sheet_hbond.sheet_id                A
      _pdbx_struct_sheet_hbond.range_id_1              1
      _pdbx_struct_sheet_hbond.range_id_2              2
      _pdbx_struct_sheet_hbond.range_1_label_atom_id   N
      _pdbx_struct_sheet_hbond.range_1_label_seq_id    4
      _pdbx_struct_sheet_hbond.range_2_label_atom_id   O
      _pdbx_struct_sheet_hbond.range_2_label_seq_id    25
      #
      """
  pdb_answer_minimal_sheet = """\
SHEET    1   A 2 SER A   4  ALA A   7  0
SHEET    2   A 2 MET A  22  GLU A  25  0  O    . .  25   N    . .   4"""
  answer_struct_sheet_loops = ["""\
loop_
  _struct_sheet.id
  _struct_sheet.type
  _struct_sheet.number_strands
  _struct_sheet.details
  A  ?  2  ?
""","""\
loop_
  _struct_sheet_order.sheet_id
  _struct_sheet_order.range_id_1
  _struct_sheet_order.range_id_2
  _struct_sheet_order.offset
  _struct_sheet_order.sense
  A  1  2  ?  ?
""","""\
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
  A  1  SER  A   4  ?  ALA  A   7  ?
  A  2  MET  A  22  ?  GLU  A  25  ?
""","""\
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
  A  1  2  N  .  .  4  ?  O  .  .  25  ?
"""]

  cif_model = iotbx.cif.reader(input_string=cif_minimal_sheet).model()
  cif_block = list(cif_model.values())[0]
  ann = annotation.from_cif_block(cif_block)
  # print ann.as_pdb_str()
  assert not show_diff(ann.as_pdb_str(), pdb_answer_minimal_sheet)

  cif_loops = ann.as_cif_loops()
  assert len(cif_loops) == 4, len(cif_loops)
  for i, sheet_loop in enumerate(cif_loops):
    out = StringIO()
    sheet_loop.show(out)
    v = out.getvalue()
    # print "\"%s\"" % v
    assert not show_diff(out.getvalue(), answer_struct_sheet_loops[i])

def tst_to_cif_helix():
  pdb_str = """\
HELIX    1 AA1 SER A   26  LYS A   30  5                                   5
HELIX    2 AA2 LEU A  196  GLU A  199  5                                   4
HELIX    3 AA3 SER B   26  LYS B   30  5                                   5
HELIX    4 AA4 LEU B  196  GLU B  199  5                                   4
"""
  answer = """\
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
  HELX_P  1  AA1  SER  A   26  ?  LYS  A   30  ?  5  ?  5
  HELX_P  2  AA2  LEU  A  196  ?  GLU  A  199  ?  5  ?  4
  HELX_P  3  AA3  SER  B   26  ?  LYS  B   30  ?  5  ?  5
  HELX_P  4  AA4  LEU  B  196  ?  GLU  B  199  ?  5  ?  4
"""

  pdb_str2 = """\
HELIX    1 AA1 SER A   26  LYS A   30  5                                   5
"""
  answer2 = """\
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
  HELX_P  1  AA1  SER  A  26  ?  LYS  A  30  ?  5  ?  5\n"""

  ann = annotation.from_records(pdb_str.split("\n"))
  cif_loops = ann.as_cif_loops()
  assert len(cif_loops) == 2
  helix_loop = cif_loops[0]
  out = StringIO()
  helix_loop.show(out)
  v = out.getvalue()
  # print("\"%s\"" % v)
  assert not show_diff(out.getvalue(), answer)

  # hmmm... when there's only one chain, there is one less 'space'
  # between resname and resnum. Guess: in previous case extra space is
  # needed to make nice columns.
  ann = annotation.from_records(pdb_str2.split("\n"))
  cif_loops = ann.as_cif_loops()
  assert len(cif_loops) == 2
  helix_loop = cif_loops[0]
  out = StringIO()
  helix_loop.show(out)
  v = out.getvalue()
  # print "\"%s\"" % v
  assert not show_diff(out.getvalue(), answer2)


def tst_to_cif_sheet():
  pdb_str1 = """\
SHEET    1 AA1 4 ILE A   7  PRO A  10  0
SHEET    2 AA1 4 MET A  92  THR A  99  1  O  LEU A  95   N  ILE A   7
SHEET    3 AA1 4 THR A  73  SER A  82 -1  N  LEU A  76   O  ILE A  94
SHEET    4 AA1 4 VAL A  34  THR A  39 -1  N  PHE A  35   O  VAL A  81
SHEET    1 AA2 3 LYS A  19  GLN A  23  0
SHEET    2 AA2 3 TRP A  59  VAL A  62 -1  O  LEU A  60   N  VAL A  22
SHEET    3 AA2 3 PHE A  51  ILE A  53 -1  N  ILE A  52   O  LYS A  61"""

  ann = annotation.from_records(pdb_str1.split("\n"))
  cif_loops = ann.as_cif_loops()
  assert len(cif_loops) == 4, len(cif_loops)

  out = StringIO()
  cif_block = iotbx.cif.model.block()
  for loop in cif_loops:
    loop.show(out)
    cif_block.add_loop(loop)
  v = out.getvalue()
  # print "\"%s\"" % v
  ann_back = annotation.from_cif_block(cif_block)
  assert ann_back.get_n_helices() == 0
  assert ann_back.get_n_sheets() == 2
  assert [len(x.strands) for x in ann_back.sheets] == [4, 3]
  # print ann_back.as_pdb_str()
  assert not show_diff(ann_back.as_pdb_str(), pdb_str1)

  pdb_str2 = """\
SHEET    1 AA1 4 ILE A   7  PRO A  10  0
SHEET    2 AA1 4 MET A  92  THR A  99  1  O  LEU A  95   N  ILE A   7
SHEET    3 AA1 4 THR A  73  SER A  82 -1  N  LEU A  76   O  ILE A  94
SHEET    4 AA1 4 VAL A  34  THR A  39 -1  N  PHE A  35   O  VAL A  81"""

  ann = annotation.from_records(pdb_str2.split("\n"))
  cif_loops = ann.as_cif_loops()
  assert len(cif_loops) == 4, len(cif_loops)

  out = StringIO()
  cif_block = iotbx.cif.model.block()
  for loop in cif_loops:
    loop.show(out)
    cif_block.add_loop(loop)
  v = out.getvalue()
  # print "\"%s\"" % v
  ann_back = annotation.from_cif_block(cif_block)
  assert ann_back.get_n_helices() == 0
  assert ann_back.get_n_sheets() == 1
  assert [len(x.strands) for x in ann_back.sheets] == [4]
  # print ann_back.as_pdb_str()
  assert not show_diff(ann_back.as_pdb_str(), pdb_str2)


pdb_records_2 = """\
HELIX   10  10 MET B    3  TYR B   15  1                                  13
HELIX   11  11 TYR B  220  LEU B  241  1                                  22
HELIX   12  12 GLU B  247  LEU B  275  1                                  29
HELIX   13  13 ILE B  281  ARG B  307  1                                  27
HELIX   14  14 GLU B  403  TYR B  430  1                                  28
HELIX   15  15 ALA B  432  VAL B  462  1                                  31
HELIX   16  16 PRO B  463  ASN B  466  5                                   4
SHEET    1  BA 6 ASP B  77  PRO B  81  0
SHEET    2  BA 6 ASN B 107  GLN B 111 -1  O  VAL B 108   N  ILE B  80
SHEET    3  BA 6 ALA B 115  TRP B 118 -1  O  ALA B 115   N  GLN B 111
SHEET    4  BA 6 ASN B  57  THR B  61 -1  O  LEU B  58   N  TRP B 118
SHEET    5  BA 6 THR B  30  GLY B  34 -1  O  ARG B  32   N  ALA B  59
SHEET    6  BA 6 ILE B 157  GLN B 159  1  O  ILE B 157   N  VAL B  31
SHEET    1  BB 4 LEU B  37  ASN B  44  0
SHEET    2  BB 4 GLU B  49  PHE B  55 -1  O  GLU B  49   N  ASN B  44
SHEET    3  BB 4 SER B 121  SER B 126 -1  O  ALA B 122   N  VAL B  54
SHEET    4  BB 4 PHE B 100  GLU B 101 -1  O  GLU B 101   N  ILE B 123
SHEET    1  BC 4 ILE B  90  LEU B  92  0
SHEET    2  BC 4 THR B 143  SER B 148 -1  O  LYS B 147   N  VAL B  91
SHEET    3  BC 4 VAL B 207  TYR B 210 -1  O  VAL B 207   N  PHE B 146
SHEET    4  BC 4 SER B 193  LYS B 195 -1  O  ARG B 194   N  THR B 208
SHEET    1  BD 2 SER B 187  ILE B 188  0
SHEET    2  BD 2 ILE B 213  GLN B 214 -1  O  GLN B 214   N  SER B 187
"""
pdb_string_2 = """\
ATOM   4705  N   TYR B 220     236.756 243.825 262.336  1.00 10.00           N
ATOM   4706  CA  TYR B 220     236.027 243.470 261.100  1.00 10.00           C
ATOM   4707  C   TYR B 220     235.677 244.591 260.121  1.00 10.00           C
ATOM   4708  O   TYR B 220     235.299 244.306 258.994  1.00 10.00           O
ATOM   4709  CB  TYR B 220     234.714 242.742 261.381  1.00 10.00           C
ATOM   4710  CG  TYR B 220     234.735 241.529 262.275  1.00 10.00           C
ATOM   4711  CD1 TYR B 220     235.783 241.294 263.193  1.00 10.00           C
ATOM   4712  CD2 TYR B 220     233.619 240.688 262.302  1.00 10.00           C
ATOM   4713  CE1 TYR B 220     235.700 240.257 264.107  1.00 10.00           C
ATOM   4714  CE2 TYR B 220     233.517 239.671 263.204  1.00 10.00           C
ATOM   4715  CZ  TYR B 220     234.545 239.451 264.103  1.00 10.00           C
ATOM   4716  OH  TYR B 220     234.389 238.416 264.996  1.00 10.00           O
ATOM   4717  N   ILE B 221     235.740 245.846 260.556  1.00 10.00           N
ATOM   4718  CA  ILE B 221     235.458 246.970 259.704  1.00 10.00           C
ATOM   4719  C   ILE B 221     236.494 246.879 258.537  1.00 10.00           C
ATOM   4720  O   ILE B 221     236.314 247.422 257.427  1.00 10.00           O
ATOM   4721  CB  ILE B 221     235.531 248.227 260.635  1.00 10.00           C
ATOM   4722  CG1 ILE B 221     236.909 248.327 261.313  1.00 10.00           C
ATOM   4723  CG2 ILE B 221     234.481 248.053 261.770  1.00 10.00           C
ATOM   4724  CD1 ILE B 221     236.886 248.507 262.830  1.00 10.00           C
ATOM   4725  N   VAL B 222     237.527 246.072 258.779  1.00 10.00           N
ATOM   4726  CA  VAL B 222     238.611 245.866 257.840  1.00 10.00           C
ATOM   4727  C   VAL B 222     238.076 245.079 256.678  1.00 10.00           C
ATOM   4728  O   VAL B 222     238.190 245.501 255.525  1.00 10.00           O
ATOM   4729  CB  VAL B 222     239.767 245.051 258.480  1.00 10.00           C
ATOM   4730  CG1 VAL B 222     241.030 245.163 257.617  1.00 10.00           C
ATOM   4731  CG2 VAL B 222     239.993 245.505 259.937  1.00 10.00           C
ATOM   4732  N   TYR B 223     237.498 243.916 256.959  1.00 10.00           N
ATOM   4733  CA  TYR B 223     236.998 243.193 255.797  1.00 10.00           C
ATOM   4734  C   TYR B 223     235.642 243.731 255.357  1.00 10.00           C
ATOM   4735  O   TYR B 223     235.112 243.268 254.374  1.00 10.00           O
ATOM   4736  CB  TYR B 223     236.908 241.695 256.094  1.00 10.00           C
ATOM   4737  CG  TYR B 223     235.804 241.324 257.057  1.00 10.00           C
ATOM   4738  CD1 TYR B 223     234.480 241.272 256.641  1.00 10.00           C
ATOM   4739  CD2 TYR B 223     236.085 241.027 258.384  1.00 10.00           C
ATOM   4740  CE1 TYR B 223     233.468 240.932 257.518  1.00 10.00           C
ATOM   4741  CE2 TYR B 223     235.080 240.686 259.269  1.00 10.00           C
ATOM   4742  CZ  TYR B 223     233.774 240.641 258.831  1.00 10.00           C
ATOM   4743  OH  TYR B 223     232.772 240.302 259.711  1.00 10.00           O
ATOM   4744  N   THR B 224     235.127 244.742 256.068  1.00 10.00           N
ATOM   4745  CA  THR B 224     233.888 245.465 255.737  1.00 10.00           C
ATOM   4746  C   THR B 224     234.274 246.452 254.592  1.00 10.00           C
ATOM   4747  O   THR B 224     233.527 246.682 253.621  1.00 10.00           O
ATOM   4748  CB  THR B 224     233.390 246.267 256.985  1.00 10.00           C
ATOM   4749  OG1 THR B 224     232.815 245.365 257.938  1.00 10.00           O
ATOM   4750  CG2 THR B 224     232.338 247.298 256.604  1.00 10.00           C
ATOM   4751  N   ILE B 225     235.481 246.993 254.709  1.00 10.00           N
ATOM   4752  CA  ILE B 225     236.079 247.939 253.743  1.00 10.00           C
ATOM   4753  C   ILE B 225     236.203 247.150 252.435  1.00 10.00           C
ATOM   4754  O   ILE B 225     235.800 247.601 251.358  1.00 10.00           O
ATOM   4755  CB  ILE B 225     237.546 248.405 254.230  1.00 10.00           C
ATOM   4756  CG1 ILE B 225     237.487 249.823 254.873  1.00 10.00           C
ATOM   4757  CG2 ILE B 225     238.554 248.277 253.078  1.00 10.00           C
ATOM   4758  CD1 ILE B 225     238.802 250.372 255.483  1.00 10.00           C
ATOM   4759  N   VAL B 226     236.723 245.926 252.549  1.00 10.00           N
ATOM   4760  CA  VAL B 226     236.904 245.129 251.334  1.00 10.00           C
ATOM   4761  C   VAL B 226     235.606 245.025 250.571  1.00 10.00           C
ATOM   4762  O   VAL B 226     235.621 245.119 249.364  1.00 10.00           O
ATOM   4763  CB  VAL B 226     237.478 243.708 251.567  1.00 10.00           C
ATOM   4764  CG1 VAL B 226     238.787 243.804 252.297  1.00 10.00           C
ATOM   4765  CG2 VAL B 226     236.480 242.836 252.257  1.00 10.00           C
ATOM   4766  N   PRO B 227     234.471 244.784 251.253  1.00 10.00           N
ATOM   4767  CA  PRO B 227     233.225 244.717 250.484  1.00 10.00           C
ATOM   4768  C   PRO B 227     233.110 246.013 249.736  1.00 10.00           C
ATOM   4769  O   PRO B 227     233.251 246.069 248.522  1.00 10.00           O
ATOM   4770  CB  PRO B 227     232.142 244.666 251.550  1.00 10.00           C
ATOM   4771  CG  PRO B 227     232.785 243.987 252.657  1.00 10.00           C
ATOM   4772  CD  PRO B 227     234.206 244.487 252.661  1.00 10.00           C
ATOM   4773  N   CYS B 228     232.868 247.074 250.497  1.00 10.00           N
ATOM   4774  CA  CYS B 228     232.716 248.396 249.898  1.00 10.00           C
ATOM   4775  C   CYS B 228     233.629 248.546 248.681  1.00 10.00           C
ATOM   4776  O   CYS B 228     233.127 248.848 247.618  1.00 10.00           O
ATOM   4777  CB  CYS B 228     232.947 249.541 250.931  1.00 10.00           C
ATOM   4778  SG  CYS B 228     231.542 249.841 252.106  1.00 10.00           S
ATOM   4779  N   ILE B 229     234.929 248.269 248.819  1.00 10.00           N
ATOM   4780  CA  ILE B 229     235.889 248.399 247.719  1.00 10.00           C
ATOM   4781  C   ILE B 229     235.932 247.224 246.730  1.00 10.00           C
ATOM   4782  O   ILE B 229     236.051 247.422 245.518  1.00 10.00           O
ATOM   4783  CB  ILE B 229     237.313 248.607 248.276  1.00 10.00           C
ATOM   4784  CG1 ILE B 229     237.316 249.781 249.260  1.00 10.00           C
ATOM   4785  CG2 ILE B 229     238.288 248.863 247.136  1.00 10.00           C
ATOM   4786  CD1 ILE B 229     238.648 250.007 249.947  1.00 10.00           C
ATOM   4787  N   LEU B 230     235.860 246.004 247.250  1.00 10.00           N
ATOM   4788  CA  LEU B 230     235.889 244.785 246.426  1.00 10.00           C
ATOM   4789  C   LEU B 230     234.553 244.686 245.700  1.00 10.00           C
ATOM   4790  O   LEU B 230     234.483 244.299 244.538  1.00 10.00           O
ATOM   4791  CB  LEU B 230     236.102 243.542 247.297  1.00 10.00           C
ATOM   4792  CG  LEU B 230     237.259 243.504 248.306  1.00 10.00           C
ATOM   4793  CD1 LEU B 230     237.222 242.152 248.975  1.00 10.00           C
ATOM   4794  CD2 LEU B 230     238.606 243.767 247.657  1.00 10.00           C
ATOM   4795  N   ILE B 231     233.483 245.047 246.402  1.00 10.00           N
ATOM   4796  CA  ILE B 231     232.175 245.042 245.773  1.00 10.00           C
ATOM   4797  C   ILE B 231     232.270 246.087 244.707  1.00 10.00           C
ATOM   4798  O   ILE B 231     231.634 245.986 243.664  1.00 10.00           O
ATOM   4799  CB  ILE B 231     231.097 245.478 246.713  1.00 10.00           C
ATOM   4800  CG1 ILE B 231     231.185 244.615 247.972  1.00 10.00           C
ATOM   4801  CG2 ILE B 231     229.721 245.466 245.969  1.00 10.00           C
ATOM   4802  CD1 ILE B 231     230.710 243.195 247.766  1.00 10.00           C
ATOM   4803  N   SER B 232     233.034 247.118 245.052  1.00 10.00           N
ATOM   4804  CA  SER B 232     233.308 248.151 244.061  1.00 10.00           C
ATOM   4805  C   SER B 232     234.054 247.577 242.862  1.00 10.00           C
ATOM   4806  O   SER B 232     233.835 247.929 241.753  1.00 10.00           O
ATOM   4807  CB  SER B 232     234.101 249.298 244.687  1.00 10.00           C
ATOM   4808  OG  SER B 232     233.261 250.136 245.462  1.00 10.00           O
ATOM   4809  N   ILE B 233     234.905 246.580 243.194  1.00 10.00           N
ATOM   4810  CA  ILE B 233     235.644 245.926 242.139  1.00 10.00           C
ATOM   4811  C   ILE B 233     234.665 245.443 241.097  1.00 10.00           C
ATOM   4812  O   ILE B 233     234.778 245.798 239.927  1.00 10.00           O
ATOM   4813  CB  ILE B 233     236.469 244.734 242.640  1.00 10.00           C
ATOM   4814  CG1 ILE B 233     237.195 245.075 243.948  1.00 10.00           C
ATOM   4815  CG2 ILE B 233     237.521 244.395 241.583  1.00 10.00           C
ATOM   4816  CD1 ILE B 233     238.246 244.065 244.374  1.00 10.00           C
ATOM   4817  N   LEU B 234     233.675 244.677 241.533  1.00 10.00           N
ATOM   4818  CA  LEU B 234     232.673 244.134 240.627  1.00 10.00           C
ATOM   4819  C   LEU B 234     231.791 245.340 240.186  1.00 10.00           C
ATOM   4820  O   LEU B 234     231.322 245.402 239.046  1.00 10.00           O
ATOM   4821  CB  LEU B 234     231.877 242.962 241.316  1.00 10.00           C
ATOM   4822  CG  LEU B 234     232.447 241.937 242.366  1.00 10.00           C
ATOM   4823  CD1 LEU B 234     231.422 240.857 242.578  1.00 10.00           C
ATOM   4824  CD2 LEU B 234     233.763 241.285 241.947  1.00 10.00           C
ATOM   4825  N   ALA B 235     231.631 246.334 241.067  1.00 10.00           N
ATOM   4826  CA  ALA B 235     230.861 247.525 240.769  1.00 10.00           C
ATOM   4827  C   ALA B 235     231.553 248.274 239.623  1.00 10.00           C
ATOM   4828  O   ALA B 235     231.035 248.323 238.503  1.00 10.00           O
ATOM   4829  CB  ALA B 235     230.785 248.393 242.002  1.00 10.00           C
ATOM   4830  N   ILE B 236     232.691 248.903 239.910  1.00 10.00           N
ATOM   4831  CA  ILE B 236     233.475 249.607 238.866  1.00 10.00           C
ATOM   4832  C   ILE B 236     233.551 248.710 237.664  1.00 10.00           C
ATOM   4833  O   ILE B 236     233.594 249.167 236.534  1.00 10.00           O
ATOM   4834  CB  ILE B 236     234.979 249.942 239.255  1.00 10.00           C
ATOM   4835  CG1 ILE B 236     235.762 250.300 237.985  1.00 10.00           C
ATOM   4836  CG2 ILE B 236     235.693 248.763 239.889  1.00 10.00           C
ATOM   4837  CD1 ILE B 236     236.838 251.361 238.197  1.00 10.00           C
ATOM   4838  N   LEU B 237     233.558 247.415 237.957  1.00 10.00           N
ATOM   4839  CA  LEU B 237     233.632 246.343 236.959  1.00 10.00           C
ATOM   4840  C   LEU B 237     232.360 246.328 236.119  1.00 10.00           C
ATOM   4841  O   LEU B 237     232.416 246.266 234.898  1.00 10.00           O
ATOM   4842  CB  LEU B 237     233.818 244.973 237.656  1.00 10.00           C
ATOM   4843  CG  LEU B 237     235.119 244.111 237.601  1.00 10.00           C
ATOM   4844  CD1 LEU B 237     236.403 244.933 237.384  1.00 10.00           C
ATOM   4845  CD2 LEU B 237     235.234 243.311 238.886  1.00 10.00           C
ATOM   4846  N   VAL B 238     231.212 246.422 236.778  1.00 10.00           N
ATOM   4847  CA  VAL B 238     229.915 246.397 236.074  1.00 10.00           C
ATOM   4848  C   VAL B 238     229.713 247.710 235.345  1.00 10.00           C
ATOM   4849  O   VAL B 238     228.903 247.810 234.436  1.00 10.00           O
ATOM   4850  CB  VAL B 238     228.702 246.194 237.050  1.00 10.00           C
ATOM   4851  CG1 VAL B 238     227.400 246.551 236.389  1.00 10.00           C
ATOM   4852  CG2 VAL B 238     228.631 244.739 237.510  1.00 10.00           C
ATOM   4853  N   PHE B 239     230.476 248.722 235.747  1.00 10.00           N
ATOM   4854  CA  PHE B 239     230.405 250.061 235.171  1.00 10.00           C
ATOM   4855  C   PHE B 239     231.509 250.132 234.127  1.00 10.00           C
ATOM   4856  O   PHE B 239     231.467 250.917 233.179  1.00 10.00           O
ATOM   4857  CB  PHE B 239     230.634 251.028 236.310  1.00 10.00           C
ATOM   4858  CG  PHE B 239     229.961 250.592 237.595  1.00 10.00           C
ATOM   4859  CD1 PHE B 239     228.755 249.886 237.561  1.00 10.00           C
ATOM   4860  CD2 PHE B 239     230.531 250.863 238.836  1.00 10.00           C
ATOM   4861  CE1 PHE B 239     228.132 249.453 238.735  1.00 10.00           C
ATOM   4862  CE2 PHE B 239     229.914 250.434 240.016  1.00 10.00           C
ATOM   4863  CZ  PHE B 239     228.711 249.724 239.965  1.00 10.00           C
ATOM   4864  N   TYR B 240     232.465 249.220 234.293  1.00 10.00           N
ATOM   4865  CA  TYR B 240     233.634 249.123 233.423  1.00 10.00           C
ATOM   4866  C   TYR B 240     233.399 248.048 232.375  1.00 10.00           C
ATOM   4867  O   TYR B 240     233.905 248.138 231.261  1.00 10.00           O
ATOM   4868  CB  TYR B 240     234.847 248.758 234.275  1.00 10.00           C
ATOM   4869  CG  TYR B 240     236.188 249.041 233.699  1.00 10.00           C
ATOM   4870  CD1 TYR B 240     236.368 249.237 232.358  1.00 10.00           C
ATOM   4871  CD2 TYR B 240     237.291 249.069 234.516  1.00 10.00           C
ATOM   4872  CE1 TYR B 240     237.606 249.449 231.854  1.00 10.00           C
ATOM   4873  CE2 TYR B 240     238.527 249.277 234.016  1.00 10.00           C
ATOM   4874  CZ  TYR B 240     238.675 249.463 232.688  1.00 10.00           C
ATOM   4875  OH  TYR B 240     239.923 249.641 232.191  1.00 10.00           O
ATOM   4876  N   LEU B 241     232.624 247.038 232.756  1.00 10.00           N
ATOM   4877  CA  LEU B 241     232.280 245.893 231.900  1.00 10.00           C
ATOM   4878  C   LEU B 241     231.473 246.333 230.662  1.00 10.00           C
ATOM   4879  O   LEU B 241     231.580 245.705 229.611  1.00 10.00           O
ATOM   4880  CB  LEU B 241     231.532 244.785 232.730  1.00 10.00           C
ATOM   4881  CG  LEU B 241     230.971 243.450 232.158  1.00 10.00           C
ATOM   4882  CD1 LEU B 241     232.007 242.657 231.452  1.00 10.00           C
ATOM   4883  CD2 LEU B 241     230.426 242.604 233.258  1.00 10.00           C
ATOM   4884  N   PRO B 242     230.679 247.430 230.757  1.00 10.00           N
ATOM   4885  CA  PRO B 242     229.911 247.858 229.579  1.00 10.00           C
ATOM   4886  C   PRO B 242     230.801 248.243 228.376  1.00 10.00           C
ATOM   4887  O   PRO B 242     230.582 247.785 227.252  1.00 10.00           O
ATOM   4888  CB  PRO B 242     229.064 249.002 230.127  1.00 10.00           C
ATOM   4889  CG  PRO B 242     229.728 249.370 231.472  1.00 10.00           C
ATOM   4890  CD  PRO B 242     230.222 248.099 231.985  1.00 10.00           C
ATOM   4891  N   PRO B 243     231.814 249.099 228.594  1.00 10.00           N
ATOM   4892  CA  PRO B 243     232.577 249.351 227.367  1.00 10.00           C
ATOM   4893  C   PRO B 243     233.753 248.371 227.292  1.00 10.00           C
ATOM   4894  O   PRO B 243     234.662 248.514 226.477  1.00 10.00           O
ATOM   4895  CB  PRO B 243     233.035 250.797 227.536  1.00 10.00           C
ATOM   4896  CG  PRO B 243     233.181 250.938 228.993  1.00 10.00           C
ATOM   4897  CD  PRO B 243     231.933 250.254 229.505  1.00 10.00           C
ATOM   4898  N   ASP B 244     233.711 247.382 228.180  1.00 10.00           N
ATOM   4899  CA  ASP B 244     234.715 246.317 228.295  1.00 10.00           C
ATOM   4900  C   ASP B 244     234.177 245.111 227.532  1.00 10.00           C
ATOM   4901  O   ASP B 244     234.866 244.491 226.708  1.00 10.00           O
ATOM   4902  CB  ASP B 244     234.922 245.955 229.765  1.00 10.00           C
ATOM   4903  CG  ASP B 244     236.381 245.800 230.130  1.00 10.00           C
ATOM   4904  OD1 ASP B 244     237.220 245.632 229.218  1.00 10.00           O
ATOM   4905  OD2 ASP B 244     236.685 245.845 231.339  1.00 10.00           O
ATOM   4906  N   ALA B 245     232.915 244.822 227.773  1.00 10.00           N
ATOM   4907  CA  ALA B 245     232.240 243.722 227.116  1.00 10.00           C
ATOM   4908  C   ALA B 245     230.748 244.083 227.062  1.00 10.00           C
ATOM   4909  O   ALA B 245     229.918 243.550 227.807  1.00 10.00           O
ATOM   4910  CB  ALA B 245     232.453 242.495 227.884  1.00 10.00           C
ATOM   4911  N   GLY B 246     230.420 245.050 226.199  1.00 10.00           N
ATOM   4912  CA  GLY B 246     229.039 245.509 226.032  1.00 10.00           C
ATOM   4913  C   GLY B 246     228.022 244.400 226.236  1.00 10.00           C
ATOM   4914  O   GLY B 246     226.971 244.623 226.816  1.00 10.00           O
ATOM   4915  N   GLU B 247     228.338 243.213 225.704  1.00 10.00           N
ATOM   4916  CA  GLU B 247     227.477 242.051 225.888  1.00 10.00           C
ATOM   4917  C   GLU B 247     227.187 241.806 227.364  1.00 10.00           C
ATOM   4918  O   GLU B 247     226.485 240.913 227.745  1.00 10.00           O
ATOM   4919  CB  GLU B 247     228.111 240.810 225.256  1.00 10.00           C
ATOM   4920  CG  GLU B 247     228.267 240.893 223.745  1.00 10.00           C
ATOM   4921  CD  GLU B 247     229.511 241.652 223.330  1.00 10.00           C
ATOM   4922  OE1 GLU B 247     230.289 242.054 224.223  1.00 10.00           O
ATOM   4923  OE2 GLU B 247     229.711 241.846 222.113  1.00 10.00           O
ATOM   4924  N   LYS B 248     227.836 242.656 228.158  1.00 10.00           N
ATOM   4925  CA  LYS B 248     227.703 242.635 229.601  1.00 10.00           C
ATOM   4926  C   LYS B 248     226.199 242.758 229.827  1.00 10.00           C
ATOM   4927  O   LYS B 248     225.647 242.319 230.840  1.00 10.00           O
ATOM   4928  CB  LYS B 248     228.497 243.831 230.186  1.00 10.00           C
ATOM   4929  CG  LYS B 248     227.758 245.140 230.358  1.00 10.00           C
ATOM   4930  CD  LYS B 248     227.205 245.217 231.753  1.00 10.00           C
ATOM   4931  CE  LYS B 248     228.297 245.156 232.785  1.00 10.00           C
ATOM   4932  NZ  LYS B 248     227.742 245.174 234.142  1.00 10.00           N
ATOM   4933  N   MET B 249     225.551 243.310 228.807  1.00 10.00           N
ATOM   4934  CA  MET B 249     224.123 243.560 228.763  1.00 10.00           C
ATOM   4935  C   MET B 249     223.315 242.787 229.811  1.00 10.00           C
ATOM   4936  O   MET B 249     222.738 243.375 230.729  1.00 10.00           O
ATOM   4937  CB  MET B 249     223.631 243.270 227.335  1.00 10.00           C
ATOM   4938  CG  MET B 249     224.404 242.140 226.624  1.00 10.00           C
ATOM   4939  SD  MET B 249     224.187 240.510 227.326  1.00 10.00           S
ATOM   4940  CE  MET B 249     222.349 240.281 227.080  1.00 10.00           C
ATOM   4941  N   SER B 250     223.305 241.467 229.682  1.00 10.00           N
ATOM   4942  CA  SER B 250     222.578 240.589 230.594  1.00 10.00           C
ATOM   4943  C   SER B 250     223.411 240.238 231.807  1.00 10.00           C
ATOM   4944  O   SER B 250     222.980 240.431 232.947  1.00 10.00           O
ATOM   4945  CB  SER B 250     222.169 239.291 229.891  1.00 10.00           C
ATOM   4946  OG  SER B 250     223.305 238.557 229.461  1.00 10.00           O
ATOM   4947  N   LEU B 251     224.606 239.715 231.551  1.00 10.00           N
ATOM   4948  CA  LEU B 251     225.518 239.325 232.623  1.00 10.00           C
ATOM   4949  C   LEU B 251     225.545 240.419 233.656  1.00 10.00           C
ATOM   4950  O   LEU B 251     225.767 240.166 234.844  1.00 10.00           O
ATOM   4951  CB  LEU B 251     226.933 239.153 232.107  1.00 10.00           C
ATOM   4952  CG  LEU B 251     227.965 240.219 232.443  1.00 10.00           C
ATOM   4953  CD1 LEU B 251     228.397 240.185 233.875  1.00 10.00           C
ATOM   4954  CD2 LEU B 251     229.133 239.974 231.566  1.00 10.00           C
ATOM   4955  N   SER B 252     225.330 241.645 233.189  1.00 10.00           N
ATOM   4956  CA  SER B 252     225.347 242.803 234.075  1.00 10.00           C
ATOM   4957  C   SER B 252     224.316 242.658 235.188  1.00 10.00           C
ATOM   4958  O   SER B 252     224.646 242.763 236.359  1.00 10.00           O
ATOM   4959  CB  SER B 252     225.100 244.087 233.281  1.00 10.00           C
ATOM   4960  OG  SER B 252     223.934 243.978 232.483  1.00 10.00           O
ATOM   4961  N   ILE B 253     223.079 242.406 234.777  1.00 10.00           N
ATOM   4962  CA  ILE B 253     221.984 242.215 235.701  1.00 10.00           C
ATOM   4963  C   ILE B 253     222.318 241.079 236.674  1.00 10.00           C
ATOM   4964  O   ILE B 253     221.966 241.155 237.861  1.00 10.00           O
ATOM   4965  CB  ILE B 253     220.637 241.967 234.909  1.00 10.00           C
ATOM   4966  CG1 ILE B 253     219.873 243.313 234.717  1.00 10.00           C
ATOM   4967  CG2 ILE B 253     219.812 240.846 235.566  1.00 10.00           C
ATOM   4968  CD1 ILE B 253     218.952 243.807 235.863  1.00 10.00           C
ATOM   4969  N   SER B 254     223.006 240.035 236.210  1.00 10.00           N
ATOM   4970  CA  SER B 254     223.376 238.942 237.101  1.00 10.00           C
ATOM   4971  C   SER B 254     224.230 239.445 238.260  1.00 10.00           C
ATOM   4972  O   SER B 254     224.153 238.954 239.385  1.00 10.00           O
ATOM   4973  CB  SER B 254     224.117 237.850 236.327  1.00 10.00           C
ATOM   4974  OG  SER B 254     225.519 237.973 236.484  1.00 10.00           O
ATOM   4975  N   ALA B 255     225.086 240.395 237.894  1.00 10.00           N
ATOM   4976  CA  ALA B 255     226.037 241.048 238.788  1.00 10.00           C
ATOM   4977  C   ALA B 255     225.307 241.739 239.912  1.00 10.00           C
ATOM   4978  O   ALA B 255     225.642 241.558 241.078  1.00 10.00           O
ATOM   4979  CB  ALA B 255     226.895 242.082 238.010  1.00 10.00           C
ATOM   4980  N   LEU B 256     224.309 242.534 239.561  1.00 10.00           N
ATOM   4981  CA  LEU B 256     223.576 243.237 240.619  1.00 10.00           C
ATOM   4982  C   LEU B 256     222.698 242.283 241.429  1.00 10.00           C
ATOM   4983  O   LEU B 256     222.629 242.451 242.637  1.00 10.00           O
ATOM   4984  CB  LEU B 256     222.714 244.357 240.023  1.00 10.00           C
ATOM   4985  CG  LEU B 256     223.454 245.455 239.256  1.00 10.00           C
ATOM   4986  CD1 LEU B 256     222.475 246.489 238.720  1.00 10.00           C
ATOM   4987  CD2 LEU B 256     224.499 246.119 240.141  1.00 10.00           C
ATOM   4988  N   LEU B 257     222.022 241.325 240.800  1.00 10.00           N
ATOM   4989  CA  LEU B 257     221.351 240.362 241.696  1.00 10.00           C
ATOM   4990  C   LEU B 257     222.312 239.918 242.800  1.00 10.00           C
ATOM   4991  O   LEU B 257     222.010 239.941 243.992  1.00 10.00           O
ATOM   4992  CB  LEU B 257     220.902 239.113 240.928  1.00 10.00           C
ATOM   4993  CG  LEU B 257     219.677 239.334 240.039  1.00 10.00           C
ATOM   4994  CD1 LEU B 257     218.840 238.062 239.860  1.00 10.00           C
ATOM   4995  CD2 LEU B 257     218.716 240.394 240.581  1.00 10.00           C
ATOM   4996  N   ALA B 258     223.421 239.380 242.242  1.00 10.00           N
ATOM   4997  CA  ALA B 258     224.562 239.012 243.063  1.00 10.00           C
ATOM   4998  C   ALA B 258     225.004 240.169 243.948  1.00 10.00           C
ATOM   4999  O   ALA B 258     225.149 240.011 245.189  1.00 10.00           O
ATOM   5000  CB  ALA B 258     225.713 238.543 242.190  1.00 10.00           C
ATOM   5001  N   LEU B 259     225.248 241.342 243.377  1.00 10.00           N
ATOM   5002  CA  LEU B 259     225.691 242.450 244.206  1.00 10.00           C
ATOM   5003  C   LEU B 259     224.678 242.839 245.284  1.00 10.00           C
ATOM   5004  O   LEU B 259     225.062 243.145 246.415  1.00 10.00           O
ATOM   5005  CB  LEU B 259     226.049 243.676 243.341  1.00 10.00           C
ATOM   5006  CG  LEU B 259     227.151 244.631 243.850  1.00 10.00           C
ATOM   5007  CD1 LEU B 259     227.592 245.577 242.736  1.00 10.00           C
ATOM   5008  CD2 LEU B 259     226.652 245.413 245.048  1.00 10.00           C
ATOM   5009  N   THR B 260     223.399 242.829 244.948  1.00 10.00           N
ATOM   5010  CA  THR B 260     222.387 243.200 245.924  1.00 10.00           C
ATOM   5011  C   THR B 260     222.189 242.157 247.018  1.00 10.00           C
ATOM   5012  O   THR B 260     222.041 242.522 248.180  1.00 10.00           O
ATOM   5013  CB  THR B 260     221.034 243.494 245.252  1.00 10.00           C
ATOM   5014  OG1 THR B 260     221.208 244.520 244.268  1.00 10.00           O
ATOM   5015  CG2 THR B 260     220.017 243.972 246.288  1.00 10.00           C
ATOM   5016  N   VAL B 261     222.180 240.871 246.671  1.00 10.00           N
ATOM   5017  CA  VAL B 261     222.004 239.850 247.696  1.00 10.00           C
ATOM   5018  C   VAL B 261     223.159 239.978 248.701  1.00 10.00           C
ATOM   5019  O   VAL B 261     222.934 240.017 249.908  1.00 10.00           O
ATOM   5020  CB  VAL B 261     221.950 238.444 247.085  1.00 10.00           C
ATOM   5021  CG1 VAL B 261     223.254 237.708 247.343  1.00 10.00           C
ATOM   5022  CG2 VAL B 261     220.760 237.676 247.653  1.00 10.00           C
ATOM   5023  N   PHE B 262     224.382 240.077 248.191  1.00 10.00           N
ATOM   5024  CA  PHE B 262     225.569 240.246 249.018  1.00 10.00           C
ATOM   5025  C   PHE B 262     225.380 241.384 250.006  1.00 10.00           C
ATOM   5026  O   PHE B 262     225.689 241.251 251.207  1.00 10.00           O
ATOM   5027  CB  PHE B 262     226.797 240.499 248.145  1.00 10.00           C
ATOM   5028  CG  PHE B 262     227.078 239.399 247.169  1.00 10.00           C
ATOM   5029  CD1 PHE B 262     226.127 238.427 246.917  1.00 10.00           C
ATOM   5030  CD2 PHE B 262     228.290 239.340 246.505  1.00 10.00           C
ATOM   5031  CE1 PHE B 262     226.382 237.410 246.017  1.00 10.00           C
ATOM   5032  CE2 PHE B 262     228.554 238.326 245.601  1.00 10.00           C
ATOM   5033  CZ  PHE B 262     227.597 237.358 245.358  1.00 10.00           C
ATOM   5034  N   LEU B 263     225.005 242.562 249.524  1.00 10.00           N
ATOM   5035  CA  LEU B 263     224.853 243.689 250.414  1.00 10.00           C
ATOM   5036  C   LEU B 263     223.675 243.543 251.357  1.00 10.00           C
ATOM   5037  O   LEU B 263     223.652 244.130 252.450  1.00 10.00           O
ATOM   5038  CB  LEU B 263     224.762 244.980 249.601  1.00 10.00           C
ATOM   5039  CG  LEU B 263     226.093 245.735 249.527  1.00 10.00           C
ATOM   5040  CD1 LEU B 263     226.330 246.117 248.086  1.00 10.00           C
ATOM   5041  CD2 LEU B 263     226.123 246.953 250.472  1.00 10.00           C
ATOM   5042  N   LEU B 264     222.686 242.760 250.956  1.00 10.00           N
ATOM   5043  CA  LEU B 264     221.542 242.579 251.816  1.00 10.00           C
ATOM   5044  C   LEU B 264     221.924 241.605 252.925  1.00 10.00           C
ATOM   5045  O   LEU B 264     221.463 241.743 254.065  1.00 10.00           O
ATOM   5046  CB  LEU B 264     220.334 242.083 251.014  1.00 10.00           C
ATOM   5047  CG  LEU B 264     219.755 243.159 250.078  1.00 10.00           C
ATOM   5048  CD1 LEU B 264     218.653 242.565 249.216  1.00 10.00           C
ATOM   5049  CD2 LEU B 264     219.225 244.338 250.898  1.00 10.00           C
ATOM   5050  N   LEU B 265     222.774 240.622 252.620  1.00 10.00           N
ATOM   5051  CA  LEU B 265     223.222 239.669 253.669  1.00 10.00           C
ATOM   5052  C   LEU B 265     224.032 240.505 254.652  1.00 10.00           C
ATOM   5053  O   LEU B 265     223.988 240.274 255.858  1.00 10.00           O
ATOM   5054  CB  LEU B 265     224.166 238.569 253.144  1.00 10.00           C
ATOM   5055  CG  LEU B 265     224.020 237.888 251.787  1.00 10.00           C
ATOM   5056  CD1 LEU B 265     225.443 237.601 251.239  1.00 10.00           C
ATOM   5057  CD2 LEU B 265     223.153 236.632 251.913  1.00 10.00           C
ATOM   5058  N   LEU B 266     224.807 241.449 254.109  1.00 10.00           N
ATOM   5059  CA  LEU B 266     225.613 242.341 254.936  1.00 10.00           C
ATOM   5060  C   LEU B 266     224.650 243.101 255.827  1.00 10.00           C
ATOM   5061  O   LEU B 266     224.959 243.415 256.978  1.00 10.00           O
ATOM   5062  CB  LEU B 266     226.400 243.325 254.070  1.00 10.00           C
ATOM   5063  CG  LEU B 266     227.671 242.791 253.407  1.00 10.00           C
ATOM   5064  CD1 LEU B 266     228.218 243.833 252.445  1.00 10.00           C
ATOM   5065  CD2 LEU B 266     228.705 242.448 254.474  1.00 10.00           C
ATOM   5066  N   ALA B 267     223.483 243.402 255.267  1.00 10.00           N
ATOM   5067  CA  ALA B 267     222.436 244.089 255.989  1.00 10.00           C
ATOM   5068  C   ALA B 267     222.023 243.185 257.154  1.00 10.00           C
ATOM   5069  O   ALA B 267     221.667 243.655 258.244  1.00 10.00           O
ATOM   5070  CB  ALA B 267     221.285 244.310 255.082  1.00 10.00           C
ATOM   5071  N   ASP B 268     222.131 241.880 256.924  1.00 10.00           N
ATOM   5072  CA  ASP B 268     221.777 240.836 257.862  1.00 10.00           C
ATOM   5073  C   ASP B 268     222.713 240.754 259.077  1.00 10.00           C
ATOM   5074  O   ASP B 268     222.260 240.620 260.224  1.00 10.00           O
ATOM   5075  CB  ASP B 268     221.764 239.492 257.131  1.00 10.00           C
ATOM   5076  CG  ASP B 268     220.507 239.255 256.337  1.00 10.00           C
ATOM   5077  OD1 ASP B 268     219.399 239.480 256.867  1.00 10.00           O
ATOM   5078  OD2 ASP B 268     220.644 238.809 255.181  1.00 10.00           O
ATOM   5079  N   LYS B 269     224.017 240.864 258.839  1.00 10.00           N
ATOM   5080  CA  LYS B 269     224.982 240.732 259.921  1.00 10.00           C
ATOM   5081  C   LYS B 269     225.111 242.032 260.702  1.00 10.00           C
ATOM   5082  O   LYS B 269     225.405 241.979 261.910  1.00 10.00           O
ATOM   5083  CB  LYS B 269     226.344 240.301 259.377  1.00 10.00           C
ATOM   5084  CG  LYS B 269     226.919 241.244 258.332  1.00 10.00           C
ATOM   5085  CD  LYS B 269     227.543 242.471 258.979  1.00 10.00           C
ATOM   5086  CE  LYS B 269     227.861 243.539 257.944  1.00 10.00           C
ATOM   5087  NZ  LYS B 269     229.289 243.957 257.999  1.00 10.00           N
ATOM   5088  N   VAL B 270     224.951 243.183 260.068  1.00 10.00           N
ATOM   5089  CA  VAL B 270     225.061 244.484 260.753  1.00 10.00           C
ATOM   5090  C   VAL B 270     224.046 244.754 261.883  1.00 10.00           C
ATOM   5091  O   VAL B 270     224.421 245.251 262.949  1.00 10.00           O
ATOM   5092  CB  VAL B 270     224.953 245.662 259.753  1.00 10.00           C
ATOM   5093  CG1 VAL B 270     225.245 246.971 260.475  1.00 10.00           C
ATOM   5094  CG2 VAL B 270     225.921 245.468 258.594  1.00 10.00           C
ATOM   5095  N   PRO B 271     222.758 244.423 261.663  1.00 10.00           N
ATOM   5096  CA  PRO B 271     221.782 244.685 262.718  1.00 10.00           C
ATOM   5097  C   PRO B 271     221.942 243.734 263.890  1.00 10.00           C
ATOM   5098  O   PRO B 271     221.722 244.111 265.044  1.00 10.00           O
ATOM   5099  CB  PRO B 271     220.463 244.476 262.003  1.00 10.00           C
ATOM   5100  CG  PRO B 271     220.768 243.256 261.202  1.00 10.00           C
ATOM   5101  CD  PRO B 271     222.173 243.534 260.649  1.00 10.00           C
ATOM   5102  N   GLU B 272     222.323 242.506 263.582  1.00 10.00           N
ATOM   5103  CA  GLU B 272     222.482 241.519 264.642  1.00 10.00           C
ATOM   5104  C   GLU B 272     223.812 241.696 265.365  1.00 10.00           C
ATOM   5105  O   GLU B 272     224.014 241.153 266.456  1.00 10.00           O
ATOM   5106  CB  GLU B 272     222.365 240.101 264.078  1.00 10.00           C
ATOM   5107  CG  GLU B 272     222.283 239.016 265.139  1.00 10.00           C
ATOM   5108  CD  GLU B 272     222.192 237.626 264.543  1.00 10.00           C
ATOM   5109  OE1 GLU B 272     223.250 237.045 264.224  1.00 10.00           O
ATOM   5110  OE2 GLU B 272     221.061 237.115 264.397  1.00 10.00           O
ATOM   5111  N   THR B 273     224.746 242.448 264.783  1.00 10.00           N
ATOM   5112  CA  THR B 273     226.022 242.664 265.441  1.00 10.00           C
ATOM   5113  C   THR B 273     225.966 243.896 266.325  1.00 10.00           C
ATOM   5114  O   THR B 273     226.650 243.953 267.352  1.00 10.00           O
ATOM   5115  CB  THR B 273     227.158 242.804 264.429  1.00 10.00           C
ATOM   5116  OG1 THR B 273     226.789 243.755 263.422  1.00 10.00           O
ATOM   5117  CG2 THR B 273     227.462 241.462 263.782  1.00 10.00           C
ATOM   5118  N   SER B 274     225.184 244.913 265.951  1.00 10.00           N
ATOM   5119  CA  SER B 274     225.097 246.091 266.795  1.00 10.00           C
ATOM   5120  C   SER B 274     224.192 245.616 267.923  1.00 10.00           C
ATOM   5121  O   SER B 274     224.236 246.155 269.001  1.00 10.00           O
ATOM   5122  CB  SER B 274     224.482 247.284 266.043  1.00 10.00           C
ATOM   5123  OG  SER B 274     223.149 247.015 265.626  1.00 10.00           O
ATOM   5124  N   LEU B 275     223.413 244.582 267.657  1.00 10.00           N
ATOM   5125  CA  LEU B 275     222.480 243.986 268.618  1.00 10.00           C
ATOM   5126  C   LEU B 275     223.102 243.867 270.011  1.00 10.00           C
ATOM   5127  O   LEU B 275     222.364 243.617 270.970  1.00 10.00           O
ATOM   5128  CB  LEU B 275     222.019 242.605 268.137  1.00 10.00           C
ATOM   5129  CG  LEU B 275     221.032 242.583 266.968  1.00 10.00           C
ATOM   5130  CD1 LEU B 275     220.269 241.268 266.933  1.00 10.00           C
ATOM   5131  CD2 LEU B 275     220.068 243.757 267.056  1.00 10.00           C
ATOM   5132  N   SER B 276     224.404 243.994 270.159  1.00 10.00           N
ATOM   5133  CA  SER B 276     224.994 243.862 271.497  1.00 10.00           C
ATOM   5134  C   SER B 276     225.986 244.956 271.827  1.00 10.00           C
ATOM   5135  O   SER B 276     226.479 245.038 272.953  1.00 10.00           O
ATOM   5136  CB  SER B 276     225.691 242.501 271.672  1.00 10.00           C
ATOM   5137  OG  SER B 276     226.767 242.591 272.594  1.00 10.00           O
ATOM   5138  N   VAL B 277     226.289 245.806 270.856  1.00 10.00           N
ATOM   5139  CA  VAL B 277     227.258 246.856 271.114  1.00 10.00           C
ATOM   5140  C   VAL B 277     227.014 248.142 270.273  1.00 10.00           C
ATOM   5141  O   VAL B 277     227.956 248.785 269.847  1.00 10.00           O
ATOM   5142  CB  VAL B 277     228.758 246.256 270.941  1.00 10.00           C
ATOM   5143  CG1 VAL B 277     228.913 244.967 271.745  1.00 10.00           C
ATOM   5144  CG2 VAL B 277     229.100 245.960 269.475  1.00 10.00           C
ATOM   5145  N   PRO B 278     225.737 248.571 270.094  1.00 10.00           N
ATOM   5146  CA  PRO B 278     225.400 249.757 269.312  1.00 10.00           C
ATOM   5147  C   PRO B 278     225.702 251.104 269.972  1.00 10.00           C
ATOM   5148  O   PRO B 278     225.755 251.223 271.199  1.00 10.00           O
ATOM   5149  CB  PRO B 278     223.904 249.567 269.050  1.00 10.00           C
ATOM   5150  CG  PRO B 278     223.440 249.030 270.354  1.00 10.00           C
ATOM   5151  CD  PRO B 278     224.540 248.034 270.745  1.00 10.00           C
ATOM   5152  N   ILE B 279     225.950 252.080 269.159  1.00 10.00           N
ATOM   5153  CA  ILE B 279     226.182 253.473 269.507  1.00 10.00           C
ATOM   5154  C   ILE B 279     225.797 254.381 268.330  1.00 10.00           C
ATOM   5155  O   ILE B 279     225.022 253.864 267.510  1.00 10.00           O
ATOM   5156  CB  ILE B 279     227.633 253.732 269.931  1.00 10.00           C
ATOM   5157  CG1 ILE B 279     227.947 252.979 271.225  1.00 10.00           C
ATOM   5158  CG2 ILE B 279     227.875 255.224 270.100  1.00 10.00           C
ATOM   5159  CD1 ILE B 279     229.127 252.039 271.115  1.00 10.00           C
ATOM   5160  N   ILE B 280     226.353 255.538 268.171  1.00 10.00           N
ATOM   5161  CA  ILE B 280     225.965 256.355 267.010  1.00 10.00           C
ATOM   5162  C   ILE B 280     226.157 255.647 265.686  1.00 10.00           C
ATOM   5163  O   ILE B 280     227.025 254.789 265.543  1.00 10.00           O
ATOM   5164  CB  ILE B 280     226.687 257.771 266.906  1.00 10.00           C
ATOM   5165  CG1 ILE B 280     228.143 257.730 267.425  1.00 10.00           C
ATOM   5166  CG2 ILE B 280     225.779 258.833 267.568  1.00 10.00           C
ATOM   5167  CD1 ILE B 280     228.263 257.823 268.900  1.00 10.00           C
"""

def tst_to_cif_annotation():
  pass

def tst_remove_empty_annotations():
  ann = annotation.from_records(pdb_records_2.split("\n"))
  pdb_h = iotbx.pdb.input(source_info=None, lines=pdb_string_2.split('\n')).\
      construct_hierarchy()
  ann.remove_empty_annotations(hierarchy=pdb_h)
  # print ann
  assert ann.get_n_helices() == 2
  assert ann.get_n_sheets() == 0

def tst_split_helices_with_prolines():
  ann = annotation.from_records(pdb_records_2.split("\n"))
  pdb_h = iotbx.pdb.input(source_info=None, lines=pdb_string_2.split('\n')).\
      construct_hierarchy()
  # pdb_h.write_pdb_file(file_name='1.pdb')
  ann.remove_empty_annotations(hierarchy=pdb_h)
  # print ann
  ann.split_helices_with_prolines(hierarchy=pdb_h)
  # print "="*50
  # print ann
  assert ann.get_n_helices() == 4
  h_sizes = [x.length for x in ann.helices]
  assert h_sizes == [7, 15, 24, 5], h_sizes

def tst_split_helices_with_prolines_2():
  """ Corner case when the helix begins with PRO, 2nd residue is also PRO"""
  pro_ann = "HELIX   16  16 PRO B  463  ALA B  467  5                                   5"
  pdb_str = """\
ATOM   5910  N   PRO B 463     239.715 263.725 264.122  1.00 10.00           N
ATOM   5911  CA  PRO B 463     238.974 264.849 264.695  1.00 10.00           C
ATOM   5912  C   PRO B 463     238.267 264.611 266.069  1.00 10.00           C
ATOM   5913  O   PRO B 463     238.533 265.320 267.055  1.00 10.00           O
ATOM   5914  CB  PRO B 463     237.997 265.183 263.567  1.00 10.00           C
ATOM   5915  CG  PRO B 463     238.911 265.060 262.343  1.00 10.00           C
ATOM   5916  CD  PRO B 463     239.700 263.801 262.637  1.00 10.00           C
ATOM   5917  N   PRO B 464     237.359 263.619 266.160  1.00 10.00           N
ATOM   5918  CA  PRO B 464     236.693 263.407 267.459  1.00 10.00           C
ATOM   5919  C   PRO B 464     237.446 262.481 268.415  1.00 10.00           C
ATOM   5920  O   PRO B 464     236.830 261.688 269.109  1.00 10.00           O
ATOM   5921  CB  PRO B 464     235.345 262.820 267.063  1.00 10.00           C
ATOM   5922  CG  PRO B 464     235.738 261.949 265.934  1.00 10.00           C
ATOM   5923  CD  PRO B 464     236.751 262.765 265.124  1.00 10.00           C
ATOM   5924  N   ASP B 465     238.768 262.633 268.474  1.00 10.00           N
ATOM   5925  CA  ASP B 465     239.620 261.794 269.304  1.00 10.00           C
ATOM   5926  C   ASP B 465     240.116 262.429 270.606  1.00 10.00           C
ATOM   5927  O   ASP B 465     240.198 261.761 271.642  1.00 10.00           O
ATOM   5928  CB  ASP B 465     240.822 261.329 268.465  1.00 10.00           C
ATOM   5929  CG  ASP B 465     241.926 262.387 268.361  1.00 10.00           C
ATOM   5930  OD1 ASP B 465     241.627 263.557 268.042  1.00 10.00           O
ATOM   5931  OD2 ASP B 465     243.106 262.040 268.593  1.00 10.00           O
ATOM   5932  N   ASN B 466     240.420 263.721 270.572  1.00 10.00           N
ATOM   5933  CA  ASN B 466     240.960 264.382 271.760  1.00 10.00           C
ATOM   5934  C   ASN B 466     240.023 264.651 272.935  1.00 10.00           C
ATOM   5935  O   ASN B 466     240.482 265.044 273.998  1.00 10.00           O
ATOM   5936  CB  ASN B 466     241.654 265.686 271.349  1.00 10.00           C
ATOM   5937  CG  ASN B 466     240.831 266.505 270.376  1.00 10.00           C
ATOM   5938  OD1 ASN B 466     239.714 266.921 270.685  1.00 10.00           O
ATOM   5939  ND2 ASN B 466     241.380 266.739 269.189  1.00 10.00           N
ATOM   5940  N   ALA B 467     238.694 264.465 272.767  1.00 10.00           N
ATOM   5941  CA  ALA B 467     237.910 264.755 273.973  1.00 10.00           C
ATOM   5942  C   ALA B 467     238.085 263.704 275.100  1.00 10.00           C
ATOM   5943  O   ALA B 467     237.157 263.430 275.874  1.00 10.00           O
ATOM   5944  CB  ALA B 467     236.468 264.889 273.432  1.00 10.00           C
  """
  ann = annotation.from_records(pro_ann.split("\n"))
  pdb_h = iotbx.pdb.input(source_info=None, lines=pdb_str.split('\n')).\
      construct_hierarchy()
  assert ann.get_n_helices() == 1
  ann.split_helices_with_prolines(hierarchy=pdb_h)
  # print ann
  assert ann.get_n_helices() == 1
  assert ann.helices[0].length == 4

def tst_split_helices_with_prolines_3():
  """ Corner case when the helix ends with PRO"""
  pro_ann = "HELIX   16  16 THR B  223  PRO B  227  5                                   5"
  pdb_str = """\
ATOM   4732  N   TYR B 223     237.498 243.916 256.959  1.00 10.00           N
ATOM   4733  CA  TYR B 223     236.998 243.193 255.797  1.00 10.00           C
ATOM   4734  C   TYR B 223     235.642 243.731 255.357  1.00 10.00           C
ATOM   4735  O   TYR B 223     235.112 243.268 254.374  1.00 10.00           O
ATOM   4736  CB  TYR B 223     236.908 241.695 256.094  1.00 10.00           C
ATOM   4737  CG  TYR B 223     235.804 241.324 257.057  1.00 10.00           C
ATOM   4738  CD1 TYR B 223     234.480 241.272 256.641  1.00 10.00           C
ATOM   4739  CD2 TYR B 223     236.085 241.027 258.384  1.00 10.00           C
ATOM   4740  CE1 TYR B 223     233.468 240.932 257.518  1.00 10.00           C
ATOM   4741  CE2 TYR B 223     235.080 240.686 259.269  1.00 10.00           C
ATOM   4742  CZ  TYR B 223     233.774 240.641 258.831  1.00 10.00           C
ATOM   4743  OH  TYR B 223     232.772 240.302 259.711  1.00 10.00           O
ATOM   4744  N   THR B 224     235.127 244.742 256.068  1.00 10.00           N
ATOM   4745  CA  THR B 224     233.888 245.465 255.737  1.00 10.00           C
ATOM   4746  C   THR B 224     234.274 246.452 254.592  1.00 10.00           C
ATOM   4747  O   THR B 224     233.527 246.682 253.621  1.00 10.00           O
ATOM   4748  CB  THR B 224     233.390 246.267 256.985  1.00 10.00           C
ATOM   4749  OG1 THR B 224     232.815 245.365 257.938  1.00 10.00           O
ATOM   4750  CG2 THR B 224     232.338 247.298 256.604  1.00 10.00           C
ATOM   4751  N   ILE B 225     235.481 246.993 254.709  1.00 10.00           N
ATOM   4752  CA  ILE B 225     236.079 247.939 253.743  1.00 10.00           C
ATOM   4753  C   ILE B 225     236.203 247.150 252.435  1.00 10.00           C
ATOM   4754  O   ILE B 225     235.800 247.601 251.358  1.00 10.00           O
ATOM   4755  CB  ILE B 225     237.546 248.405 254.230  1.00 10.00           C
ATOM   4756  CG1 ILE B 225     237.487 249.823 254.873  1.00 10.00           C
ATOM   4757  CG2 ILE B 225     238.554 248.277 253.078  1.00 10.00           C
ATOM   4758  CD1 ILE B 225     238.802 250.372 255.483  1.00 10.00           C
ATOM   4759  N   VAL B 226     236.723 245.926 252.549  1.00 10.00           N
ATOM   4760  CA  VAL B 226     236.904 245.129 251.334  1.00 10.00           C
ATOM   4761  C   VAL B 226     235.606 245.025 250.571  1.00 10.00           C
ATOM   4762  O   VAL B 226     235.621 245.119 249.364  1.00 10.00           O
ATOM   4763  CB  VAL B 226     237.478 243.708 251.567  1.00 10.00           C
ATOM   4764  CG1 VAL B 226     238.787 243.804 252.297  1.00 10.00           C
ATOM   4765  CG2 VAL B 226     236.480 242.836 252.257  1.00 10.00           C
ATOM   4766  N   PRO B 227     234.471 244.784 251.253  1.00 10.00           N
ATOM   4767  CA  PRO B 227     233.225 244.717 250.484  1.00 10.00           C
ATOM   4768  C   PRO B 227     233.110 246.013 249.736  1.00 10.00           C
ATOM   4769  O   PRO B 227     233.251 246.069 248.522  1.00 10.00           O
ATOM   4770  CB  PRO B 227     232.142 244.666 251.550  1.00 10.00           C
ATOM   4771  CG  PRO B 227     232.785 243.987 252.657  1.00 10.00           C
ATOM   4772  CD  PRO B 227     234.206 244.487 252.661  1.00 10.00           C
"""
  ann = annotation.from_records(pro_ann.split("\n"))
  pdb_h = iotbx.pdb.input(source_info=None, lines=pdb_str.split('\n')).\
      construct_hierarchy()
  assert ann.get_n_helices() == 1
  ann.split_helices_with_prolines(hierarchy=pdb_h)
  # print ann
  assert ann.get_n_helices() == 1
  assert ann.helices[0].length == 4

def tst_remove_short_annotations():
  ann = annotation.from_records(pdb_records_2.split("\n"))
  assert ann.get_n_helices() == 7
  assert ann.get_n_sheets() == 4
  ann.remove_short_annotations()
  assert ann.get_n_helices() == 6, ann.get_n_helices()
  assert ann.get_n_sheets() == 3, ann.get_n_sheets()
  nstr = [x.n_strands for x in ann.sheets]
  # print ann
  assert nstr == [5,3,2], nstr

def tst_remove_3_10_helices():
  ann = annotation.from_records(pdb_records_2.split("\n"))
  assert ann.get_n_helices() == 7
  assert ann.get_n_sheets() == 4
  ann.remove_3_10_helices()
  assert ann.get_n_helices() == 6, ann.get_n_helices()
  assert ann.get_n_sheets() == 4, ann.get_n_sheets()
  # print ann

def tst_concatenate_consecutive_helices():
  ann_str = """\
HELIX  233 233 PHE Z   12  CYS Z   43  1                                  32
HELIX  234 234 LEU Z   49  ILE Z   54  1                                   6
HELIX  235 235 ILE Z   54  ILE Z   62  1                                   9
HELIX  236 236 ILE Z   62  SER Z   77  1                                  16
HELIX  237 237 ALA Z   83  GLN Z  122  1                                  40
HELIX  238 238 LEU Z  125  ALA Z  136  1                                  12
  """
  ann = annotation.from_records(ann_str.split("\n"))
  assert ann.get_n_helices() == 6
  assert ann.get_n_sheets() == 0
  ann.concatenate_consecutive_helices()
  assert ann.get_n_helices() == 4
  # print ann
  h_sizes = [x.length for x in ann.helices]
  assert h_sizes == [32, 29, 40, 12], h_sizes

def tst_concatenate_consecutive_helices2():
  """ Very tight turn between two consecutive helices results in their
  concatenation. """
  if (not libtbx.env.has_module(name="ksdssp")):
    print("Skipped: required module ksdssp not present")
    return
  file_path = libtbx.env.find_in_repositories(
    relative_path="cctbx_project/iotbx/regression/secondary_structure/5a63_chainCp.pdb",
    test=os.path.isfile)
  if (file_path is None):
    print('WARNING: Skipping tst_concatenate_consecutive_helices2("%s"): input file not available' % file_path)
    return
  inp = iotbx.pdb.input(file_name=file_path)
  original_ann = inp.extract_secondary_structure()
  assert original_ann.get_n_helices() == 2
  assert original_ann.get_n_sheets() == 0
  assert original_ann.helices[0].get_start_resseq_as_int() == 155
  assert original_ann.helices[0].get_end_resseq_as_int() == 185
  assert original_ann.helices[1].get_start_resseq_as_int() == 186
  assert original_ann.helices[1].get_end_resseq_as_int() == 206
  h = inp.construct_hierarchy()
  ann = original_ann.deep_copy()
  ann.concatenate_consecutive_helices(hierarchy=h)
  assert ann.get_n_helices() == 2
  assert ann.get_n_sheets() == 0
  # and make sure they got 'spreaded'
  assert ann.helices[0].get_start_resseq_as_int() == 155
  assert ann.helices[0].get_end_resseq_as_int() == 184
  assert ann.helices[1].get_start_resseq_as_int() == 187
  assert ann.helices[1].get_end_resseq_as_int() == 206
  # print ann

def tst_concatenate_consecutive_helices3():
  """ Very tight turn between two consecutive helices. They are separated
  by PRO residue, but this is not enough, they should be spreaded like in
  previous test. """
  if (not libtbx.env.has_module(name="ksdssp")):
    print("Skipped: required module ksdssp not present")
    return
  file_path = libtbx.env.find_in_repositories(
    relative_path="cctbx_project/iotbx/regression/secondary_structure/5a63_chainBp.pdb",
    test=os.path.isfile)
  if (file_path is None):
    print('WARNING: Skipping tst_concatenate_consecutive_helices2("%s"): input file not available' % file_path)
    return
  inp = iotbx.pdb.input(file_name=file_path)
  original_ann = inp.extract_secondary_structure()
  assert original_ann.get_n_helices() == 2
  assert original_ann.get_n_sheets() == 0
  assert original_ann.helices[0].get_start_resseq_as_int() == 218
  assert original_ann.helices[0].get_end_resseq_as_int() == 241
  assert original_ann.helices[1].get_start_resseq_as_int() == 242
  assert original_ann.helices[1].get_end_resseq_as_int() == 260
  h = inp.construct_hierarchy()
  ann = original_ann.deep_copy()
  ann.concatenate_consecutive_helices(hierarchy=h)
  assert ann.get_n_helices() == 2
  assert ann.get_n_sheets() == 0
  # and make sure they got 'spreaded'
  assert ann.helices[0].get_start_resseq_as_int() == 218
  assert ann.helices[0].get_end_resseq_as_int() == 240
  assert ann.helices[1].get_start_resseq_as_int() == 243
  assert ann.helices[1].get_end_resseq_as_int() == 260

def tst_concatenate_consecutive_helices4():
  pdb_str = """\
HELIX  130 AO4 ASP i  421  GLY i  433  1                                  25
HELIX  131 AO5 GLY i  434  UNK i  441  1                                   8
ATOM  29887  N   ASP i 421     299.903 227.017 234.954  1.00 60.00           N
ATOM  29888  CA  ASP i 421     298.473 226.790 234.820  1.00 60.00           C
ATOM  29889  C   ASP i 421     298.155 226.535 233.351  1.00 60.00           C
ATOM  29890  O   ASP i 421     297.473 225.562 233.002  1.00 60.00           O
ATOM  29891  CB  ASP i 421     297.680 227.986 235.352  1.00 60.00           C
ATOM  29892  CG  ASP i 421     296.231 227.646 235.639  1.00 60.00           C
ATOM  29893  OD1 ASP i 421     295.733 226.637 235.095  1.00 60.00           O
ATOM  29894  OD2 ASP i 421     295.586 228.390 236.409  1.00 60.00           O1-
ATOM  29895  N   ASP i 422     298.673 227.411 232.492  1.00 60.00           N
ATOM  29896  CA  ASP i 422     298.476 227.272 231.053  1.00 60.00           C
ATOM  29897  C   ASP i 422     299.031 225.943 230.539  1.00 60.00           C
ATOM  29898  O   ASP i 422     298.395 225.273 229.722  1.00 60.00           O
ATOM  29899  CB  ASP i 422     299.132 228.437 230.309  1.00 60.00           C
ATOM  29900  CG  ASP i 422     298.563 229.783 230.718  1.00 60.00           C
ATOM  29901  OD1 ASP i 422     297.410 229.825 231.194  1.00 60.00           O
ATOM  29902  OD2 ASP i 422     299.274 230.798 230.560  1.00 60.00           O1-
ATOM  29903  N   MET i 423     300.208 225.559 231.027  1.00 60.00           N
ATOM  29904  CA  MET i 423     300.834 224.314 230.589  1.00 60.00           C
ATOM  29905  C   MET i 423     300.122 223.092 231.164  1.00 60.00           C
ATOM  29906  O   MET i 423     300.167 222.012 230.579  1.00 60.00           O
ATOM  29907  CB  MET i 423     302.314 224.291 230.974  1.00 60.00           C
ATOM  29908  CG  MET i 423     303.096 225.488 230.468  1.00 60.00           C
ATOM  29909  SD  MET i 423     302.799 225.820 228.722  1.00 60.00           S
ATOM  29910  CE  MET i 423     302.434 227.573 228.771  1.00 60.00           C
ATOM  29911  N   GLN i 424     299.477 223.265 232.315  1.00 60.00           N
ATOM  29912  CA  GLN i 424     298.638 222.217 232.893  1.00 60.00           C
ATOM  29913  C   GLN i 424     297.452 221.966 231.973  1.00 60.00           C
ATOM  29914  O   GLN i 424     297.170 220.820 231.589  1.00 60.00           O
ATOM  29915  CB  GLN i 424     298.157 222.611 234.295  1.00 60.00           C
ATOM  29916  CG  GLN i 424     297.207 221.609 234.945  1.00 60.00           C
ATOM  29917  CD  GLN i 424     296.714 222.061 236.312  1.00 60.00           C
ATOM  29918  OE1 GLN i 424     297.059 221.470 237.336  1.00 60.00           O
ATOM  29919  NE2 GLN i 424     295.895 223.107 236.333  1.00 60.00           N
ATOM  29920  N   ARG i 425     296.767 223.051 231.620  1.00 60.00           N
ATOM  29921  CA  ARG i 425     295.650 222.977 230.687  1.00 60.00           C
ATOM  29922  C   ARG i 425     296.091 222.293 229.399  1.00 60.00           C
ATOM  29923  O   ARG i 425     295.457 221.346 228.945  1.00 60.00           O
ATOM  29924  CB  ARG i 425     295.092 224.370 230.377  1.00 60.00           C
ATOM  29925  CG  ARG i 425     293.590 224.512 230.605  1.00 60.00           C
ATOM  29926  CD  ARG i 425     292.980 225.574 229.700  1.00 60.00           C
ATOM  29927  NE  ARG i 425     291.868 226.275 230.339  1.00 60.00           N
ATOM  29928  CZ  ARG i 425     291.310 227.385 229.870  1.00 60.00           C
ATOM  29929  NH1 ARG i 425     291.756 227.932 228.749  1.00 60.00           N1+
ATOM  29930  NH2 ARG i 425     290.306 227.953 230.524  1.00 60.00           N
ATOM  29931  N   MET i 426     297.200 222.769 228.837  1.00 60.00           N
ATOM  29932  CA  MET i 426     297.726 222.249 227.576  1.00 60.00           C
ATOM  29933  C   MET i 426     298.063 220.762 227.672  1.00 60.00           C
ATOM  29934  O   MET i 426     297.832 220.003 226.728  1.00 60.00           O
ATOM  29935  CB  MET i 426     298.965 223.045 227.158  1.00 60.00           C
ATOM  29936  CG  MET i 426     299.272 222.986 225.672  1.00 60.00           C
ATOM  29937  SD  MET i 426     297.903 223.582 224.659  1.00 60.00           S
ATOM  29938  CE  MET i 426     297.647 225.215 225.350  1.00 60.00           C
ATOM  29939  N   MET i 427     298.608 220.357 228.815  1.00 60.00           N
ATOM  29940  CA  MET i 427     298.910 218.957 229.074  1.00 60.00           C
ATOM  29941  C   MET i 427     297.629 218.130 229.045  1.00 60.00           C
ATOM  29942  O   MET i 427     297.599 217.034 228.481  1.00 60.00           O
ATOM  29943  CB  MET i 427     299.623 218.797 230.420  1.00 60.00           C
ATOM  29944  CG  MET i 427     299.974 217.361 230.776  1.00 60.00           C
ATOM  29945  SD  MET i 427     301.027 216.576 229.543  1.00 60.00           S
ATOM  29946  CE  MET i 427     302.516 217.559 229.687  1.00 60.00           C
ATOM  29947  N   LYS i 428     296.566 218.668 229.639  1.00 60.00           N
ATOM  29948  CA  LYS i 428     295.287 217.960 229.648  1.00 60.00           C
ATOM  29949  C   LYS i 428     294.684 217.856 228.244  1.00 60.00           C
ATOM  29950  O   LYS i 428     294.204 216.790 227.846  1.00 60.00           O
ATOM  29951  CB  LYS i 428     294.293 218.638 230.591  1.00 60.00           C
ATOM  29952  CG  LYS i 428     292.995 217.862 230.757  1.00 60.00           C
ATOM  29953  CD  LYS i 428     293.268 216.449 231.255  1.00 60.00           C
ATOM  29954  CE  LYS i 428     292.556 215.408 230.404  1.00 60.00           C
ATOM  29955  NZ  LYS i 428     292.980 214.020 230.745  1.00 60.00           N1+
ATOM  29956  N   LYS i 429     294.712 218.960 227.504  1.00 60.00           N
ATOM  29957  CA  LYS i 429     294.172 218.989 226.148  1.00 60.00           C
ATOM  29958  C   LYS i 429     294.929 218.009 225.261  1.00 60.00           C
ATOM  29959  O   LYS i 429     294.346 217.366 224.388  1.00 60.00           O
ATOM  29960  CB  LYS i 429     294.242 220.400 225.555  1.00 60.00           C
ATOM  29961  CG  LYS i 429     293.617 221.490 226.417  1.00 60.00           C
ATOM  29962  CD  LYS i 429     292.188 221.161 226.813  1.00 60.00           C
ATOM  29963  CE  LYS i 429     291.587 222.262 227.676  1.00 60.00           C
ATOM  29964  NZ  LYS i 429     290.229 221.908 228.178  1.00 60.00           N1+
ATOM  29965  N   MET i 430     296.234 217.901 225.494  1.00 60.00           N
ATOM  29966  CA  MET i 430     297.067 216.954 224.765  1.00 60.00           C
ATOM  29967  C   MET i 430     296.716 215.527 225.165  1.00 60.00           C
ATOM  29968  O   MET i 430     296.753 214.611 224.342  1.00 60.00           O
ATOM  29969  CB  MET i 430     298.549 217.227 225.033  1.00 60.00           C
ATOM  29970  CG  MET i 430     299.494 216.622 224.008  1.00 60.00           C
ATOM  29971  SD  MET i 430     299.571 217.597 222.493  1.00 60.00           S
ATOM  29972  CE  MET i 430     300.171 219.160 223.133  1.00 60.00           C
ATOM  29973  N   LYS i 431     296.369 215.348 226.436  1.00 60.00           N
ATOM  29974  CA  LYS i 431     296.049 214.027 226.968  1.00 60.00           C
ATOM  29975  C   LYS i 431     294.744 213.476 226.401  1.00 60.00           C
ATOM  29976  O   LYS i 431     294.664 212.306 226.030  1.00 60.00           O
ATOM  29977  CB  LYS i 431     295.967 214.077 228.496  1.00 60.00           C
ATOM  29978  N   LYS i 432     293.725 214.326 226.335  1.00 60.00           N
ATOM  29979  CA  LYS i 432     292.386 213.879 225.955  1.00 60.00           C
ATOM  29980  C   LYS i 432     292.198 213.671 224.449  1.00 60.00           C
ATOM  29981  O   LYS i 432     291.551 212.713 224.030  1.00 60.00           O
ATOM  29982  CB  LYS i 432     291.343 214.877 226.461  1.00 60.00           C
ATOM  29983  CG  LYS i 432     291.590 216.291 225.986  1.00 60.00           C
ATOM  29984  CD  LYS i 432     290.369 217.161 226.134  1.00 60.00           C
ATOM  29985  CE  LYS i 432     290.549 218.436 225.337  1.00 60.00           C
ATOM  29986  NZ  LYS i 432     289.514 219.436 225.680  1.00 60.00           N1+
ATOM  29987  N   GLY i 433     292.750 214.573 223.643  1.00 60.00           N
ATOM  29988  CA  GLY i 433     292.514 214.550 222.211  1.00 60.00           C
ATOM  29989  C   GLY i 433     293.553 213.806 221.400  1.00 60.00           C
ATOM  29990  O   GLY i 433     294.685 213.618 221.845  1.00 60.00           O
ATOM  29991  N   GLY i 434     293.168 213.376 220.205  1.00 60.00           N
ATOM  29992  CA  GLY i 434     294.107 212.724 219.309  1.00 60.00           C
ATOM  29993  C   GLY i 434     295.180 213.712 218.931  1.00 60.00           C
ATOM  29994  O   GLY i 434     294.922 214.909 218.843  1.00 60.00           O
TER
END
  """
  pdb_inp = iotbx.pdb.input(lines=pdb_str, source_info=None)
  h = pdb_inp.construct_hierarchy()
  ss_ann = pdb_inp.extract_secondary_structure()
  assert ss_ann.get_n_helices() == 2
  ss_ann.concatenate_consecutive_helices(hierarchy=h)
  assert ss_ann.get_n_helices() == 1
  # print ss_ann
  assert ss_ann.helices[0].get_start_resseq_as_int() == 421
  assert ss_ann.helices[0].get_end_resseq_as_int() == 433

def tst_simple_elements():
  ann_str = """\
HELIX    3 AA3 SER B   26  LYS B   30  5                                   5
HELIX    4 AA4 LEU B  196  GLU B  199  5                                   4
SHEET    1 AA1 4 ILE A   7  PRO A  10  0
SHEET    2 AA1 4 MET A  92  THR A  99  1  O  LEU A  95   N  ILE A   7
SHEET    3 AA1 4 THR A  73  SER A  82 -1  N  LEU A  76   O  ILE A  94
SHEET    4 AA1 4 VAL A  34  THR A  39 -1  N  PHE A  35   O  VAL A  81
  """
  ann = annotation.from_records(ann_str.split("\n"))
  assert ann.get_n_helices() == 2
  assert ann.get_n_sheets() == 1
  n_h = 0
  n_st = 0
  se = ann.simple_elements()
  assert len(se) == 6
  for e in se:
    if isinstance(e, pdb_helix):
      n_h += 1
    if isinstance(e, pdb_strand):
      n_st += 1
  assert n_h == 2
  assert n_st == 4
  s = ann.as_pdb_str()

def tst_filter_sheets_with_long_hbonds():
  """ Attention. In this case there are 2 almost identical SHEET annotations."""
  if (not libtbx.env.has_module(name="ksdssp")):
    print("Skipped: required module ksdssp not present")
    return
  file_path = libtbx.env.find_in_repositories(
    relative_path="cctbx_project/iotbx/regression/secondary_structure/3jd6_noh.pdb",
    test=os.path.isfile)
  if (file_path is None):
    print('WARNING: Skipping tst_filter_sheets_with_long_hbonds("%s"): input file not available' % file_path)
    return
  inp = iotbx.pdb.input(file_name=file_path)
  original_ann = inp.extract_secondary_structure()
  assert original_ann.get_n_helices() == 1
  assert original_ann.get_n_sheets() == 2
  h = inp.construct_hierarchy()
  ann = original_ann.deep_copy()
  ann.filter_sheets_with_long_hbonds(hierarchy=h)
  assert ann.get_n_helices() == 1, ann.get_n_helices()
  assert ann.get_n_sheets() == 3, ann.get_n_sheets() # with filtering out short ones

def tst_filter_sheets_with_long_hbonds2():
  """ bug found in 4a7h where ksdssp defines sheet 200 residues long sheet for
  rather small selection.
  A 2 ILE C 384  ARG C 586  1  N  LEU C 385   O  VAL C  97
  correct definition can be provided by any other method in ss_manager. """
  if (not libtbx.env.has_module(name="ksdssp")):
    print("Skipped: required module ksdssp not present")
    return
  file_path = libtbx.env.find_in_repositories(
    relative_path="cctbx_project/iotbx/regression/secondary_structure/4a7h_chainC.pdb",
    test=os.path.isfile)
  if (file_path is None):
    print('WARNING: Skipping tst_filter_sheets_with_long_hbonds2("%s"): input file not available' % file_path)
    return
  inp = iotbx.pdb.input(file_name=file_path)
  original_ann = inp.extract_secondary_structure()
  assert original_ann.get_n_helices() == 0
  assert original_ann.get_n_sheets() == 3
  h = inp.construct_hierarchy()
  ann = original_ann.deep_copy()
  try:
    ann.filter_sheets_with_long_hbonds(hierarchy=h)
  except Exception as e:
    if e.args[0].startswith("It is 4a7h"):
      pass
    else:
      raise e
  else:
    # disable for now
    pass
    # assert 0

def tst_filter_sheets_with_long_hbonds3():
  """ bug found in 1qmo: want to remove all sheets, but they appear not sorted."""
  if (not libtbx.env.has_module(name="ksdssp")):
    print("Skipped: required module ksdssp not present")
    return
  file_path = libtbx.env.find_in_repositories(
    relative_path="cctbx_project/iotbx/regression/secondary_structure/1ubf_cutted.pdb",
    test=os.path.isfile)
  if (file_path is None):
    print('WARNING: Skipping tst_filter_sheets_with_long_hbonds3("%s"): input file not available' % file_path)
    return
  inp = iotbx.pdb.input(file_name=file_path)
  original_ann = inp.extract_secondary_structure()
  assert original_ann.get_n_helices() == 11
  assert original_ann.get_n_sheets() == 2
  h = inp.construct_hierarchy()
  ann = original_ann.deep_copy()
  ann.filter_sheets_with_long_hbonds(hierarchy=h)
  # print(ann)
  assert ann.get_n_helices() == 11, ann.get_n_helices()
  assert ann.get_n_sheets() == 2, ann.get_n_sheets()

def tst_reset_sheet_ids():
  ann_str = """\
SHEET    1 AA1 4 ILE A   7  PRO A  10  0
SHEET    2 AA1 4 MET A  92  THR A  99  1  O  LEU A  95   N  ILE A   7
  """
  ann_str2 = """\
SHEET    1 AA1 4 ILE A   7  PRO A  10  0
SHEET    2 AA2 4 MET A  92  THR A  99  0
  """
  # One sheet
  ann = annotation.from_records(ann_str.split("\n"))
  assert ann.get_n_helices() == 0
  assert ann.get_n_sheets() == 1
  assert ann.sheets[0].sheet_id == 'AA1'
  ann.reset_sheet_ids()
  assert ann.get_n_helices() == 0
  assert ann.get_n_sheets() == 1
  assert ann.sheets[0].sheet_id == '0', ann.sheets[0].sheet_id

  # Two sheets
  ann = annotation.from_records(ann_str2.split("\n"))
  assert ann.get_n_helices() == 0
  assert ann.get_n_sheets() == 2
  assert ann.sheets[0].sheet_id == 'AA1'
  assert ann.sheets[1].sheet_id == 'AA2'
  ann.reset_sheet_ids()
  assert ann.get_n_helices() == 0
  assert ann.get_n_sheets() == 2
  assert ann.sheets[0].sheet_id == '0', ann.sheets[0].sheet_id
  assert ann.sheets[1].sheet_id == '1', ann.sheets[0].sheet_id

  # 11 sheets
  ann = annotation.from_records(ann_str.split("\n"))
  for x in range(10):
    ann.sheets.append(ann.sheets[0].deep_copy())
  assert ann.get_n_sheets() == 11
  ann.reset_sheet_ids()
  sheet_ids = [x.sheet_id for x in ann.sheets]
  assert sheet_ids ==['0', '1', '2', '3', '4', '5', '6', '7', '8', '9', 'A'], sheet_ids

  # 51 sheets
  ann = annotation.from_records(ann_str.split("\n"))
  for x in range(50):
    ann.sheets.append(ann.sheets[0].deep_copy())
  assert ann.get_n_sheets() == 51
  ann.reset_sheet_ids()
  sheet_ids = [x.sheet_id for x in ann.sheets]
  assert sheet_ids ==['00', '01', '02', '03', '04', '05', '06', '07', '08', '09',
      '0A', '0B', '0C', '0D', '0E', '0F', '0G', '0H', '0I', '0J', '0K', '0L',
      '0M', '0N', '0O', '0P', '0Q', '0R', '0S', '0T', '0U', '0V', '0W', '0X',
      '0Y', '0Z', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19',
      '1A', '1B', '1C', '1D', '1E'], sheet_ids

def tst_struct_conf_various():
  struct_conf_str = """
data_test
loop_
_struct_conf.beg_auth_asym_id
_struct_conf.beg_auth_comp_id
_struct_conf.beg_auth_seq_id
_struct_conf.beg_label_asym_id
_struct_conf.beg_label_comp_id
_struct_conf.beg_label_seq_id
_struct_conf.conf_type_id
_struct_conf.end_auth_asym_id
_struct_conf.end_auth_comp_id
_struct_conf.end_auth_seq_id
_struct_conf.end_label_asym_id
_struct_conf.end_label_comp_id
_struct_conf.end_label_seq_id
_struct_conf.id
_struct_conf.pdbx_beg_PDB_ins_code
_struct_conf.pdbx_end_PDB_ins_code
A LEU 8    A LEU 8    HELX_RH_AL_P A ARG 14   A ARG 14   HELX_RH_AL_P1  ? ?
A ASP 16   A ASP 16   HELX_RH_3T_P A GLU 18   A GLU 18   HELX_RH_3T_P1  ? ?
A PRO 20   A PRO 20   TURN_TY1_P   A GLY 21   A GLY 21   TURN_TY1_P1    ? ?
A THR 23   A THR 23   BEND         A THR 23   A THR 23   BEND1          ? ?
A GLN 24   A GLN 24   STRN         A LEU 25   A LEU 25   STRN1          ? ?
A ASN 26   A ASN 26   BEND         A ASN 26   A ASN 26   BEND2          ? ?
A ARG 27   A ARG 27   HELX_LH_PP_P A ALA 30   A ALA 30   HELX_LH_PP_P1  ? ?
A PRO 32   A PRO 32   BEND         A PRO 32   A PRO 32   BEND3          ? ?
A TRP 37   A TRP 37   STRN         A TRP 37   A TRP 37   STRN2          ? ?
A ARG 38   A ARG 38   BEND         A ASN 39   A ASN 39   BEND4          ? ?
A SER 40   A SER 40   HELX_RH_AL_P A THR 45   A THR 45   HELX_RH_AL_P2  ? ?
A ASP 46   A ASP 46   TURN_TY1_P   A ASP 46   A ASP 46   TURN_TY1_P2    ? ?
A ARG 47   A ARG 47   HELX_LH_PP_P A SER 49   A SER 49   HELX_LH_PP_P2  ? ?
A GLN 50   A GLN 50   TURN_TY1_P   A GLN 51   A GLN 51   TURN_TY1_P3    ? ?
A LEU 52   A LEU 52   STRN         A SER 54   A SER 54   STRN3          ? ?
A GLY 57   A GLY 57   STRN         A PHE 64   A PHE 64   STRN4          ? ?
A PRO 65   A PRO 65   BEND         A ALA 66   A ALA 66   BEND5          ? ?
A PRO 67   A PRO 67   HELX_RH_3T_P A ALA 69   A ALA 69   HELX_RH_3T_P2  ? ?
A VAL 70   A VAL 70   BEND         A VAL 70   A VAL 70   BEND6          ? ?
A GLU 72   A GLU 72   HELX_RH_3T_P A LEU 75   A LEU 75   HELX_RH_3T_P3  ? ?
A GLU 76   A GLU 76   TURN_TY1_P   A GLU 76   A GLU 76   TURN_TY1_P4    ? ?
A CYS 77   A CYS 77   BEND         A CYS 77   A CYS 77   BEND7          ? ?
A PRO 80   A PRO 80   TURN_TY1_P   A GLU 81   A GLU 81   TURN_TY1_P5    ? ?
A ALA 82   A ALA 82   HELX_LH_PP_P A ALA 82   A ALA 82   HELX_LH_PP_P3  ? ?
A ASP 83   A ASP 83   STRN         A VAL 87   A VAL 87   STRN5          ? ?
A PRO 88   A PRO 88   BEND         A PRO 88   A PRO 88   BEND8          ? ?
A TRP 91   A TRP 91   HELX_RH_3T_P A HIS 94   A HIS 94   HELX_RH_3T_P4  ? ?
A GLY 95   A GLY 95   TURN_TY1_P   A GLY 95   A GLY 95   TURN_TY1_P6    ? ?
A TYR 96   A TYR 96   BEND         A ASP 97   A ASP 97   BEND9          ? ?
A ILE 100  A ILE 100  STRN         A THR 102  A THR 102  STRN6          ? ?
A ASN 103  A ASN 103  BEND         A THR 105  A THR 105  BEND10         ? ?
A PRO 107  A PRO 107  BEND         A ILE 108  A ILE 108  BEND11         ? ?
A PRO 112  A PRO 112  TURN_TY1_P   A PRO 113  A PRO 113  TURN_TY1_P7    ? ?
A THR 117  A THR 117  BEND         A GLU 118  A GLU 118  BEND12         ? ?
A THR 121  A THR 121  STRN         A PHE 128  A PHE 128  STRN7          ? ?
A GLU 132  A GLU 132  HELX_RH_AL_P A GLN 136  A GLN 136  HELX_RH_AL_P3  ? ?
A GLU 137  A GLU 137  BEND         A GLY 138  A GLY 138  BEND13         ? ?
A GLN 139  A GLN 139  STRN         A PHE 144  A PHE 144  STRN8          ? ?
A ASP 145  A ASP 145  BEND         A ASP 145  A ASP 145  BEND14         ? ?
A VAL 147  A VAL 147  STRN         A ASN 148  A ASN 148  STRN9          ? ?
A SER 149  A SER 149  BEND         A SER 149  A SER 149  BEND15         ? ?
A ALA 150  A ALA 150  STRN         A CYS 155  A CYS 155  STRN10         ? ?
A ASN 156  A ASN 156  TURN_TY1_P   A GLY 157  A GLY 157  TURN_TY1_P8    ? ?
A ARG 158  A ARG 158  STRN         A GLN 164  A GLN 164  STRN11         ? ?
A SER 166  A SER 166  TURN_TY1_P   A ARG 167  A ARG 167  TURN_TY1_P9    ? ?
A LEU 168  A LEU 168  BEND         A LEU 168  A LEU 168  BEND16         ? ?
A SER 170  A SER 170  STRN         A ASP 173  A ASP 173  STRN12         ? ?
A SER 175  A SER 175  TURN_TY1_P   A PHE 177  A PHE 177  TURN_TY1_P10   ? ?
A ALA 180  A ALA 180  BEND         A GLY 181  A GLY 181  BEND17         ? ?
A ASN 183  A ASN 183  STRN         A LEU 190  A LEU 190  STRN13         ? ?
A ARG 191  A ARG 191  BEND         A ARG 191  A ARG 191  BEND18         ? ?
A ASP 194  A ASP 194  HELX_RH_3T_P A LEU 198  A LEU 198  HELX_RH_3T_P5  ? ?
A GLU 199  A GLU 199  STRN         A GLU 199  A GLU 199  STRN14         ? ?
A ASP 202  A ASP 202  BEND         A ASP 202  A ASP 202  BEND19         ? ?
A MET 203  A MET 203  STRN         A ARG 205  A ARG 205  STRN15         ? ?
A GLY 208  A GLY 208  STRN         A ILE 209  A ILE 209  STRN16         ? ?
A ARG 211  A ARG 211  BEND         A ARG 211  A ARG 211  BEND20         ? ?
A VAL 213  A VAL 213  STRN         A LYS 218  A LYS 218  STRN17         ? ?
A THR 220  A THR 220  BEND         A THR 221  A THR 221  BEND21         ? ?
A GLN 222  A GLN 222  STRN         A PHE 232  A PHE 232  STRN18         ? ?
A ASP 234  A ASP 234  TURN_TY1_P   A ASP 235  A ASP 235  TURN_TY1_P11   ? ?
A PHE 236  A PHE 236  BEND         A SER 237  A SER 237  BEND22         ? ?
A ARG 238  A ARG 238  STRN         A CYS 248  A CYS 248  STRN19         ? ?
A GLY 249  A GLY 249  BEND         A GLY 249  A GLY 249  BEND23         ? ?
A ASP 253  A ASP 253  TURN_TY1_P   A TYR 254  A TYR 254  TURN_TY1_P12   ? ?
A LEU 255  A LEU 255  STRN         A GLN 263  A GLN 263  STRN20         ? ?
A GLY 264  A GLY 264  TURN_TY1_P   A GLU 265  A GLU 265  TURN_TY1_P13   ? ?
A THR 266  A THR 266  STRN         A ALA 273  A ALA 273  STRN21         ? ?
A PHE 275  A PHE 275  BEND         A PHE 275  A PHE 275  BEND24         ? ?
A GLY 277  A GLY 277  HELX_LH_PP_P A ILE 279  A ILE 279  HELX_LH_PP_P4  ? ?
A ASP 281  A ASP 281  STRN         A ASP 281  A ASP 281  STRN22         ? ?
A GLU 282  A GLU 282  TURN_TY1_P   A ARG 283  A ARG 283  TURN_TY1_P14   ? ?
A GLY 284  A GLY 284  STRN         A GLY 284  A GLY 284  STRN23         ? ?
A ALA 287  A ALA 287  TURN_TY1_P   A ASP 288  A ASP 288  TURN_TY1_P15   ? ?
A ARG 289  A ARG 289  STRN         A GLU 297  A GLU 297  STRN24         ? ?
A ASN 298  A ASN 298  BEND         A ASN 298  A ASN 298  BEND25         ? ?
A LEU 301  A LEU 301  STRN         A TRP 302  A TRP 302  STRN25         ? ?
A ALA 304  A ALA 304  BEND         A ILE 306  A ILE 306  BEND26         ? ?
A PRO 307  A PRO 307  STRN         A PRO 307  A PRO 307  STRN26         ? ?
A TYR 310  A TYR 310  STRN         A THR 318  A THR 318  STRN27         ? ?
A ALA 319  A ALA 319  TURN_TY1_P   A ASP 320  A ASP 320  TURN_TY1_P16   ? ?
A GLY 321  A GLY 321  BEND         A GLY 321  A GLY 321  BEND27         ? ?
A LEU 323  A LEU 323  STRN         A VAL 331  A VAL 331  STRN28         ? ?
A PHE 333  A PHE 333  STRN         A PHE 333  A PHE 333  STRN29         ? ?
A VAL 336  A VAL 336  STRN         A GLU 339  A GLU 339  STRN30         ? ?
A ASN 340  A ASN 340  TURN_TY1_P   A GLY 341  A GLY 341  TURN_TY1_P17   ? ?
A LEU 342  A LEU 342  STRN         A LEU 345  A LEU 345  STRN31         ? ?
A ASN 346  A ASN 346  TURN_TY1_P   A GLY 347  A GLY 347  TURN_TY1_P18   ? ?
A LYS 348  A LYS 348  STRN         A PRO 349  A PRO 349  STRN32         ? ?
A ILE 352  A ILE 352  STRN         A ASN 356  A ASN 356  STRN33         ? ?
A HIS 361  A HIS 361  STRN         A HIS 361  A HIS 361  STRN34         ? ?
A PRO 362  A PRO 362  TURN_TY1_P   A HIS 364  A HIS 364  TURN_TY1_P19   ? ?
A GLY 365  A GLY 365  STRN         A GLY 365  A GLY 365  STRN35         ? ?
A GLN 366  A GLN 366  TURN_TY1_P   A GLN 366  A GLN 366  TURN_TY1_P20   ? ?
A GLU 370  A GLU 370  HELX_RH_AL_P A GLN 382  A GLN 382  HELX_RH_AL_P4  ? ?
A ASN 383  A ASN 383  TURN_TY1_P   A ASN 384  A ASN 384  TURN_TY1_P21   ? ?
A ALA 387  A ALA 387  STRN         A ARG 389  A ARG 389  STRN36         ? ?
A SER 391  A SER 391  TURN_TY1_P   A HIS 392  A HIS 392  TURN_TY1_P22   ? ?
A TYR 393  A TYR 393  BEND         A TYR 393  A TYR 393  BEND28         ? ?
A PRO 397  A PRO 397  HELX_RH_AL_P A TYR 406  A TYR 406  HELX_RH_AL_P5  ? ?
A GLY 407  A GLY 407  TURN_TY1_P   A GLY 407  A GLY 407  TURN_TY1_P23   ? ?
A TYR 409  A TYR 409  STRN         A GLU 413  A GLU 413  STRN37         ? ?
A ASN 415  A ASN 415  BEND         A ASN 415  A ASN 415  BEND29         ? ?
A GLU 417  A GLU 417  STRN         A GLU 417  A GLU 417  STRN38         ? ?
A HIS 419  A HIS 419  TURN_TY1_P   A GLY 420  A GLY 420  TURN_TY1_P24   ? ?
A MET 421  A MET 421  BEND         A PRO 423  A PRO 423  BEND30         ? ?
A MET 424  A MET 424  HELX_RH_3T_P A ARG 426  A ARG 426  HELX_RH_3T_P6  ? ?
A LEU 427  A LEU 427  TURN_TY1_P   A ASP 429  A ASP 429  TURN_TY1_P25   ? ?
A PRO 431  A PRO 431  TURN_TY1_P   A TRP 433  A TRP 433  TURN_TY1_P26   ? ?
A LEU 434  A LEU 434  HELX_RH_AL_P A ASP 448  A ASP 448  HELX_RH_AL_P6  ? ?
A ARG 449  A ARG 449  TURN_TY1_P   A ASN 450  A ASN 450  TURN_TY1_P27   ? ?
A PRO 452  A PRO 452  BEND         A SER 453  A SER 453  BEND31         ? ?
A VAL 454  A VAL 454  STRN         A SER 458  A SER 458  STRN39         ? ?
A ASN 461  A ASN 461  STRN         A ASN 461  A ASN 461  STRN40         ? ?
A GLU 462  A GLU 462  BEND         A GLU 462  A GLU 462  BEND32         ? ?
A ALA 467  A ALA 467  HELX_RH_AL_P A VAL 479  A VAL 479  HELX_RH_AL_P7  ? ?
A PRO 481  A PRO 481  TURN_TY1_P   A SER 482  A SER 482  TURN_TY1_P28   ? ?
A ARG 483  A ARG 483  BEND         A ARG 483  A ARG 483  BEND33         ? ?
A VAL 485  A VAL 485  STRN         A GLN 486  A GLN 486  STRN41         ? ?
A GLY 489  A GLY 489  TURN_TY1_P   A GLY 491  A GLY 491  TURN_TY1_P29   ? ?
A ALA 492  A ALA 492  BEND         A THR 494  A THR 494  BEND34         ? ?
A THR 495  A THR 495  TURN_TY1_P   A ALA 496  A ALA 496  TURN_TY1_P30   ? ?
A ASP 498  A ASP 498  BEND         A ILE 499  A ILE 499  BEND35         ? ?
A ILE 500  A ILE 500  STRN         A ILE 500  A ILE 500  STRN42         ? ?
A MET 503  A MET 503  STRN         A MET 503  A MET 503  STRN43         ? ?
A TYR 504  A TYR 504  BEND         A TYR 504  A TYR 504  BEND36         ? ?
A VAL 507  A VAL 507  BEND         A ASP 508  A ASP 508  BEND37         ? ?
A GLN 511  A GLN 511  STRN         A GLN 511  A GLN 511  STRN44         ? ?
A PRO 514  A PRO 514  BEND         A VAL 516  A VAL 516  BEND38         ? ?
A LYS 518  A LYS 518  STRN         A LYS 518  A LYS 518  STRN45         ? ?
A ILE 521  A ILE 521  HELX_RH_AL_P A TRP 524  A TRP 524  HELX_RH_AL_P8  ? ?
A LEU 525  A LEU 525  TURN_TY1_P   A SER 526  A SER 526  TURN_TY1_P31   ? ?
A LEU 527  A LEU 527  BEND         A LEU 527  A LEU 527  BEND39         ? ?
A PRO 528  A PRO 528  TURN_TY1_P   A GLY 529  A GLY 529  TURN_TY1_P32   ? ?
A ARG 532  A ARG 532  BEND         A ARG 532  A ARG 532  BEND40         ? ?
A LEU 534  A LEU 534  STRN         A TYR 539  A TYR 539  STRN46         ? ?
A ALA 540  A ALA 540  BEND         A ALA 540  A ALA 540  BEND41         ? ?
A GLY 544  A GLY 544  BEND         A ASN 545  A ASN 545  BEND42         ? ?
A GLY 548  A GLY 548  TURN_TY1_P   A GLY 549  A GLY 549  TURN_TY1_P33   ? ?
A PHE 550  A PHE 550  HELX_RH_AL_P A GLN 559  A GLN 559  HELX_RH_AL_P9  ? ?
A PRO 561  A PRO 561  TURN_TY1_P   A ARG 562  A ARG 562  TURN_TY1_P34   ? ?
A LEU 563  A LEU 563  STRN         A VAL 568  A VAL 568  STRN47         ? ?
A TRP 569  A TRP 569  BEND         A TRP 569  A TRP 569  BEND43         ? ?
A TRP 571  A TRP 571  BEND         A TRP 571  A TRP 571  BEND44         ? ?
A VAL 572  A VAL 572  STRN         A VAL 572  A VAL 572  STRN48         ? ?
A LEU 576  A LEU 576  STRN         A TYR 579  A TYR 579  STRN49         ? ?
A GLU 581  A GLU 581  TURN_TY1_P   A ASN 582  A ASN 582  TURN_TY1_P35   ? ?
A PRO 585  A PRO 585  STRN         A ALA 588  A ALA 588  STRN50         ? ?
A GLY 590  A GLY 590  TURN_TY1_P   A PHE 593  A PHE 593  TURN_TY1_P36   ? ?
A GLY 594  A GLY 594  BEND         A GLY 594  A GLY 594  BEND45         ? ?
A THR 596  A THR 596  BEND         A PRO 597  A PRO 597  BEND46         ? ?
A ARG 600  A ARG 600  HELX_RH_3T_P A CYS 603  A CYS 603  HELX_RH_3T_P7  ? ?
A LEU 607  A LEU 607  BEND         A LEU 607  A LEU 607  BEND47         ? ?
A VAL 608  A VAL 608  STRN         A VAL 608  A VAL 608  STRN51         ? ?
A ALA 610  A ALA 610  TURN_TY1_P   A ASP 611  A ASP 611  TURN_TY1_P37   ? ?
A ARG 612  A ARG 612  BEND         A ARG 612  A ARG 612  BEND48         ? ?
A THR 613  A THR 613  HELX_LH_PP_P A THR 613  A THR 613  HELX_LH_PP_P5  ? ?
A PRO 614  A PRO 614  STRN         A PRO 614  A PRO 614  STRN52         ? ?
A HIS 615  A HIS 615  HELX_LH_PP_P A HIS 615  A HIS 615  HELX_LH_PP_P6  ? ?
A PRO 616  A PRO 616  TURN_TY1_P   A PRO 616  A PRO 616  TURN_TY1_P38   ? ?
A ALA 617  A ALA 617  HELX_RH_AL_P A GLN 624  A GLN 624  HELX_RH_AL_P10 ? ?
A GLN 625  A GLN 625  TURN_TY1_P   A GLN 625  A GLN 625  TURN_TY1_P39   ? ?
A PHE 627  A PHE 627  BEND         A PHE 627  A PHE 627  BEND49         ? ?
A PHE 628  A PHE 628  STRN         A SER 633  A SER 633  STRN53         ? ?
A GLY 634  A GLY 634  TURN_TY1_P   A GLN 635  A GLN 635  TURN_TY1_P40   ? ?
A THR 636  A THR 636  STRN         A SER 641  A SER 641  STRN54         ? ?
A LEU 644  A LEU 644  BEND         A PHE 645  A PHE 645  BEND50         ? ?
A SER 648  A SER 648  BEND         A SER 648  A SER 648  BEND51         ? ?
A ASP 649  A ASP 649  TURN_TY1_P   A ASN 650  A ASN 650  TURN_TY1_P41   ? ?
A LEU 652  A LEU 652  STRN         A LEU 659  A LEU 659  STRN55         ? ?
A ASP 660  A ASP 660  TURN_TY1_P   A GLY 661  A GLY 661  TURN_TY1_P42   ? ?
A LYS 662  A LYS 662  STRN         A PRO 670  A PRO 670  STRN56         ? ?
A PRO 675  A PRO 675  TURN_TY1_P   A GLN 676  A GLN 676  TURN_TY1_P43   ? ?
A LYS 678  A LYS 678  STRN         A GLU 682  A GLU 682  STRN57         ? ?
A LEU 683  A LEU 683  HELX_LH_PP_P A PRO 687  A PRO 687  HELX_LH_PP_P7  ? ?
A GLU 690  A GLU 690  BEND         A SER 691  A SER 691  BEND52         ? ?
A GLY 693  A GLY 693  STRN         A GLN 703  A GLN 703  STRN58         ? ?
A PRO 704  A PRO 704  BEND         A ASN 705  A ASN 705  BEND53         ? ?
A ALA 708  A ALA 708  BEND         A SER 710  A SER 710  BEND54         ? ?
A ALA 712  A ALA 712  TURN_TY1_P   A GLY 713  A GLY 713  TURN_TY1_P44   ? ?
A HIS 714  A HIS 714  STRN         A ASN 726  A ASN 726  STRN59         ? ?
A THR 730  A THR 730  HELX_LH_PP_P A SER 735  A SER 735  HELX_LH_PP_P8  ? ?
A HIS 736  A HIS 736  BEND         A ALA 737  A ALA 737  BEND55         ? ?
A PRO 739  A PRO 739  HELX_LH_PP_P A PRO 739  A PRO 739  HELX_LH_PP_P9  ? ?
A HIS 740  A HIS 740  STRN         A THR 743  A THR 743  STRN60         ? ?
A GLU 745  A GLU 745  BEND         A MET 746  A MET 746  BEND56         ? ?
A ASP 747  A ASP 747  STRN         A LEU 752  A LEU 752  STRN61         ? ?
A GLY 753  A GLY 753  TURN_TY1_P   A ASN 754  A ASN 754  TURN_TY1_P45   ? ?
A LYS 755  A LYS 755  STRN         A ASN 760  A ASN 760  STRN62         ? ?
A ARG 761  A ARG 761  TURN_TY1_P   A SER 763  A SER 763  TURN_TY1_P46   ? ?
A PHE 765  A PHE 765  BEND         A PHE 765  A PHE 765  BEND57         ? ?
A LEU 766  A LEU 766  STRN         A ILE 771  A ILE 771  STRN63         ? ?
A GLY 772  A GLY 772  TURN_TY1_P   A ASP 773  A ASP 773  TURN_TY1_P47   ? ?
A LYS 774  A LYS 774  STRN         A LYS 775  A LYS 775  STRN64         ? ?
A LEU 777  A LEU 777  STRN         A GLN 784  A GLN 784  STRN65         ? ?
A ALA 788  A ALA 788  HELX_LH_PP_P A LEU 790  A LEU 790  HELX_LH_PP_P10 ? ?
A ASP 791  A ASP 791  HELX_RH_AL_P A ILE 794  A ILE 794  HELX_RH_AL_P11 ? ?
A GLY 795  A GLY 795  TURN_TY1_P   A VAL 796  A VAL 796  TURN_TY1_P48   ? ?
A ALA 799  A ALA 799  BEND         A ARG 801  A ARG 801  BEND58         ? ?
A PRO 804  A PRO 804  TURN_TY1_P   A ASN 805  A ASN 805  TURN_TY1_P49   ? ?
A ALA 806  A ALA 806  BEND         A ALA 806  A ALA 806  BEND59         ? ?
A TRP 807  A TRP 807  HELX_RH_AL_P A ALA 813  A ALA 813  HELX_RH_AL_P12 ? ?
A ALA 814  A ALA 814  TURN_TY1_P   A GLY 815  A GLY 815  TURN_TY1_P50   ? ?
A HIS 816  A HIS 816  HELX_RH_3T_P A GLN 818  A GLN 818  HELX_RH_3T_P8  ? ?
A GLU 820  A GLU 820  STRN         A THR 830  A THR 830  STRN66         ? ?
A ALA 832  A ALA 832  BEND         A ASP 833  A ASP 833  BEND60         ? ?
A ALA 834  A ALA 834  STRN         A HIS 845  A HIS 845  STRN67         ? ?
A GLN 846  A GLN 846  TURN_TY1_P   A GLY 847  A GLY 847  TURN_TY1_P51   ? ?
A LYS 848  A LYS 848  STRN         A ASP 860  A ASP 860  STRN68         ? ?
A GLY 861  A GLY 861  TURN_TY1_P   A SER 862  A SER 862  TURN_TY1_P52   ? ?
A GLY 863  A GLY 863  BEND         A GLY 863  A GLY 863  BEND61         ? ?
A MET 865  A MET 865  STRN         A VAL 873  A VAL 873  STRN69         ? ?
A SER 875  A SER 875  TURN_TY1_P   A ASP 876  A ASP 876  TURN_TY1_P53   ? ?
A THR 877  A THR 877  BEND         A THR 877  A THR 877  BEND62         ? ?
A PRO 878  A PRO 878  HELX_LH_PP_P A PRO 880  A PRO 880  HELX_LH_PP_P11 ? ?
A ALA 881  A ALA 881  BEND         A ALA 881  A ALA 881  BEND63         ? ?
A ARG 882  A ARG 882  STRN         A LEU 889  A LEU 889  STRN70         ? ?
A ALA 890  A ALA 890  BEND         A ALA 890  A ALA 890  BEND64         ? ?
A GLU 894  A GLU 894  BEND         A GLU 894  A GLU 894  BEND65         ? ?
A ARG 895  A ARG 895  STRN         A GLY 902  A GLY 902  STRN71         ? ?
A PRO 903  A PRO 903  BEND         A GLN 904  A GLN 904  BEND66         ? ?
A TYR 907  A TYR 907  STRN         A TYR 907  A TYR 907  STRN72         ? ?
A PRO 908  A PRO 908  TURN_TY1_P   A ASP 909  A ASP 909  TURN_TY1_P54   ? ?
A ARG 910  A ARG 910  STRN         A ARG 910  A ARG 910  STRN73         ? ?
A THR 912  A THR 912  BEND         A ALA 913  A ALA 913  BEND67         ? ?
A CYS 915  A CYS 915  STRN         A PRO 922  A PRO 922  STRN74         ? ?
A LEU 923  A LEU 923  HELX_RH_3T_P A MET 926  A MET 926  HELX_RH_3T_P9  ? ?
A TYR 927  A TYR 927  BEND         A TYR 927  A TYR 927  BEND68         ? ?
A VAL 931  A VAL 931  BEND         A PHE 932  A PHE 932  BEND69         ? ?
A GLU 935  A GLU 935  BEND         A GLU 935  A GLU 935  BEND70         ? ?
A ARG 939  A ARG 939  STRN         A TYR 947  A TYR 947  STRN75         ? ?
A GLY 948  A GLY 948  TURN_TY1_P   A PRO 949  A PRO 949  TURN_TY1_P55   ? ?
A HIS 950  A HIS 950  STRN         A SER 961  A SER 961  STRN76         ? ?
A ARG 962  A ARG 962  BEND         A TYR 963  A TYR 963  BEND71         ? ?
A GLN 965  A GLN 965  HELX_RH_AL_P A GLU 970  A GLU 970  HELX_RH_AL_P13 ? ?
A SER 972  A SER 972  BEND         A HIS 973  A HIS 973  BEND72         ? ?
A ARG 974  A ARG 974  HELX_RH_3T_P A LEU 976  A LEU 976  HELX_RH_3T_P10 ? ?
A ALA 979  A ALA 979  HELX_LH_PP_P A GLU 980  A GLU 980  HELX_LH_PP_P12 ? ?
A GLU 981  A GLU 981  BEND         A GLY 982  A GLY 982  BEND73         ? ?
A THR 983  A THR 983  STRN         A HIS 991  A HIS 991  STRN77         ? ?
A ASP 997  A ASP 997  BEND         A ASP 997  A ASP 997  BEND74         ? ?
A SER 999  A SER 999  BEND         A SER 1001 A SER 1001 BEND75         ? ?
A ALA 1006 A ALA 1006 HELX_RH_3T_P A PHE 1008 A PHE 1008 HELX_RH_3T_P11 ? ?
A GLN 1009 A GLN 1009 BEND         A GLN 1009 A GLN 1009 BEND76         ? ?
A ALA 1012 A ALA 1012 BEND         A GLY 1013 A GLY 1013 BEND77         ? ?
A ARG 1014 A ARG 1014 STRN         A GLN 1023 A GLN 1023 STRN78         ? ?
#
    """

  pdb_answer = """\
HELIX    1   1 LEU A    8  ARG A   14  1                                   6
HELIX    2   2 ASP A   16  GLU A   18  5                                   2
HELIX    3   3 SER A   40  THR A   45  1                                   5
HELIX    4   4 PRO A   67  ALA A   69  5                                   2
HELIX    5   5 GLU A   72  LEU A   75  5                                   3
HELIX    6   6 TRP A   91  HIS A   94  5                                   3
HELIX    7   7 GLU A  132  GLN A  136  1                                   4
HELIX    8   8 ASP A  194  LEU A  198  5                                   4
HELIX    9   9 GLU A  370  GLN A  382  1                                  12
HELIX   10  10 PRO A  397  TYR A  406  1                                   9
HELIX   11  11 MET A  424  ARG A  426  5                                   2
HELIX   12  12 LEU A  434  ASP A  448  1                                  14
HELIX   13  13 ALA A  467  VAL A  479  1                                  12
HELIX   14  14 ILE A  521  TRP A  524  1                                   3
HELIX   15  15 PHE A  550  GLN A  559  1                                   9
HELIX   16  16 ARG A  600  CYS A  603  5                                   3
HELIX   17  17 ALA A  617  GLN A  624  1                                   7
HELIX   18  18 ASP A  791  ILE A  794  1                                   3
HELIX   19  19 TRP A  807  ALA A  813  1                                   6
HELIX   20  20 HIS A  816  GLN A  818  5                                   2
HELIX   21  21 LEU A  923  MET A  926  5                                   3
HELIX   22  22 GLN A  965  GLU A  970  1                                   5
HELIX   23  23 ARG A  974  LEU A  976  5                                   2
HELIX   24  24 ALA A 1006  PHE A 1008  5                                   2"""

  cif_model = iotbx.cif.reader(input_string=struct_conf_str).model()
  cif_block = list(cif_model.values())[0]
  ann = annotation.from_cif_block(cif_block)
  # print(ann.as_pdb_str())
  assert not show_diff(ann.as_pdb_str(), pdb_answer)

if (__name__ == "__main__"):
  t0 = time.time()
  test_helix_interface()
  test_sheet_interface()
  tst_from_cif_block()
  tst_from_cif_block_2()
  tst_from_cif_block_3()
  tst_from_cif_block_4()
  tst_from_minimal_cif_helix()
  tst_from_minimal_cif_sheet()
  tst_to_cif_helix()
  tst_to_cif_sheet()
  tst_to_cif_annotation()
  tst_remove_empty_annotations()
  tst_split_helices_with_prolines()
  tst_split_helices_with_prolines_2()
  tst_split_helices_with_prolines_3()
  tst_remove_short_annotations()
  tst_remove_3_10_helices()
  tst_concatenate_consecutive_helices()
  tst_concatenate_consecutive_helices2()
  tst_concatenate_consecutive_helices3()
  tst_concatenate_consecutive_helices4()
  tst_filter_sheets_with_long_hbonds()
  tst_filter_sheets_with_long_hbonds2()
  tst_filter_sheets_with_long_hbonds3()
  tst_reset_sheet_ids()
  tst_simple_elements()
  tst_struct_conf_various()
  print("OK time =%8.3f"%(time.time() - t0))
