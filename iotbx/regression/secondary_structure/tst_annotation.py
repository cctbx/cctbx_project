from __future__ import division
import time
from iotbx.pdb.secondary_structure import annotation, pdb_helix
import iotbx.cif
from libtbx.test_utils import show_diff

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
  assert h.serial == 1
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
  assert st.strand_id == 2, st.strand_id
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

def tst_from_cif_block():
  test_str_1 = """\
data_4ZTE
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
  cif_block = cif_model.values()[0]
  ann = annotation.from_cif_block(cif_block)
  assert len(ann.helices) == 4
  resnames = [x.start_resname for x in ann.helices]
  assert resnames == ["SER","LEU","SER","LEU"]
  resnames = [x.end_resname for x in ann.helices]
  assert resnames == ["LYS","GLU","LYS","GLU"]
  assert len(ann.sheets) == 10
  assert [len(x.strands) for x in ann.sheets] == [4, 3, 2, 4, 3, 4, 3, 2, 4, 3]
  # print ann.as_pdb_str()
  assert not show_diff(pdb_str, ann.as_pdb_str())

if (__name__ == "__main__"):
  t0 = time.time()
  test_helix_interface()
  test_sheet_interface()
  tst_from_cif_block()
  print "OK time =%8.3f"%(time.time() - t0)
