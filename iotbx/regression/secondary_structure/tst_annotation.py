from __future__ import division
import time
from iotbx.pdb.secondary_structure import annotation, pdb_helix

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

if (__name__ == "__main__"):
  t0 = time.time()
  test_helix_interface()
  test_sheet_interface()
  print "OK time =%8.3f"%(time.time() - t0)
