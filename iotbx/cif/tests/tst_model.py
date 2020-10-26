from __future__ import absolute_import, division, print_function
from cctbx.array_family import flex
from libtbx.test_utils import Exception_expected, show_diff
from libtbx.utils import Sorry
from libtbx.containers import OrderedDict
from six.moves import cStringIO as StringIO
import copy

import iotbx.cif
from iotbx.cif import model

def exercise_cif_model():
  loop = model.loop()
  loop["_loop_a"] = flex.double((1,2,3))
  loop.add_columns({'_loop_c': [4,5,6],
                    '_loop_b': ['7','8','9']})
  loop.add_row((4,7,'0'))
  try: loop["_loop_invalid"] = (7,8)
  except AssertionError: pass
  else: raise Exception_expected
  assert len(loop) == 3 # the number of columns (keys)
  assert loop.size() == 4 # the number of rows (loop iterations)
  assert list(loop.keys()) == ['_loop_a', '_loop_c', '_loop_b']
  try: loop["no_leading_underscore"] = 3
  except Sorry: pass
  else: raise Exception_expected
  loop2 = model.loop(header=("_loop2_a", "_loop2_b"), data=(1,2,3,4,5,6))
  assert list(loop2.keys()) == ["_loop2_a", "_loop2_b"]
  assert list(loop2.values()) == [flex.std_string(['1', '3', '5']),
                            flex.std_string(['2', '4', '6'])]
  assert list(loop2.iterrows()) == [
    {'_loop2_a': '1', '_loop2_b': '2'},
    {'_loop2_a': '3', '_loop2_b': '4'},
    {'_loop2_a': '5', '_loop2_b': '6'}]
  loop3 = model.loop(
    data={
      "_loop3_a": flex.int((-1, 2, 3)),
      "_loop3_b": flex.double((1.1, 2.2, 3.3)),
      "_loop3_c": flex.size_t((1, 2, 3)),
    }
  )
  for k in "abc":
    assert isinstance(loop3["_loop3_%s" % k], flex.std_string)
  #
  block = model.block()
  block["_tag"] = 3
  block["_tag1"] = "'a string'"
  block["_another_tag"] = 3.142
  assert "_tag" in block
  assert "_tag1" in block
  assert "_another_tag" in block
  assert block["_tag"] == '3'
  assert block["_tag1"] == "'a string'"
  assert block["_another_tag"] == "3.142"
  assert list(block.keys()) == ['_tag', '_tag1', '_another_tag']
  assert list(block.values()) == ["3", "'a string'", "3.142"]
  try: block["no_leading_underscore"] = 3
  except Sorry: pass
  else: raise Exception_expected
  block.add_loop(loop)
  assert len(block) == 6
  assert list(block.items()) == [
    ('_tag', '3'), ('_tag1', "'a string'"), ('_another_tag', '3.142'),
    ('_loop_a', flex.std_string(['1', '2', '3', '4'])),
    ('_loop_c', flex.std_string(['4', '5', '6', '7'])),
    ('_loop_b', flex.std_string(['7', '8', '9', '0']))]
  block['_loop_c'] = [11, 12, 13, 14]
  assert '_loop_c' in list(block.loops['_loop'].keys())
  assert list(block['_loop_c']) == ['11', '12', '13', '14']
  #
  block1 = model.block()
  block1["_tag"] = 2
  block1["_tag2"] = 1.2
  loop3 = model.loop(header=("_loop_a", "_loop_b"), data=(6,5,4,3,2,1))
  block1.add_loop(loop2)
  block1.add_loop(loop3)
  block.update(block1)
  for key in block._items.keys():
    assert key in ['_another_tag', '_tag2', '_tag', '_tag1']
  for value in block._items.values():
    assert value in ['3.142', '1.2', '2', "'a string'"]
  assert list(block.loops.keys()) == ['_loop', '_loop2']
  assert list(block.keys()) == ['_tag', '_tag1', '_another_tag', '_loop_a',
                          '_loop_b','_tag2', '_loop2_a', '_loop2_b']
  assert list(block['_loop_a']) == ['6', '4', '2']
  assert list(block['_loop_b']) == ['5', '3', '1']
  assert list(block['_loop2_a']) == ['1', '3', '5']
  assert list(block['_loop2_b']) == ['2', '4', '6']
  try: block['_loop_c']
  except KeyError: pass
  else: raise Exception_expected
  bad_loop = model.loop(header=("_a", "_b"), data=(1,2,3,4,5,6))
  block1.add_loop(bad_loop)
  assert "_a" in block1
  assert "_b" in block1
  assert list(block.get_looped_item("_loop_a")) == ['6', '4', '2']
  try: block.get_looped_item("_tag", value_error=ValueError)
  except ValueError: pass
  else: raise Exception_expected
  assert list(block.get_looped_item("_tag", value_error=None)) == ['2']
  try: block.get_looped_item("_none_existent")
  except KeyError: pass
  else: raise Exception_expected
  assert block.get_looped_item(
    "_none_existent", key_error=None, default="my_default") == "my_default"
  assert block.get_single_item("_tag") == "2"
  try: block.get_single_item("_loop_a")
  except ValueError: pass
  else: raise Exception_expected
  assert block.get_single_item(
    "_loop_a", value_error=None, default="default") == "default"
  try: block.get_single_item("_none_existent")
  except KeyError: pass
  else: raise Exception_expected
  assert block.get_single_item("_none_existent", key_error=None) is None
  #
  cif_model = model.cif()
  cif_model["fred"] = block
  assert "fred" in cif_model
  assert cif_model["frEd"] is block
  assert cif_model["fred"]["_Tag"] == '2'
  cif_model["fred"]["_tag"] = 4
  assert cif_model["fred"]["_tag"] == '4'
  del cif_model["fred"]["_tAg"]
  try: cif_model["fred"]["_tag"]
  except KeyError: pass
  else: raise Exception_expected
  cm = cif_model.deepcopy()
  l = cm["fred"]["_loop"]
  del cm["Fred"]["_loop_B"]
  assert "_loop_b" not in cm["fred"]
  assert "_loop_b" not in l
  assert "_loop" in cm["fred"].loops
  del cm["fred"]["_loop_a"]
  assert "_loop" not in cm["fred"].loops
  del cm["fred"]["_loop2"]
  assert "_loop2" not in cm["fred"].loops
  s = StringIO()
  print(cm, file=s)
  assert not show_diff(s.getvalue(),
"""\
data_fred
_tag1                             'a string'
_another_tag                      3.142
_tag2                             1.2

""")
  #
  cm2 = cif_model.copy()
  cm3 = cif_model.deepcopy()
  assert cm2['fred']['_loop_a'] is cif_model ['fred']['_loop_a']
  assert cm3['fred']['_loop_a'] is not cif_model ['fred']['_loop_a']
  b2 = copy.copy(block)
  b3 = copy.deepcopy(block)
  assert b2['_loop_b'] is block['_loop_b']
  assert b3['_loop_b'] is not block['_loop_b']
  l2 = loop.copy()
  l3 = loop.deepcopy()
  assert l2['_loop_b'] is loop['_loop_b']
  assert l3['_loop_b'] is not loop['_loop_b']
  #
  s = StringIO()
  cif_model.show(out=s)
  assert not show_diff(s.getvalue(),
"""\
data_fred
_tag1                             'a string'
_another_tag                      3.142
loop_
  _loop_a
  _loop_b
  6  5
  4  3
  2  1

_tag2                             1.2
loop_
  _loop2_a
  _loop2_b
  1  2
  3  4
  5  6

""")
  s = StringIO()
  cif_model.show(out=s, indent="    ", data_name_field_width=0)
  assert not show_diff(s.getvalue(),
"""\
data_fred
_tag1 'a string'
_another_tag 3.142
loop_
    _loop_a
    _loop_b
    6  5
    4  3
    2  1

_tag2 1.2
loop_
    _loop2_a
    _loop2_b
    1  2
    3  4
    5  6

""")
  s = StringIO()
  cif_model.show(out=s, indent="", indent_row="   ", data_name_field_width=0)
  assert not show_diff(s.getvalue(),
"""\
data_fred
_tag1 'a string'
_another_tag 3.142
loop_
_loop_a
_loop_b
   6  5
   4  3
   2  1

_tag2 1.2
loop_
_loop2_a
_loop2_b
   1  2
   3  4
   5  6

""")
  cif_model.sort(recursive=True)
  s = StringIO()
  cif_model.show(out=s)
  assert not show_diff(s.getvalue(),
"""\
data_fred
_another_tag                      3.142
_tag1                             'a string'
_tag2                             1.2
loop_
  _loop_a
  _loop_b
  6  5
  4  3
  2  1

loop_
  _loop2_a
  _loop2_b
  1  2
  3  4
  5  6

""")
  save = model.save()
  save.add_loop(l3)
  save['_tag1'] = 3
  block = model.block()
  block['bob'] = save
  cm = model.cif({'fred': block})
  s = StringIO()
  cm.show(out=s)
  assert not show_diff(s.getvalue(),
"""data_fred

save_bob
   loop_
    _loop_a
    _loop_c
    _loop_b
    1  11  7
    2  12  8
    3  13  9
    4  14  0

  _tag1                             3
  save_

""")
  b1 = model.block()
  b1['_a'] = 1
  b1['_b'] = 2
  b1['_c'] = 3
  b2 = model.block()
  b2['_a'] = 2
  b2['_c'] = 3
  b2['_d'] = 4
  b3 = b1.difference(b2)
  b4 = b2.difference(b1)
  for item in b3.items():
    assert item in [('_b', '2'), ('_a', '2')]
  for item in b4.items():
    assert item in [('_d', '4'), ('_a', '1')]
  l = model.loop(data=dict(_loop_d=(1,2),_loop_e=(3,4),_loop_f=(5,6)))
  assert l == l
  assert l == l.deepcopy()
  assert l != l2
  assert l != l3
  l2 = model.loop(data=dict(_loop_d=(1,2,3),_loop_e=(3,4,5),_loop_f=(5,6,7)))
  b1.add_loop(l)
  b2.add_loop(l2)
  b5 = b1.difference(b2)
  assert b5['_loop'] == l2
  l = model.loop(data=OrderedDict((('_loop_a',(1,21,-13)),
                                   ('_loop_b',(-221.3,3.01,4.246)),
                                   ('_loop_c',("a","b","cc")))))
  b = model.block()
  b.add_loop(l)
  cm = model.cif({'fred':b})
  s = StringIO()
  cm.show(out=s, loop_format_strings={'_loop':'% 4i% 8.2f %s'})
  assert not show_diff(s.getvalue(),"""\
data_fred
loop_
  _loop_a
  _loop_b
  _loop_c
   1 -221.30 a
  21    3.01 b
 -13    4.25 cc

""")
  s = StringIO()
  cm.show(out=s)
  assert not show_diff(s.getvalue(),"""\
data_fred
loop_
  _loop_a
  _loop_b
  _loop_c
    1  -221.3  a
   21    3.01  b
  -13   4.246  cc

""")
  l.add_row((".", "?", "."))
  s = StringIO()
  cm.show(out=s)
  assert not show_diff(s.getvalue(),"""\
data_fred
loop_
  _loop_a
  _loop_b
  _loop_c
    1  -221.3  a
   21    3.01  b
  -13   4.246  cc
    .       ?  .

""")
  l.delete_row(index=1)
  s = StringIO()
  cm.show(out=s)
  assert not show_diff(s.getvalue(),"""\
data_fred
loop_
  _loop_a
  _loop_b
  _loop_c
    1  -221.3  a
  -13   4.246  cc
    .       ?  .

""")
  l2 = l.deepcopy()
  l2.delete_row(index=0)
  l2.delete_row(index=0)
  l2.delete_row(index=0)
  try: l2.show(out=s)
  except AssertionError as e: pass
  else: raise Exception_expected
  l.clear()
  try: l.show(out=s)
  except AssertionError as e: pass
  else: raise Exception_expected
  #
  loop = model.loop(data={"_a_1": ('string with spaces','nospaces'),
                          "_a_2": ('a', 'b')})
  s = StringIO()
  loop.show(out=s, align_columns=True)
  assert not show_diff(s.getvalue(), """\
loop_
  _a_1
  _a_2
  'string with spaces'  a
  nospaces              b
""")
  #
  cb = model.block()
  cm = model.cif()
  cm["a"] = cb
  cb["_b"] = ""
  s = StringIO()
  cm.show(out=s)
  assert not show_diff(s.getvalue(), """\
data_a
_b                                ''
""")
  #
  loop = model.loop(data=OrderedDict((
    ("_entity_poly.entity_id", ('1', '2', '3')),
    ("_entity_poly.pdbx_seq_one_letter_code", (
      "TFGSGEADCGLRPLFEKKSLEDKTERELLESYIDGR",
      """\
IVEGSDAEIGMSPWQVMLFRKSPQELLCGASLISDRWVLTAAHCLLYPPWDKNFTENDLLVRIGKHSRTRYERNIEKISM
THVFRLKKWIQKVIDQFGE""",
      "NGDFEEIPEE(TYS)LQ",
    )),
    ("_entity_poly.pdbx_seq_one_letter_code_can", (
      "TFGSGEADCGLRPLFEKKSLEDKTERELLESYIDGR",
      """\
IVEGSDAEIGMSPWQVMLFRKSPQELLCGASLISDRWVLTAAHCLLYPPWDKNFTENDLLVRIGKHSRTRYERNIEKISM
THVFRLKKWIQKVIDQFGE""",
      "NGDFEEIPEEYLQ",
    )),
    ("_entity_poly.pdbx_strand_id", ('L', 'H', 'I'))
  )))
  s = StringIO()
  loop.show(out=s, align_columns=True)
  s.seek(0)
  assert not show_diff("\n".join(l.rstrip() for l in s.readlines()),"""\
loop_
  _entity_poly.entity_id
  _entity_poly.pdbx_seq_one_letter_code
  _entity_poly.pdbx_seq_one_letter_code_can
  _entity_poly.pdbx_strand_id
  1  TFGSGEADCGLRPLFEKKSLEDKTERELLESYIDGR  TFGSGEADCGLRPLFEKKSLEDKTERELLESYIDGR  L
  2
;
IVEGSDAEIGMSPWQVMLFRKSPQELLCGASLISDRWVLTAAHCLLYPPWDKNFTENDLLVRIGKHSRTRYERNIEKISM
THVFRLKKWIQKVIDQFGE
;

;
IVEGSDAEIGMSPWQVMLFRKSPQELLCGASLISDRWVLTAAHCLLYPPWDKNFTENDLLVRIGKHSRTRYERNIEKISM
THVFRLKKWIQKVIDQFGE
;
  H
  3  NGDFEEIPEE(TYS)LQ                     NGDFEEIPEEYLQ                         I\
""")
  #
  cb = model.block()
  cm = model.cif()
  cm["a"] = cb
  cb["_a"] = '1 "a" 2'
  cb["_b"] = "1 'b' 3"
  cb["_c"] = "O1'"
  cb["_d"] = 'O2"'
  cb["_e"] = """1 'a' "b" 3"""
  s = StringIO()
  print(cm, file=s)
  s.seek(0)
  assert not show_diff("\n".join(l.rstrip() for l in s.readlines()), """\
data_a
_a                                '1 "a" 2'
_b                                "1 'b' 3"
_c                                O1'
_d                                O2"
_e
;
1 'a' "b" 3
;

""")
  # verify that what we wrote out above is valid CIF and we can read it back in
  cm2 = iotbx.cif.reader(input_string=s.getvalue()).model()
  cb2 = cm2["a"]
  assert cb2["_a"] == cb["_a"]
  assert cb2["_b"] == cb["_b"]
  assert cb2["_c"] == cb["_c"]
  assert cb2["_d"] == cb["_d"]
  assert cb2["_e"].strip() == cb["_e"]
  #
  cm = iotbx.cif.reader(input_string="""\
data_a
loop_
  _pdbx_refine_tls_group.id
  _pdbx_refine_tls_group.refine_tls_id
  _pdbx_refine_tls_group.selection
  _pdbx_refine_tls_group.selection_details
  1  1  ?  "chain 'A' and (resid    2  through   15 )"
  2  2  ?  "chain 'A' and (resid   16  through   26 )"
  3  3  ?  "chain 'A' and (resid   27  through   43 )"
  4  4  ?  "chain 'B' and (resid    1  through   14 )"
  5  5  ?  "chain 'B' and (resid   15  through   20 )"
""").model()
  print(cm)
  #
  cif_block = model.block()
  loop_a = model.loop(header=("_a.1", "_a.2"), data=(1,2,3,4,5,6))
  cif_block.add_loop(loop_a)
  assert cif_block.get_loop("_a") is loop_a
  assert cif_block.get_loop_or_row("_a") is loop_a
  assert cif_block.get_loop("_b") is None
  assert cif_block.get_loop_or_row("_b") is None
  assert cif_block.get_loop("_b", default=loop_a) is loop_a
  assert cif_block.get_loop_or_row("_b", default=loop_a) is loop_a
  loop_a = cif_block.get_loop_with_defaults(
    "_a", default_dict={"_a.2":".", "_a.3":"?", "_a.4":"."})
  assert list(cif_block["_a.1"]) == ['1', '3', '5']
  assert list(cif_block["_a.2"]) == ['2', '4', '6']
  assert list(cif_block["_a.3"]) == ['?', '?', '?']
  assert list(cif_block["_a.4"]) == ['.', '.', '.']
  loop_a.add_row({"_a.3":"a", "_a.4":"b"})
  loop_a.add_row({"_a.3":"c", "_a.4":"d"}, default_value=".")
  assert list(cif_block["_a.1"]) == ['1', '3', '5', '?', '.']
  assert list(cif_block["_a.2"]) == ['2', '4', '6', '?', '.']
  assert list(cif_block["_a.3"]) == ['?', '?', '?', 'a', 'c']
  assert list(cif_block["_a.4"]) == ['.', '.', '.', 'b', 'd']
  loop_B = model.loop(header=("_B.1", "_B.2", "_B.3"), data=(1,2,3,4,5,6))
  cif_block.add_loop(loop_B)
  assert cif_block.get_loop("_B") is loop_B
  assert cif_block.get_loop_or_row("_B") is loop_B
  assert cif_block.get_loop("_b") is loop_B
  assert cif_block.get_loop_or_row("_b") is loop_B
  #
  cif_block = model.block()
  cif_block['_a'] = """\
123
456"""
  s = StringIO()
  cif_block.show(out=s)
  s.seek(0)
  assert not show_diff("\n".join([l.strip() for l in s.readlines()]), """\
_a
;
123
456
;
""")


  cm = iotbx.cif.reader(input_string="""\
data_a
  _test_row.id 1
  _test_row.data2 2
  _test_row.data3 3
  _test_row.data4 44
#
loop_
_test_row_range.sheet_id
_test_row_range.id
_test_row_range.beg_label_comp_id
_test_row_range.beg_label_asym_id
A 1 SER A
A 2 MET A
#
""").model()
  #
  cif_block = list(cm.values())[0]
  loop_or_row = cif_block.get_loop_or_row('_test_row')
  assert loop_or_row.n_rows() == 1
  assert loop_or_row.n_columns() == 4
  assert list(loop_or_row['_test_row.id']) == ['1']
  assert list(loop_or_row['_test_row.data2']) == ['2']
  assert list(loop_or_row['_test_row.data3']) == ['3']
  assert list(loop_or_row['_test_row.data4']) == ['44']
  for r in loop_or_row.iterrows():
    assert list(r['_test_row.id']) == ['1']
    assert list(r['_test_row.data2']) == ['2']
    assert list(r['_test_row.data3']) == ['3']
    assert list(r['_test_row.data4']) == ['4','4']

def test_301():
  cif_model = iotbx.cif.reader(input_string="""\
data_test
loop_
  _symmetry_symop_operation_xyz
 x,y,z
 -x,y+1/2,-z
""").model()
  del cif_model["test"]["_symmetry_symop_operation_xyz"]
  s = str(cif_model)
  assert s.strip() == 'data_test'

def test_show_not_modify():
  """
  Test that content of the object is not changed during .show()
  """
  cif_str = """\
data_r3p4rsf
#
_symmetry.entry_id               3p4r
_symmetry.space_group_name_H-M   'P 21 21 21'
_symmetry.Int_Tables_number      19
#
loop_
_symmetry_equiv.id
_symmetry_equiv.pos_as_xyz
1 'X,  Y,  Z'
2 '-X+1/2,  -Y,  Z+1/2'
3 'X+1/2,  -Y+1/2,  -Z'
4 '-X,  Y+1/2,  -Z+1/2'
"""
  cif_reader = iotbx.cif.reader(input_string=cif_str)
  model = cif_reader.model()
  print (model['r3p4rsf']['_symmetry.space_group_name_H-M'])
  print (list(model['r3p4rsf']['_symmetry_equiv.pos_as_xyz']))
  assert list(model['r3p4rsf']['_symmetry_equiv.pos_as_xyz']) ==\
      ['X,  Y,  Z', '-X+1/2,  -Y,  Z+1/2', 'X+1/2,  -Y+1/2,  -Z', '-X,  Y+1/2,  -Z+1/2']
  s = StringIO()
  model.show(out=s)
  print (model['r3p4rsf']['_symmetry.space_group_name_H-M'])
  print (list(model['r3p4rsf']['_symmetry_equiv.pos_as_xyz']))
  assert list(model['r3p4rsf']['_symmetry_equiv.pos_as_xyz']) ==\
      ['X,  Y,  Z', '-X+1/2,  -Y,  Z+1/2', 'X+1/2,  -Y+1/2,  -Z', '-X,  Y+1/2,  -Z+1/2']

if __name__ == '__main__':
  exercise_cif_model()
  test_301()
  test_show_not_modify()
  print("OK")
