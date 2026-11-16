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

def read_back(text):
  return iotbx.cif.reader(input_string="data_x\n" + text).model()["x"]

def exercise_loop_show_quoting():
  # one column per quoting rule of format_value, aligned and unaligned; the
  # printed text is pinned and must read back
  loop = model.loop(data=OrderedDict((
    ("_t.number", ("1.5", "-2", "1(2)", "nan", ".5", "1e5")),
    ("_t.sparse", ("1", "?", ".", "-0", "?", "22")),
    ("_t.unknown", ("?",) * 6),
    ("_t.text", ("abc", "o", "f", "'quoted'", "a,b", "x")),
    ("_t.leader", ("", "_under", "#hash", "$dollar", "[br", "]br")),
    ("_t.blank", ("a b", "tab\tx", "it's x", "'q r'", "x y z", "q")),
  )))
  def check_read_back(text):
    block = read_back(text)
    assert list(block["_t.number"]) == ["1.5", "-2", "1(2)", "nan", ".5", "1e5"]
    assert list(block["_t.sparse"]) == ["1", "?", ".", "-0", "?", "22"]
    assert list(block["_t.text"]) == ["abc", "o", "f", "quoted", "a,b", "x"]
    assert list(block["_t.leader"]) == [
      "", "_under", "#hash", "$dollar", "[br", "]br"]
    assert list(block["_t.blank"]) == [
      "a b", "tab\tx", "it's x", "q r", "x y z", "q"]
  s = StringIO()
  loop.show(out=s)
  assert not show_diff(s.getvalue(),
    "loop_\n  _t.number\n  _t.sparse\n  _t.unknown\n  _t.text\n  _t.leader\n"
    "  _t.blank\n"
    "   1.5   1  ?  abc       ''         'a b'\n"
    "    -2   ?  ?  o         '_under'   'tab\tx'\n"
    "  1(2)   .  ?  f         '#hash'    'it's x'\n"
    "   nan  -0  ?  'quoted'  '$dollar'  'q r'\n"
    "    .5   ?  ?  a,b       '[br'      'x y z'\n"
    "   1e5  22  ?  x         ']br'      q\n")
  check_read_back(s.getvalue())
  s = StringIO()
  loop.show(out=s, align_columns=False)
  assert not show_diff(s.getvalue(),
    "loop_\n  _t.number\n  _t.sparse\n  _t.unknown\n  _t.text\n  _t.leader\n"
    "  _t.blank\n"
    "   1.5 1 ? abc '' 'a b'\n"
    "   -2 ? ? o '_under' 'tab\tx'\n"
    "   1(2) . ? f '#hash' 'it's x'\n"
    "   nan -0 ? 'quoted' '$dollar' 'q r'\n"
    "   .5 ? ? a,b '[br' 'x y z'\n"
    "   1e5 22 ? x ']br' q\n")
  check_read_back(s.getvalue())
  # text fields: a value with both kinds of quote and a blank, a value that
  # is already a text field, a value with a newline; they are left out of the
  # column width
  loop = model.loop(data=OrderedDict((
    ("_m.id", ("1", "2", "3", "4")),
    ("_m.field", ("say \"hi\" 'x'", ";\nsemi\n;", "line1\nline2", "ok")),
    ("_m.tail", ("a", "b", "c", "d")),
  )))
  s = StringIO()
  loop.show(out=s)
  assert not show_diff(s.getvalue(),
    "loop_\n  _m.id\n  _m.field\n  _m.tail\n"
    "  1  'say \"hi\" 'x''  a\n"
    "  2  \n;\nsemi\n;\n      b\n"
    "  3  \n;\nline1\nline2\n;\n  c\n"
    "  4  ok              d\n")
  block = read_back(s.getvalue())
  assert list(block["_m.id"]) == ["1", "2", "3", "4"]
  assert list(block["_m.field"]) == [
    "say \"hi\" 'x'", "\nsemi\n", "\nline1\nline2\n", "ok"]
  assert list(block["_m.tail"]) == ["a", "b", "c", "d"]

def exercise_quoting_scan():
  # the per-column scan of loop.show() flags a value exactly when
  # format_value changes it or it holds a newline (a text field, kept out
  # of the column width), and a value flex.double accepts is never changed
  probes = ["", " ", "a", "a b", " a", "a ", "a\tb", "a\nb", "\n", "a\rb",
            "a\x0bb", "a\x0cb", "'", "'a", "'a b'", "'a'", '"a"', '"a b"',
            "'a b' x", "'a'b", "'ab\"", "'a\x01b", "'a\x01 b", "'\xe9'",
            "'\xe9 '", "#a", "$a", "[a", "]a", "_a", ";a", ";a;", ";\na\n;",
            " ;a;", "\xa0;a;", "\xa0;a", "\xa0a", "a;b;", "a'b", 'a"b',
            "a#b", "a$b", "a_b", "a;b", "a[b]", "?", ".", "..", "1", "-1",
            "+1", "1.5", ".5", "1.", "1e5", "1E-5", "0x10", "nan", "inf",
            "-inf", "NaN", "Infinity", "-0", "1,2", "1_000", "abc", "abc's",
            "o", "f", "x", "\xe9", "x\xa0y", "\u2003", "a\x1c",
            # newline values format_value leaves alone
            "\n;a;\n", "'a'\nb"]
  for value in probes:
    column = flex.std_string([value])
    changed = model.format_value(value) != value
    flagged = model._column_needs_format_value(column)
    assert flagged == (changed or "\n" in value), repr(value)
    try:
      flex.double(column)
    except ValueError:
      continue
    assert not changed, repr(value)
    assert model._is_numeric_column(column), repr(value)
  assert not model._column_needs_format_value(flex.std_string())
  # the flagged value at any position of a column
  for value in ("", "a b", "#a", "a\nb", "\xa0;a;"):
    for i in range(5):
      column = flex.std_string(5, "x")
      column[i] = value
      assert model._column_needs_format_value(column), (repr(value), i)
  assert not model._column_needs_format_value(flex.std_string(5, "x"))
  # a pre-formatted text field does not widen its column
  loop = model.loop(data=OrderedDict((
    ("_p.id", ("1", "2")),
    ("_p.text", ("\n;a long preformatted text field\n;\n", "x")),
    ("_p.tail", ("A", "B")),
  )))
  s = StringIO()
  loop.show(out=s)
  assert not show_diff(s.getvalue(),
    "loop_\n  _p.id\n  _p.text\n  _p.tail\n"
    "  1  \n;a long preformatted text field\n;\n  A\n"
    "  2  x  B\n")

def exercise_loop_show_many_rows():
  # more rows than show() prints at a time
  n = 70000
  loop = model.loop(data=OrderedDict((
    ("_c.i", flex.int_range(n)),
    ("_c.s", flex.std_string(["s%d" % (i % 7) for i in range(n)])),
  )))
  expected = ["loop_", "  _c.i", "  _c.s"]
  for i in range(n):
    expected.append("  %5d  s%d" % (i, i % 7))
  s = StringIO()
  loop.show(out=s)
  assert s.getvalue() == "\n".join(expected) + "\n"

def exercise_format_rows_python():
  cols = [flex.std_string(["1", "22", "?"]),
          flex.std_string(["a", "bb", "ccc"])]
  out = model._format_rows_python(cols, 0, 3, "  ", [2, -3])
  assert out == "   1  a\n  22  bb\n   ?  ccc\n", repr(out)
  out = model._format_rows_python(cols, 1, 3, "  ", [])
  assert out == "  22 bb\n  ? ccc\n", repr(out)

class python_references(object):
  """While active, loop.show() runs through the retained Python references
  of its C++ helpers."""
  references = (("_format_rows", "_format_rows_python"),
                ("_is_numeric_column", "_flex_double_accepts"))
  def __enter__(self):
    self.saved = [(n, getattr(model, n)) for n, r in self.references]
    for n, r in self.references:
      setattr(model, n, getattr(model, r))
  def __exit__(self, *args):
    for n, f in self.saved:
      setattr(model, n, f)

def exercise_format_rows_parity():
  from iotbx.cif import ext
  cases = {
    "numeric":   ([["1.5", "22", ".", "?"], ["a", "bb", "ccc", "d"]], [4, -3]),
    "trailing":  ([["x", "y"], ["a  ", "b"]], [-1, -3]),
    "textfield": ([["1", "2"], ["\n;long text\n;\n", "b"]], [1, -1]),
    "utf8":      ([["\xe9", "ab"], ["1", "2"]], [-2, 1]),
    "unaligned": ([["1", "22"], ["a", "bb"]], []),
    # the strip removes what str.rstrip() removes, so an unquoted last
    # value loses a trailing Unicode space or separator: parity with the
    # Python path, not a feature
    "pyspace":   ([["x", "y", "z"], ["a\x1c", "b\xa0", "c\u2003"]], [1, -4]),
  }
  for name, (raw, widths) in cases.items():
    cols = [flex.std_string(c) for c in raw]
    n = cols[0].size()
    for prefix in ("", "  ", "\t"):
      ref = model._format_rows_python(cols, 0, n, prefix, widths)
      got = ext.format_loop_rows(cols, 0, n, prefix, widths)
      assert got == ref, (name, prefix, repr(got), repr(ref))
  # chunks: any split gives the same bytes
  cols = [flex.std_string([str(i) for i in range(7)]),
          flex.std_string(["v%d" % i for i in range(7)])]
  w = [2, -3]
  whole = ext.format_loop_rows(cols, 0, 7, " ", w)
  for k in range(8):
    assert ext.format_loop_rows(cols, 0, k, " ", w) + \
           ext.format_loop_rows(cols, k, 7, " ", w) == whole
  assert ext.format_loop_rows(cols, 3, 3, " ", w) == ""
  for bad in ((0, 8), (5, 3)):
    try: ext.format_loop_rows(cols, bad[0], bad[1], "", [])
    except ValueError: pass
    else: raise AssertionError(bad)
  # loop.show() writes the same bytes through C++ as through the references
  s1 = StringIO(); s2 = StringIO()
  lp = model.loop(data={"_a.x": flex.std_string(["1", "2"]),
                        "_a.y": flex.std_string(["p q", "r"])})
  lp.show(out=s1)
  with python_references():
    lp.show(out=s2)
  assert s1.getvalue() == s2.getvalue()
  assert model._format_rows is not model._format_rows_python

def exercise_numeric_column_parity():
  # _is_numeric_column decides the alignment of a column; the C++ scan of
  # plain numbers plus flex.double on the rest must give the answer of
  # flex.double on the whole column
  import random
  from iotbx.cif import ext
  probes = ["1", "-1", "+1", "1.5", ".5", "1.", "1e5", "1E-5", "0x10",
            "nan", "inf", "-inf", "NaN", "Infinity", "-0", "1,2", "1_000",
            "1(2)", "1.5(12)", "1e400", "1e-400", " 1", "1 ", "", "1d5",
            "0x1p3", "nan(1)", "1e", "e5", "--1", "1.2.3", "1e+5", "+.5",
            "-.5e-3", "5.", "1(", "1)", "1()", "(2)", "1(2", "1.0e+00",
            "00012", "1e0005", "12345678901234567890", "1e308", "1e309",
            "4.9e-324", "1\t", "1\n", "Inf", "-nan", "1E", "1e+", "1.e5",
            ".e5", ".", "?", "..", "1.5+03", "1.5-03", "+", "-", "e",
            "1e99", "1e-99", "0" * 39 + "1", "1" + "0" * 39, "1" * 41,
            "1e100", "1e-100", "\xe9", "1\xe9", "abc"]
  for value in probes:
    for column in (flex.std_string([value]),
                   flex.std_string([".", "1", value, "?"]),
                   flex.std_string(["1.5", value, "nan"])):
      assert model._is_numeric_column(column) == \
        model._flex_double_accepts(column), (repr(value), list(column))
  assert model._is_numeric_column(flex.std_string())
  assert model._is_numeric_column(flex.std_string([".", "?"]))
  # the scan: index of the first value that is neither ., ? nor a plain
  # number, from a start index; the size when there is none
  column = flex.std_string(["1", ".", "nan", "2", "?", "x"])
  assert ext.first_non_plain_number(column, 0) == 2
  assert ext.first_non_plain_number(column, 3) == 5
  assert ext.first_non_plain_number(column, 6) == 6
  assert ext.first_non_plain_number(flex.std_string(["1", "2"]), 0) == 2
  # every plain number the scan accepts is accepted by flex.double, and
  # a mutation of one is judged alike by both paths
  random.seed(0)
  def digits():
    return "".join(random.choice("0123456789")
                   for i in range(random.randrange(1, 12)))
  def plain():
    s = random.choice(["", "+", "-"])
    kind = random.randrange(3)
    if kind == 0: s += digits()
    elif kind == 1: s += digits() + "." + random.choice(["", digits()])
    else: s += "." + digits()
    if random.randrange(2):
      s += random.choice("eE") + random.choice(["", "+", "-"]) \
         + "".join(random.choice("0123456789")
                   for i in range(random.randrange(1, 3)))
    return s
  for i in range(3000):
    value = plain()
    column = flex.std_string([value])
    assert ext.first_non_plain_number(column, 0) == 1, repr(value)
    flex.double(column)
    j = random.randrange(len(value) + 1)
    mutant = value[:j] + random.choice("e.+- xg()0") + value[j:]
    column = flex.std_string([value, mutant, "?"])
    assert model._is_numeric_column(column) == \
      model._flex_double_accepts(column), repr(mutant)
  # values the scan leaves to flex.double go one at a time up to
  # _numeric_single_values, then the rest of the column in one call
  sizes = []
  class counting_flex(object):
    @staticmethod
    def double(values):
      sizes.append(values.size())
      return flex.double(values)
  model.flex = counting_flex
  try:
    column = flex.std_string(5000, "nan")
    assert model._is_numeric_column(column)
    assert sizes == [1] * model._numeric_single_values + [5000], sizes
    column[4999] = "x"
    assert not model._is_numeric_column(column)
    column = flex.std_string(5000, "1.5(3)")
    column[3000] = "1.5("
    assert not model._is_numeric_column(column)
  finally:
    model.flex = flex
  assert model._is_numeric_column is not model._flex_double_accepts

if __name__ == '__main__':
  exercise_cif_model()
  test_301()
  test_show_not_modify()
  exercise_loop_show_quoting()
  exercise_quoting_scan()
  exercise_loop_show_many_rows()
  exercise_format_rows_python()
  exercise_format_rows_parity()
  exercise_numeric_column_parity()
  print("OK")
