from cctbx.array_family import flex
from libtbx.test_utils import Exception_expected, show_diff
from cStringIO import StringIO

def exercise_cif_model():
  import iotbx.cif
  from iotbx.cif import model
  cif_model = model.cif()
  #
  loop = model.loop()
  loop["_loop_a"] = (1,2,3)
  loop.add_columns({'_loop_c': [4,5,6],
                    '_loop_b': ['7','8','9']})
  loop.add_row((4,7,'0'))
  try: loop["_loop_invalid"] = (7,8)
  except AssertionError: pass
  else: raise Exception_expected
  assert len(loop) == 3 # the number of columns (keys)
  assert loop.size() == 4 # the number of rows (loop iterations)
  assert loop.keys() == ['_loop_a', '_loop_c', '_loop_b']
  loop2 = model.loop(header=("_loop2_a", "_loop2_b"), data=(1,2,3,4,5,6))
  assert loop2.keys() == ["_loop2_a", "_loop2_b"]
  assert loop2.values() == [flex.std_string(['1', '3', '5']),
                            flex.std_string(['2', '4', '6'])]
  assert list(loop2.iterrows()) == [['1', '2'], ['3', '4'], ['5', '6']]
  #
  block = model.block()
  block["_tag"] = 3
  block["_tag1"] = "a string"
  block["_another_tag"] = 3.142
  assert "_tag" in block
  assert "_tag1" in block
  assert "_another_tag" in block
  assert block["_tag"] == '3'
  assert block["_tag1"] == "a string"
  assert block["_another_tag"] == "3.142"
  assert block.keys() == ['_tag', '_tag1', '_another_tag']
  assert block.values() == ["3", 'a string', "3.142"]
  block.add_loop(loop)
  assert len(block) == 6
  assert block.items() == [
    ('_tag', '3'), ('_tag1', 'a string'), ('_another_tag', '3.142'),
    ('_loop_a', flex.std_string(['1', '2', '3', '4'])),
    ('_loop_c', flex.std_string(['4', '5', '6', '7'])),
    ('_loop_b', flex.std_string(['7', '8', '9', '0']))]
  #
  block1 = model.block()
  block1["_tag"] = 2
  block1["_tag2"] = 1.2
  loop3 = model.loop(header=("_loop_a", "_loop_b"), data=(6,5,4,3,2,1))
  block1.add_loop(loop2)
  block1.add_loop(loop3)
  block.update(block1)
  assert block._items.keys() == ['_another_tag', '_tag2', '_tag', '_tag1']
  assert block._items.values() == ['3.142', '1.2', '2', 'a string']
  assert block.loops.keys() == ['_loop', '_loop2']
  assert block.keys() == ['_tag', '_tag1', '_another_tag', '_loop_a',
                          '_loop_b','_tag2', '_loop2_a', '_loop2_b']
  assert list(block['_loop_a']) == ['6', '4', '2']
  assert list(block['_loop_b']) == ['5', '3', '1']
  assert list(block['_loop2_a']) == ['1', '3', '5']
  assert list(block['_loop2_b']) == ['2', '4', '6']
  try: block['_loop_c']
  except KeyError: pass
  else: raise Exception_expected
  #
  cif_model["fred"] = block
  assert "fred" in cif_model
  assert cif_model["fred"] is block
  assert cif_model["fred"]["_tag"] == '2'
  cif_model["fred"]["_tag"] = 4
  assert cif_model["fred"]["_tag"] == '4'
  del cif_model["fred"]["_tag"]
  try: cif_model["fred"]["_tag"]
  except KeyError: pass
  else: raise Exception_expected
  cm = cif_model.deepcopy()
  l = cm["fred"]["_loop"]
  del cm["fred"]["_loop_b"]
  assert not cm["fred"].has_key("_loop_b")
  assert not l.has_key("_loop_b")
  assert cm["fred"].loops.has_key("_loop")
  del cm["fred"]["_loop_a"]
  assert not cm["fred"].loops.has_key("_loop")
  del cm["fred"]["_loop2"]
  assert not cm["fred"].loops.has_key("_loop2")
  s = StringIO()
  print >> s, cm
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
  b2 = block.copy()
  b3 = block.deepcopy()
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
   6 5
   4 3
   2 1

_tag2                             1.2
loop_
  _loop2_a
  _loop2_b
   1 2
   3 4
   5 6

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
     6 5
     4 3
     2 1

_tag2 1.2
loop_
    _loop2_a
    _loop2_b
     1 2
     3 4
     5 6

""")

if __name__ == '__main__':
  exercise_cif_model()
  print "OK"
