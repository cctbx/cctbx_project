from __future__ import absolute_import, division, print_function
from scitbx.array_family import flex
from libtbx.test_utils import show_diff
from six.moves import cStringIO as StringIO
try:
  from six.moves import cPickle as pickle
except ImportError:
  import pickle
import sys

def run(args):
  assert len(args) == 0
  #
  s = flex.smart_selection()
  assert s.all_size is None
  assert s.selected_size is None
  assert s.flags is None
  assert s.indices is None
  assert s == s
  assert not (s != s)
  d = pickle.dumps(s)
  l = pickle.loads(d)
  assert l == s
  assert s.format_summary() == "None"
  #
  p = s
  s = flex.smart_selection(flags=flex.bool())
  assert s.all_size == 0
  assert s.selected_size == 0
  assert s.flags.size() == 0
  assert list(s.indices) == []
  assert s == s
  assert not (s == p)
  assert s != p
  assert pickle.loads(pickle.dumps(s)) == s
  assert s.format_summary() == "None (empty array)"
  #
  p = s
  s = flex.smart_selection(flags=flex.bool([False,True]))
  assert s.all_size == 2
  assert s.selected_size == 1
  assert s.flags.all_eq(flex.bool([False,True]))
  assert list(s.indices) == [1]
  assert s == s
  assert s != p
  assert pickle.loads(pickle.dumps(s)) == s
  assert s.format_summary() == "1 of 2"
  #
  p = s
  s = flex.smart_selection(flags=flex.bool([True]))
  assert s.all_size == 1
  assert s.selected_size == 1
  assert s.flags.all_eq(True)
  assert list(s.indices) == [0]
  assert s == s
  assert s != p
  assert pickle.loads(pickle.dumps(s)) == s
  assert s.format_summary() == "all (1)"
  #
  p = s
  s = flex.smart_selection(indices=flex.size_t())
  assert s.all_size is None
  assert s.selected_size == 0
  assert s.flags is None
  assert list(s.indices) == []
  assert s == s
  assert s != p
  assert pickle.loads(pickle.dumps(s)) == s
  assert s.format_summary() == "0"
  #
  p = s
  s = flex.smart_selection(indices=flex.size_t([2,4]))
  assert s.all_size is None
  assert s.selected_size == 2
  assert s.flags is None
  assert list(s.indices) == [2,4]
  assert s == s
  assert s != p
  assert pickle.loads(pickle.dumps(s)) == s
  assert s.format_summary() == "2"
  #
  p = s
  s = flex.smart_selection(indices=flex.size_t([0,2,4]), all_size=5)
  assert s.all_size == 5
  assert s.selected_size == 3
  assert s.flags.all_eq(flex.bool([True,False,True,False,True]))
  assert list(s.indices) == [0,2,4]
  assert s == s
  assert s != p
  assert pickle.loads(pickle.dumps(s)) == s
  assert s.format_summary() == "3 of 5"
  o = StringIO()
  s.show_summary(out=o)
  assert not show_diff(o.getvalue(), """\
selected elements: 3 of 5
""")
  o = StringIO()
  s.show_summary(out=o, prefix="$&")
  assert not show_diff(o.getvalue(), """\
$&selected elements: 3 of 5
""")
  o = StringIO()
  s.show_summary(out=o, prefix="&$", label="")
  assert not show_diff(o.getvalue(), """\
&$3 of 5
""")
  #
  s = flex.smart_selection(flags=flex.bool(50, flex.size_t([0,2,4])))
  assert s._flags is not None
  assert s._indices is None
  l = pickle.loads(pickle.dumps(s))
  assert l == s
  assert l._flags is None
  assert l._indices is not None
  #
  print("OK")

if (__name__ == "__main__"):
  run(sys.argv[1:])
