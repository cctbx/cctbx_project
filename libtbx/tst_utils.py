from libtbx import utils
from libtbx.test_utils import Exception_expected, approx_equal, show_diff
from cStringIO import StringIO

def exercise_forward_compatibility():
  import itertools
  f = itertools.izip_longest
  assert list(f([], [])) == []
  assert list(f([1,2], [3])) == [(1, 3), (2, None)]
  assert list(f([1], [2,3])) == [(1, 2), (None, 3)]

def exercise_user_plus_sys_time():
  s = StringIO()
  utils.user_plus_sys_time().show_elapsed(out=s, prefix="e: ")
  s = s.getvalue()
  assert s.startswith("e: ")
  assert s.endswith(" s")
  utils.user_plus_sys_time().show_delta(out=s, prefix="d: ")
  s = s.getvalue()
  assert s.startswith("d: ")
  assert s.endswith(" s")

def exercise_indented_display():
  out = StringIO()
  level0 = utils.buffered_indentor(file_object=out)
  print >> level0, "level0"
  level0.flush()
  level1 = level0.shift_right()
  print >> level1, "level1"
  level1.flush()
  assert out.getvalue() == ""
  level1.write_buffer()
  assert not show_diff(out.getvalue(), """\
level0
  level1
""")
  print >> level1, "abc",
  level1.write_buffer()
  assert not show_diff(out.getvalue(), """\
level0
  level1
  abc""")
  print >> level1
  level1.write_buffer()
  assert not show_diff(out.getvalue(), """\
level0
  level1
  abc
""")
  print >> level1, "def",
  level1.write_buffer()
  assert not show_diff(out.getvalue(), """\
level0
  level1
  abc
  def""")
  level1.write("")
  print >> level1, "hij"
  level1.write_buffer()
  assert not show_diff(out.getvalue(), """\
level0
  level1
  abc
  def hij
""")

def exercise_approx_equal():
  assert approx_equal(1., 1. + 1e-11)
  assert approx_equal(1+1j, 0.997+1.004j, eps=1e-2)
  assert approx_equal(1, 0.997+0.004j, eps=1e-2)
  assert approx_equal(1+0.003j, 0.997, eps=1e-2)
  assert approx_equal([ 2.5, 3.4+5.8j, 7.89],
                      [ 2.4+0.1j, 3.5+5.9j, 7.90], eps=0.2)

def exercise():
  exercise_forward_compatibility()
  assert utils.sequence_index_dict(["a", "b"]) == {"a": 0, "b": 1}
  assert utils.flat_list(0) == [0]
  assert utils.flat_list([1,2,3]) == [1,2,3]
  assert utils.flat_list([1,[2,3,4],3]) == [1,2,3,4,3]
  assert utils.flat_list([1,[2,3,4],[[3,4],[5,6]]]) == [1,2,3,4,3,4,5,6]
  try:
    raise RuntimeError("Trial")
  except KeyboardInterrupt: raise
  except:
    assert utils.format_exception() == "RuntimeError: Trial"
  else: raise Exception_expected
  try:
    assert 1 == 2
  except KeyboardInterrupt: raise
  except:
    s = utils.format_exception()
    assert s.startswith("AssertionError: ")
    assert s.find("tst_utils.py line ") >= 0
  else: raise Exception_expected
  exercise_indented_display()
  exercise_approx_equal()
  print "OK"

if (__name__ == "__main__"):
  exercise()
