from libtbx import utils
from libtbx.test_utils import show_diff
from cStringIO import StringIO

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

def exercise():
  assert utils.flat_list(0) == [0]
  assert utils.flat_list([1,2,3]) == [1,2,3]
  assert utils.flat_list([1,[2,3,4],3]) == [1,2,3,4,3]
  assert utils.flat_list([1,[2,3,4],[[3,4],[5,6]]]) == [1,2,3,4,3,4,5,6]
  try:
    raise RuntimeError("Trial")
  except KeyboardInterrupt: raise
  except:
    assert utils.format_exception() == "RuntimeError: Trial"
  exercise_indented_display()
  print "OK"

if (__name__ == "__main__"):
  exercise()
