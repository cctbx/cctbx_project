from libtbx import utils
from libtbx.test_utils import Exception_expected, approx_equal, show_diff
from cStringIO import StringIO

def exercise_forward_compatibility():
  import itertools
  f = itertools.izip_longest
  assert list(f([], [])) == []
  assert list(f([1,2], [3])) == [(1, 3), (2, None)]
  assert list(f([1], [2,3])) == [(1, 2), (None, 3)]

def exercise_misc():
  utils.host_and_user().show(prefix="### ")
  time_in_seconds = 1.1
  for i_trial in xrange(55):
    time_in_seconds = time_in_seconds**1.1
    time_units, time_unit = utils.human_readable_time(
      time_in_seconds=time_in_seconds)
    assert approx_equal(
      utils.human_readable_time_as_seconds(time_units, time_unit),
      time_in_seconds)
  #
  fts = utils.format_timestamp
  f12 = utils.format_timestamp_12_hour
  f24 = utils.format_timestamp_24_hour
  def check(string, expected):
    assert len(string) == len(expected)
  check(f12(1280007000), 'Jul 24 2010 02:30 PM')
  check(f24(1280007000), 'Jul 24 2010 14:30')
  check(f12(1280007000, True), '24-07-10 02:30 PM')
  check(f24(1280007000, True), '24-07-10 14:30')
  check(fts(1280007000), 'Jul 24 2010 02:30 PM')
  #
  nfs = utils.number_from_string
  for string in ["True", "False"]:
    try: nfs(string=string)
    except ValueError, e:
      assert str(e) == 'Error interpreting "%s" as a numeric expression.' % (
        string)
    else: raise Exception_expected
  assert nfs(string="-42") == -42
  assert approx_equal(nfs(string="3.14"), 3.14)
  assert approx_equal(nfs(string="cos(0)"), 1)
  try: nfs(string="xxx(0)")
  except ValueError, e:
    assert str(e).startswith(
      'Error interpreting "xxx(0)" as a numeric expression: ')
  else: raise Exception_expected
  #
  s = "[0.143139, -0.125121, None, -0.308607]"
  assert numstr(values=eval(s)) == s
  s = "[0.1431391, -0.1251212, None, -0.3086073]"
  assert numstr7(values=eval(s)) == s
  #
  for s,i in {"2000000" : 2000000,
              "2k" : 2048,
              "2Kb" : 2048,
              "2 Kb" : 2048,
              "5Mb" : 5*1024*1024,
              "2.5Gb" : 2.5*1024*1024*1024,
              "1T": 1024*1024*1024*1024,
              10000 : 10000,
              5.5 : 5.5,
              }.items():
    assert utils.get_memory_from_string(s) == i
  #
  assert utils.tupleize(1) == (1,)
  assert utils.tupleize("abcde") == ('a', 'b', 'c', 'd', 'e')
  assert utils.tupleize([1,2,3]) == (1,2,3)
  #
  sf = utils.search_for
  assert sf(pattern="fox", mode="==", lines=["fox", "foxes"]) \
      == ["fox"]
  assert sf(pattern="o", mode="find", lines=["fox", "bird", "mouse"]) \
      == ["fox", "mouse"]
  assert sf(pattern="fox", mode="startswith", lines=["fox", "foxes"]) \
      == ["fox", "foxes"]
  assert sf(pattern="xes", mode="endswith", lines=["fox", "foxes"]) \
      == ["foxes"]
  assert sf(pattern="es$", mode="re.search", lines=["geese", "foxes"]) \
      == ["foxes"]
  assert sf(pattern="ge", mode="re.match", lines=["geese", "angel"]) \
      == ["geese"]
  #
  nd1d = utils.n_dim_index_from_one_dim
  for size in xrange(1,5):
    for i1d in xrange(size):
      assert nd1d(i1d=i1d, sizes=(size,)) == [i1d]
  for sizes in [(1,1), (1,3), (3,1), (2,3)]:
    ni, nj = sizes
    for i in xrange(ni):
      for j in xrange(nj):
        i1d = i*nj+j
        assert nd1d(i1d=i1d, sizes=sizes) == [i,j]
  for sizes in [(1,1,1), (1,3,1), (3,2,1), (4,3,2)]:
    ni, nj, nk = sizes
    for i in xrange(ni):
      for j in xrange(nj):
        for k in xrange(nk):
          i1d = (i*nj+j)*nk+k
          assert nd1d(i1d=i1d, sizes=sizes) == [i,j,k]
  #
  from libtbx import easy_run
  b = easy_run.fully_buffered(
    command="libtbx.raise_exception_for_testing")
  for lines in [b.stdout_lines, b.stderr_lines]:
    assert lines[0].startswith("EXCEPTION_INFO: show_stack(0): ")
    assert lines[-1] == "EXCEPTION_INFO: RuntimeError: Just for testing."
  b = easy_run.fully_buffered(
    command="libtbx.raise_exception_for_testing silent")
  b.raise_if_errors_or_output()
  #
  frange = utils.frange
  samples = utils.samples
  assert approx_equal([i/10. for i in range(-2,2)], frange(-0.2,0.2,0.1))
  assert approx_equal([i/10. for i in range(-2,2+1)], samples(-0.2,0.2,0.1))
  assert approx_equal([i/10. for i in range(2,-2,-1)], frange(0.2,-0.2,-0.1))
  assert approx_equal([i/10. for i in range(2,-2-1,-1)], samples(0.2,-0.2,-0.1))
  assert approx_equal([i/4. for i in range(4,8)], frange(1, 2, 0.25))
  assert approx_equal([i/4. for i in range(4,8+1)], samples(1, 2, 0.25))
  assert approx_equal([0.2+i/3. for i in range(4)], frange(0.2, 1.3, 1./3))
  assert approx_equal([0.2+i/3. for i in range(4)], samples(0.2, 1.3, 1./3))
  assert approx_equal(range(5) , frange(5))
  assert approx_equal(range(5+1) , samples(5))
  assert approx_equal(range(-5), frange(-5))
  assert approx_equal(range(-5-1), samples(-5))
  assert approx_equal(range(1,3), frange(1, 3))
  assert approx_equal(range(1,3+1), samples(1, 3))
  assert approx_equal([i/10. for i in range(20,9,-2)], frange(2.0,0.9,-0.2))
  assert approx_equal([i/10. for i in range(20,9,-2)], samples(2.0,0.9,-0.2))
  #
  ff = utils.format_float_with_standard_uncertainty
  assert ff(21.234567, 0.0013) == "21.2346(13)"
  assert ff(21.234567, 0.0023) == "21.235(2)"
  assert ff(12345, 45) == "12350(50)"
  assert ff(12.3,1.2) == "12.3(12)"
  assert ff(-0.2451, 0.8135) == "-0.2(8)"
  assert ff(1.234, 0.196) == "1.2(2)"
  assert ff(1.234, 0.193) == "1.23(19)"
  #
  for n in xrange(4):
    assert len(utils.random_hex_code(number_of_digits=n)) == n
  #
  print "multiprocessing problem:", utils.detect_multiprocessing_problem()
  #
  print "base36_timestamp():", utils.base36_timestamp(), "now"
  print "base36_timestamp():", utils.base36_timestamp(
    seconds_since_epoch=115855*365.2425*24*60*60), "year 115855 CE"
  #
  print "get_svn_revision():", utils.get_svn_revision()
  print "get_build_tag():", utils.get_build_tag()

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

def run(args):
  assert len(args) == 0
  exercise_forward_compatibility()
  exercise_misc()
  assert utils.sequence_index_dict(["a", "b"]) == {"a": 0, "b": 1}
  assert utils.flat_list(0) == [0]
  assert utils.flat_list([1,2,3]) == [1,2,3]
  assert utils.flat_list([1,[2,3,4],3]) == [1,2,3,4,3]
  assert utils.flat_list([1,[2,3,4],[[3,4],[5,6]]]) == [1,2,3,4,3,4,5,6]
  try:
    raise RuntimeError("Trial")
  except KeyboardInterrupt: raise
  except Exception:
    assert utils.format_exception() == "RuntimeError: Trial"
  else: raise Exception_expected
  try:
    assert 1 == 2
  except KeyboardInterrupt: raise
  except Exception:
    s = utils.format_exception()
    assert s.startswith("AssertionError: ")
    assert s.find("tst_utils.py line ") >= 0
  else: raise Exception_expected
  exercise_indented_display()
  exercise_approx_equal()
  print utils.format_cpu_times()

if (__name__ == "__main__"):
  import sys
  run(args=sys.argv[1:])
