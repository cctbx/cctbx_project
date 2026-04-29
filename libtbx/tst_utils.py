from __future__ import absolute_import, division, print_function
from libtbx import utils
from libtbx.test_utils import Exception_expected, approx_equal, show_diff
from six.moves import cStringIO as StringIO
import warnings
import random
import time
import os
import stat
import tempfile
from six.moves import range

def exercise_misc():
  utils.host_and_user().show(prefix="### ")
  time_in_seconds = 1.1
  for i_trial in range(55):
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
    except ValueError as e:
      assert str(e) == 'Error interpreting "%s" as a numeric expression.' % (
        string)
    else: raise Exception_expected
  assert nfs(string="-42") == -42
  assert approx_equal(nfs(string="3.14"), 3.14)
  assert approx_equal(nfs(string="cos(0)"), 1)
  try: nfs(string="xxx(0)")
  except ValueError as e:
    assert str(e).startswith(
      'Error interpreting "xxx(0)" as a numeric expression: ')
  else: raise Exception_expected
  #
  s = "[0.143139, -0.125121, None, -0.308607]"
  assert numstr(values=eval(s)) == s
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
  for size in range(1,5):
    for i1d in range(size):
      assert nd1d(i1d=i1d, sizes=(size,)) == [i1d]
  for sizes in [(1,1), (1,3), (3,1), (2,3)]:
    ni, nj = sizes
    for i in range(ni):
      for j in range(nj):
        i1d = i*nj+j
        assert nd1d(i1d=i1d, sizes=sizes) == [i,j]
  for sizes in [(1,1,1), (1,3,1), (3,2,1), (4,3,2)]:
    ni, nj, nk = sizes
    for i in range(ni):
      for j in range(nj):
        for k in range(nk):
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
  assert approx_equal(list(range(5)) , frange(5))
  assert approx_equal(list(range(5+1)) , samples(5))
  assert approx_equal(list(range(-5)), frange(-5))
  assert approx_equal(list(range(-5-1)), samples(-5))
  assert approx_equal(list(range(1,3)), frange(1, 3))
  assert approx_equal(list(range(1,3+1)), samples(1, 3))
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
  for n in range(4):
    assert len(utils.random_hex_code(number_of_digits=n)) == n
  #
  print("multiprocessing problem:", utils.detect_multiprocessing_problem())
  #
  print("base36_timestamp():", utils.base36_timestamp(), "now")
  print("base36_timestamp():", utils.base36_timestamp(
    seconds_since_epoch=115855*365.2425*24*60*60), "year 115855 CE")
  #
  print("get_svn_revision():", utils.get_svn_revision())
  print("get_build_tag():", utils.get_build_tag())
  # concatenate_python_script
  # XXX the string concatenation here is required to trick libtbx.find_clutter,
  # which will warn about repetition of the future division import.
  script = """
from __future__ """ + """import division
import os.path

def foo():
  print "bar"
"""
  d = tempfile.mkdtemp()
  name = os.path.join(d, "tst_libtbx_utils_python_script.py")
  name2 = os.path.join(d, "tst_libtbx_utils_python_script2.py")
  with open(name, "w") as f:
    f.write(script)
  f = open(name2, "w")
  utils.concatenate_python_script(out=f, file_name=name)
  f.close()
  with open(name2) as f:
    lines = f.readlines()
  have_def = False
  for line in lines :
    assert (not "__future__" in line)
    if line.startswith("def foo"):
      have_def = True
  assert have_def

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
  print("level0", file=level0)
  level0.flush()
  level1 = level0.shift_right()
  print("level1", file=level1)
  level1.flush()
  assert out.getvalue() == ""
  level1.write_buffer()
  assert not show_diff(out.getvalue(), """\
level0
  level1
""")
  print("abc", end='', file=level1)
  level1.write_buffer()
  assert not show_diff(out.getvalue(), """\
level0
  level1
  abc""")
  print(file=level1)
  level1.write_buffer()
  assert not show_diff(out.getvalue(), """\
level0
  level1
  abc
""")
  print("def", end='', file=level1)
  level1.write_buffer()
  assert not show_diff(out.getvalue(), """\
level0
  level1
  abc
  def""")
  level1.write("")
  print("hij", file=level1)
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

def exercise_file_utils():
  dir_name = tempfile.mkdtemp()
  if (not os.path.exists(dir_name)):
    os.mkdir(dir_name)
  sorted_files = []
  for prefix in ["XYZ", "abc", "qwerty", "123"] :
    file_name = os.path.join(dir_name, "%s.txt" % prefix)
    with open(file_name, "w") as f:
      f.write(prefix)
    sorted_files.append(file_name)
    time.sleep(1) # XXX the mtime resolution is in seconds :(
  f = open(os.path.join(dir_name, "hkl.log"), "w")
  f.write("hkl")
  f.close()
  file_names = utils.find_files(dir_name, pattern=".txt$")
  sorted_files_2 = utils.sort_files_by_mtime(file_names)
  assert (sorted_files_2 == sorted_files), '''
  Files not in correct order:
    %s
    %s
  ''' % (sorted_files_2, sorted_files)

def exercise_dir_utils():
  dirs = ["tst_utils_1", "tst_utils_2", "tst_utils_45"]
  for dir_name in dirs :
    if (os.path.isdir(dir_name)) : os.rmdir(dir_name)
  dir_name = utils.create_run_directory("tst_utils")
  assert (os.path.basename(dir_name) == "tst_utils_1")
  dir_name = utils.create_run_directory("tst_utils")
  assert (os.path.basename(dir_name) == "tst_utils_2")
  dir_name = utils.create_run_directory("tst_utils", 45)
  assert (os.path.basename(dir_name) == "tst_utils_45")
  for dir_name in dirs :
    os.rmdir(dir_name)
  file_name = "/cctbx/%s/%s/XXXX.pdb" % (random.random(), random.random())
  try :
    utils.check_if_output_directory_exists(file_name)
  except utils.Sorry :
    pass
  else :
    raise Exception_expected
  dir_name = os.getcwd()
  utils.check_if_output_directory_exists(dir_name=dir_name)
  dir_created = False
  if (not os.path.exists("Dropbox")):
    os.mkdir("Dropbox")
    dir_created = True
  dir_name = os.path.join(os.getcwd(), "Dropbox")
  with warnings.catch_warnings(record=True) as w:
    warnings.simplefilter("always")
    utils.check_if_output_directory_exists(dir_name=dir_name)
    assert len(w) == 1
    assert "Dropbox directory" in str(w[-1].message)
  if (dir_created):
    os.rmdir("Dropbox")
  host_info = utils.host_and_user()
  assert not utils.allow_delete_directory(host_info.homedir)
  target_dir = os.path.join(host_info.homedir, "Downloads")
  assert not utils.allow_delete_directory(target_dir)
  target_dir = os.path.join(host_info.homedir, "data", "lysozyme")
  assert utils.allow_delete_directory(target_dir)

def exercise_retrieve_unless_exists():
  from six.moves import urllib
  filehandle, filename = tempfile.mkstemp(prefix='kings_of_france')
  # we will need to pass filename to functions which will open it
  # on Windows this causes a permission exception
  os.close(filehandle)
  with open(filename, 'w') as f:
    f.write(
      'Henri IV, Louis XIII, Louis XIV, Louis XV, Louis XVI, Louis XVIII')
  digestname = os.path.join(os.path.dirname(f.name), 'digests.txt')
  with open(digestname, 'w') as f:
    f.writelines([
      ('%s %s\n') %
      (os.path.basename(filename), utils.md5_hexdigest(filename)),
      (os.path.basename(filename), utils.sha256_hexdigest(filename)),
      'something_else  yyyyyyy',
    ])
  os.chmod(digestname,
           os.stat(digestname).st_mode | stat.S_IWGRP | stat.S_IWOTH)
  d = tempfile.mkdtemp()
  targetname = os.path.join(d, 'target')
  try: os.remove(targetname)
  except Exception: pass
  url = 'file:' + urllib.request.pathname2url(filename)
  assert (utils.retrieve_unless_exists(url=url, filename=targetname) ==
          "Downloaded")
  with open(filename) as source, open(targetname) as target:
    assert source.read() == target.read()
  assert (utils.retrieve_unless_exists(url=url, filename=targetname) ==
          "Cached")
  with open(filename) as source, open(targetname) as target:
    assert source.read() == target.read()

def exercise_str_unicode():
  # tests for to_unicode and to_str
  s = '\xc3\x85'
  u = u'\xc5'
  assert(to_unicode(s) == u)
  assert(to_str(u) == s)

def exercise_group_args():
  from libtbx import group_args
  out = StringIO()
  a = group_args(
      a=1,
      b=2,
      c=3)
  assert a.a==1
  assert a.b==2
  assert a.c==3
  b = group_args(
      d = 'd',
      e = 'e')
  assert b.d=='d'
  assert b.e=='e'
  print(a, file=out)
  v = out.getvalue()
  assert not show_diff(v, """group_args
  a                              : 1
  b                              : 2
  c                              : 3\n""")
  a.merge(b)
  assert a.a==1
  assert a.b==2
  assert a.c==3
  assert a.d=='d'
  assert a.e=='e'
  assert b.d=='d'
  assert b.e=='e'
  c = group_args(
      a = 11,
      b = 12)
  a.merge(c)
  assert a.a==11
  assert a.b==12
  assert a.c==3
  assert c.a==11
  assert c.b==12
  #
  r = group_args(x=1)
  r.stop_dynamic_attributes()
  err = None
  try:
    r.w = 0
  except TypeError as e:
    err = e
  assert str(err) == "Dynamic attributes disabled."

def exercise_round2():
  assert(2 == int(utils.round2(1.5, 0)))
  assert(3 == int(utils.round2(2.5, 0)))
  assert(-2 == int(utils.round2(-1.5, 0)))
  assert(-3 == int(utils.round2(-2.5, 0)))
  assert approx_equal(0.2, utils.round2(0.15, 1))
  assert approx_equal(0.3, utils.round2(0.25, 1))
  assert approx_equal(-0.2, utils.round2(-0.15, 1))
  assert approx_equal(-0.3, utils.round2(-0.25, 1))

def exercise_guess_total_memory():
  assert(utils.guess_total_memory() > 0)

def exercise_display_context():
  text = """
   line with word1
   another line
   another line with word2
   another line with word1
   line with word3
"""
  from libtbx.utils import display_context
  text_block_list = display_context(text = text,
     n_context = 1, search_word = 'word1', quiet = True)
  assert [text_block_list[0].text_block] == [ '     \n  **    line with word1\n        another line\n']

  text_block_list = display_context(text = text,
     n_context = 1, search_word = 'word2', quiet = True)
  assert [text_block_list[0].text_block] == [ '        another line\n  **    another line with word2\n        another line with word1\n']

  text_block_list = display_context(text = text,
     n_context = 1, search_word = 'word2', required_word ='word1', quiet = True)
  assert [text_block_list[0].text_block] == [ '        another line\n  **    another line with word2\n        another line with word1\n']


def run(args):
  assert len(args) == 0
  if '--exercise-retrieve-unless-exists' in args:
    exercise_retrieve_unless_exists()
  else:
    print('Skipping exercise_retrieve_unless_exists')
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
  exercise_file_utils()
  exercise_dir_utils()
  exercise_group_args()
  exercise_round2()
  exercise_display_context()
  print(utils.format_cpu_times())

if (__name__ == "__main__"):
  import sys
  run(args=sys.argv[1:])
