def exercise_posix_relpath(f, enable_abspath_if_through_root):
  # based on .test_relpath() in Python-2.7.2/Lib/test/test_posixpath.py
  import os
  (real_getcwd, os.getcwd) = (os.getcwd, lambda: r"/home/user/bar")
  try:
    curdir = os.path.split(os.getcwd())[-1]
    def check(args, r_expected, ra_expected=None):
      if (ra_expected is None): ra_expected = r_expected
      result = f(*args)
      assert result == r_expected
      if (enable_abspath_if_through_root):
        result = f(*(args), **{"enable_abspath_if_through_root": True})
        assert result == ra_expected
    check(["a"], "a")
    check([os.path.abspath("a")], "a")
    check(["a/b"], "a/b")
    check(["../a/b"], "../a/b")
    check(["a", "../b"], "../"+curdir+"/a")
    check(["a/b", "../c"], "../"+curdir+"/a/b")
    check(["a", "b/c"], "../../a")
    check(["a", "a"], ".")
    check(["/foo/bar/bat", "/x/y/z"], '../../../foo/bar/bat', "/foo/bar/bat")
    check(["/foo/bar/bat", "/foo/bar"], 'bat')
    check(["/foo/bar/bat", "/"], 'foo/bar/bat', "/foo/bar/bat")
    check(["/", "/foo/bar/bat"], '../../..', "/")
    check(["/foo/bar/bat", "/x"], '../foo/bar/bat', "/foo/bar/bat")
    check(["/x", "/foo/bar/bat"], '../../../x', "/x")
    check(["/", "/"], '.', "/")
    check(["/a", "/a"], '.')
    check(["/a/b", "/a/b"], '.')
    #
    # added tests, may be partially redundant
    check(["/a/b/c", "/"], "a/b/c", "/a/b/c")
    check(["/a/b/c", "/x"], "../a/b/c", "/a/b/c")
    check(["/a/b/c", "/x/y"], "../../a/b/c", "/a/b/c")
    check(["/a/b/c", "/a"], "b/c")
    check(["/a/b/c", "/a/b"], "c")
    check(["/a/b/c", "/a/b/c"], ".")
    check(["/a/b/c", "/a/b/c/d"], "..")
    check(["/conky/mountpoint/a", "/conky/mountpoint/b/c"], "../../a")
  finally:
    os.getcwd = real_getcwd

def exercise_nt_relpath(f, enable_abspath_if_through_root):
  import os
  currentdir = os.path.split(os.getcwd())[-1]
  def check(args, r_expected, ra_expected=None):
    if (ra_expected is None): ra_expected = r_expected
    result = f(*args)
    assert result == r_expected
    if (enable_abspath_if_through_root):
      result = f(*(args), **{"enable_abspath_if_through_root": True})
      assert result == ra_expected
  check(["a"], 'a')
  check([os.path.abspath("a")], 'a')
  check(["a/b"], 'a\\b')
  check(["../a/b"], '..\\a\\b')
  check(["a", "../b"], '..\\'+currentdir+'\\a')
  check(["a/b", "../c"], '..\\'+currentdir+'\\a\\b')
  check(["a", "b/c"], '..\\..\\a')
  check(["//conky/mountpoint/a", "//conky/mountpoint/b/c"], '..\\..\\a')
  check(["a", "a"], '.')
  check(["/foo/bar/bat", "/x/y/z"], '..\\..\\..\\foo\\bar\\bat',
    "\\foo\\bar\\bat")
  check(["/foo/bar/bat", "/foo/bar"], 'bat')
  check(["/foo/bar/bat", "/"], 'foo\\bar\\bat', "\\foo\\bar\\bat")
  check(["/", "/foo/bar/bat"], '..\\..\\..', "\\")
  check(["/foo/bar/bat", "/x"], '..\\foo\\bar\\bat', "\\foo\\bar\\bat")
  check(["/x", "/foo/bar/bat"], '..\\..\\..\\x', "\\x")
  check(["/", "/"], '.', "\\")
  check(["/a", "/a"], '.')
  check(["/a/b", "/a/b"], '.')
  check(["c:/foo", "C:/FOO"], '.')
  #
  if (enable_abspath_if_through_root):
    assert f("c:\\foo", "d:\\foo", True) == "c:\\foo"
    assert f("//m/d", "//n/d", True) == "\\\\m\\d"
    assert f("d:\\foo", "//n/d", True) == "d:\\foo"

def exercise_relpath():
  import sys, os
  if (os.name == "nt"):
    exercise = exercise_nt_relpath
  else:
    exercise = exercise_posix_relpath
  if (sys.version_info[:2] >= (2,7)): # Python 2.6 is slightly different
    from os.path import relpath
    exercise(relpath, False)
  from libtbx.path import relpath
  exercise(relpath, True)

def exercise_move_old_create_new_directory():
  from libtbx.path import move_old_create_new_directory as mocnd
  import os
  mocnd("tmp_mocnd")
  assert len(os.listdir("tmp_mocnd")) == 0
  for i in xrange(3):
    mocnd("tmp_mocnd/a")
  assert sorted(os.listdir("tmp_mocnd")) == ["a", "a_001", "a_002"]
  open("tmp_mocnd/a_23", "w")
  mocnd("tmp_mocnd/a")
  assert sorted(os.listdir("tmp_mocnd")) == [
    "a", "a_001", "a_002", "a_024", "a_23"]
  open("tmp_mocnd/a_log", "w")
  mocnd("tmp_mocnd/a")
  assert sorted(os.listdir("tmp_mocnd")) == [
    "a", "a_001", "a_002", "a_024", "a_025", "a_23", "a_log"]
  mocnd("tmp_mocnd/b")
  assert sorted(os.listdir("tmp_mocnd")) == [
    "a", "a_001", "a_002", "a_024", "a_025", "a_23", "a_log", "b"]
  mocnd("tmp_mocnd/b", serial_sep="", serial_fmt="%d")
  assert sorted(os.listdir("tmp_mocnd")) == [
    "a", "a_001", "a_002", "a_024", "a_025", "a_23", "a_log", "b", "b1"]

def run(args):
  assert len(args) == 0
  exercise_relpath()
  exercise_move_old_create_new_directory()
  from libtbx.path import random_new_directory_name
  assert len(random_new_directory_name()) == len("tmp_dir_00000000")
  print "OK"

if (__name__ == "__main__"):
  import sys
  run(args=sys.argv[1:])
