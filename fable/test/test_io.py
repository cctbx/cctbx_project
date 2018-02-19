from __future__ import absolute_import, division, print_function

import os

from libtbx import easy_run
import pytest

def remove_file(path):
  if os.path.exists(path):
    os.remove(path)
  assert not os.path.exists(path)

def build_cmds(tst_f):
  from fable import simple_compilation
  comp_env = simple_compilation.environment()
  exe = comp_env.exe_suffix
  result = []
  remove_file("fable_cout")
  cmd = "fable.cout %s --link" % tst_f
  easy_run.fully_buffered(command=cmd).raise_if_errors()
  assert os.path.exists("fable_cout"+exe)
  result.append("fable_cout"+exe)
  return [os.path.join(".", cmd) for cmd in result]

@pytest.mark.parametrize("status, expected", [
    ("old", [['F'], False, ['T'], True, ['T'], True, 1]),
    ("new", [['T'], True, ['F'], True, ['F'], True, 1]),
    ("unknown", [['T'], True, ['T'], True, ['T'], True, 1]),
    ("scratch", [['T'], False, ['T'], True, ['T'], True, 1]),
])
def test_open(tmpdir, status, expected):
  tmpdir.chdir()
  testfile = tmpdir.join("exercise_open_%s.f" % status)
  testfile.write("""\
      program prog
      open(1, file='exercise_open.tmp', status='%s', iostat=ios)
      write(6, '(l1)') (ios .eq. 0)
      end
""" % status)

  for cmd in build_cmds(tst_f=testfile.strpath):
    remove_file("exercise_open.tmp")
    stdout_1 = easy_run.fully_buffered(
      command=cmd).raise_if_errors().stdout_lines
    exists_1 = os.path.exists("exercise_open.tmp")

    open("exercise_open.tmp", "w")
    assert os.path.exists("exercise_open.tmp")
    stdout_2 = easy_run.fully_buffered(
      command=cmd).raise_if_errors().stdout_lines
    exists_2 = os.path.exists("exercise_open.tmp")

    open("exercise_open.tmp", "w").write("X")
    stdout_3 = easy_run.fully_buffered(
      command=cmd).raise_if_errors().stdout_lines
    exists_3 = os.path.exists("exercise_open.tmp")
    if (exists_3):
      size_3 = os.path.getsize("exercise_open.tmp")
    else:
      size_3 = None

    results = [
      stdout_1, exists_1, stdout_2, exists_2, stdout_3, exists_3, size_3]
    assert results == expected

def test_mixed_read_write(tmpdir):
  tmpdir.chdir()
  tmp = "exercise_mixed_read_write.tmp"
  tst_f = "exercise_mixed_read_write.f"
  open(tst_f, "w").write("""\
      program prog
      open(
     &  unit=1,
     &  file='%s',
     &  status='old')
      read(1, '(i2)') num
      write(6, '(i2)') num*2
      write(1, '(i2)') 78
      end
""" % tmp)
  for cmd in build_cmds(tst_f=tst_f):
    open(tmp, "w").write("""\
12
34
56
""")
    stdout = easy_run.fully_buffered(
      command=cmd).raise_if_errors().stdout_lines
    assert stdout == ["24"]
    tmp_text = open(tmp, "rb").read()
    assert tmp_text == "12\n78\n".replace("\n", os.linesep)

def test_read_from_non_existing_file(tmpdir):
  tmpdir.chdir()
  tst_f = "exercise_read_from_non_existing_file.f"
  remove_file("fem_io_unit_001")
  open(tst_f, "w").write("""\
      program prog
      read(1, *) num
      end
""")
  for cmd in build_cmds(tst_f=tst_f):
    stdout = easy_run.fully_buffered(
      command=cmd, join_stdout_stderr=True).stdout_lines
    assert stdout == ["std::exception what(): End of input during read"]
